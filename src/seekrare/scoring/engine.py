"""
Dual-Dynamic Scoring Engine

Core innovation of SeekRare.

For each variant v and annotation column c:

    score_c(v) = normalized(v[c]) × semantic_similarity(c["hpo"], relevant_HPOs) × weight[c]

The final variant score:
    score(v) = Σ_c  score_c(v)

Where:
    - weight[c] is from LLM's weight vector (changes per patient)
    - semantic_similarity is from HPO ontology (changes per HPO relevance)

This is "dual-dynamic" because BOTH dimensions are personalized.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from loguru import logger


# ── Column normalization strategies ──────────────────────────────────────────

def _normalize_minmax(series: pd.Series) -> pd.Series:
    """Min-max normalization to [0, 1]."""
    mn, mx = series.min(), series.max()
    if mx == mn:
        return pd.Series(0.0, index=series.index)
    return (series - mn) / (mx - mn)


def _normalize_rank(series: pd.Series) -> pd.Series:
    """Rank-based normalization: highest value gets 1.0."""
    return series.rank(pct=True, ascending=True)


def _normalize_clinvar(series: pd.Series) -> pd.Series:
    """
    ClinVar-specific normalization.
    Pathogenic > Likely pathogenic > VUS > Likely benign > Benign
    """
    mapping = {
        "Pathogenic": 1.0,
        "Likely pathogenic": 0.8,
        "Uncertain significance": 0.5,
        "Likely benign": 0.2,
        "Benign": 0.0,
        "Conflicting interpretations": 0.4,
        "not provided": 0.3,
        "Pathogenic/Likely Pathogenic": 0.9,
    }
    return series.map(mapping).fillna(0.3)


def _normalize_cadd(series: pd.Series) -> pd.Series:
    """CADD: already [0, 99], higher = more deleterious."""
    return _normalize_minmax(series.fillna(0))


def _normalize_gnomad(series: pd.Series) -> pd.Series:
    """gnomAD: invert (lower AF = higher score) and normalize."""
    inv = (1.0 - series.fillna(0).clip(upper=1.0)) ** 2
    return _normalize_minmax(inv)


def _normalize_impact(series: pd.Series) -> pd.Series:
    """Impact: HIGH > MODERATE > LOW."""
    mapping = {"HIGH": 1.0, "MODERATE": 0.6, "LOW": 0.3, "MODIFIER": 0.1}
    return series.map(mapping).fillna(0.1)


def _normalize_sift(series: pd.Series) -> pd.Series:
    """SIFT: 0=damaging, 1=tolerated → invert."""
    return 1.0 - _normalize_minmax(series.fillna(0.5))


def _normalize_polyphen(series: pd.Series) -> pd.Series:
    """PolyPhen: 0=benign, 1=damaging."""
    return _normalize_minmax(series.fillna(0.5))


def _normalize_stars(series: pd.Series) -> pd.Series:
    """ClinVar stars: 0-4 → normalize."""
    return _normalize_minmax(series.fillna(0))


# ── HPO semantic similarity ─────────────────────────────────────────────────

def compute_hpo_similarity(
    hpo_terms: str | list[str],
    relevant_hpos: list[dict],
    column_name: str = "hpo_terms",
) -> float:
    """
    Compute semantic similarity between a variant's HPO annotations
    and the patient's relevant HPO terms (from LLM).

    If variant has no HPO terms → 0.0 (irrelevant to patient phenotype)
    If column is not HPO-related → 1.0 (weight applies directly)

    Parameters
    ----------
    hpo_terms : str or list
        Variant's HPO annotations (pipe-separated string or list)
    relevant_hpos : list[dict]
        From LLM output: [{"hpo_id": "HP:0001250", "score": 0.95}, ...]

    Returns
    -------
    float
        Semantic similarity in [0.0, 1.0]
    """
    if not relevant_hpos:
        return 1.0  # No preference, no boost

    if not hpo_terms or (isinstance(hpo_terms, str) and hpo_terms.strip() in ("", ".")):
        return 0.0  # Variant has no HPO terms, not relevant

    if isinstance(hpo_terms, str):
        variant_hpos = set(t.strip() for t in hpo_terms.replace(";", "|").split("|") if t.strip())
    else:
        variant_hpos = set(str(t).strip() for t in hpo_terms if t)

    # Build relevance lookup: hpo_id → score
    hpo_relevance = {h["hpo_id"]: h["score"] for h in relevant_hpos}

    # If no overlap at all → 0
    overlap = variant_hpos & set(hpo_relevance.keys())
    if not overlap:
        return 0.0

    # Weighted average of overlapping HPO relevance scores
    total_score = sum(hpo_relevance[hpo_id] for hpo_id in overlap)
    max_possible = sum(hpo_relevance.values())

    return total_score / max_possible if max_possible > 0 else 0.0


# ── Main scoring engine ──────────────────────────────────────────────────────

# Annotation column → (normalizer function, applies_to_hpo_col)
_ANNOTATION_SPECS: dict[str, tuple[callable, bool]] = {
    "clinvar_significance": (_normalize_clinvar, False),
    "clinvar_stars": (_normalize_stars, False),
    "cadd_score": (_normalize_cadd, False),
    "gnomad_af": (_normalize_gnomad, False),
    "sift_score": (_normalize_sift, False),
    "polyphen_score": (_normalize_polyphen, False),
    "impact": (_normalize_impact, False),
    "pLI": (_normalize_minmax, False),
    "oeLoF": (_normalize_minmax, False),
    "hpo_terms": (None, True),  # special: handled by HPO similarity
}


class DualDynamicScorer:
    """
    Dual-Dynamic Scoring Engine.

    Applies LLM-derived weights and HPO semantic similarity to each
    variant's annotation columns, producing a personalized score.

    Parameters
    ----------
    weight_vector : dict[str, float]
        From LLM: e.g., {"clinvar_significance": 0.35, "hpo_terms": 0.30, ...}
        Weights should sum to 1.0 (will be renormalized if not).
    """

    def __init__(self, weight_vector: dict[str, float]):
        # Renormalize weights
        total = sum(weight_vector.values())
        if total <= 0:
            logger.warning("Weight vector sums to 0 or negative, using equal weights")
            n = len(weight_vector) or 1
            self.weights = {k: 1.0 / n for k in weight_vector}
        else:
            self.weights = {k: v / total for k, v in weight_vector.items()}

        logger.info(f"DualDynamicScorer initialized with weights: {self.weights}")

    def score(
        self,
        df: pd.DataFrame,
        relevant_hpos: list[dict],
    ) -> pd.DataFrame:
        """
        Score all variants with dual-dynamic scoring.

        Parameters
        ----------
        df : pd.DataFrame
            Annotated variants (from annotate_variants)
        relevant_hpos : list[dict]
            From LLM: [{"hpo_id": "...", "score": 0.95}, ...]

        Returns
        -------
        pd.DataFrame
            Original df with new columns:
            - score_{col} for each annotation column
            - hpo_similarity: semantic similarity to relevant HPOs
            - seekrare_score: final weighted sum
        """
        df = df.copy()

        # Compute per-column normalized scores and HPO similarity
        column_scores = {}
        hpo_similarities = np.zeros(len(df))

        for col, (normalizer, is_hpo) in _ANNOTATION_SPECS.items():
            if col not in df.columns:
                continue

            if is_hpo:
                # HPO similarity for each variant
                for idx, row in df.iterrows():
                    hpo_similarities[idx] = compute_hpo_similarity(
                        row.get("hpo_terms", ""), relevant_hpos
                    )
                column_scores[col] = hpo_similarities.copy()
            else:
                try:
                    col_data = pd.to_numeric(df[col], errors="coerce")
                    normalized = normalizer(col_data)
                    column_scores[col] = normalized.values
                except Exception as e:
                    logger.warning(f"Failed to normalize column '{col}': {e}")
                    continue

        # ── Compute final seekrare_score ──────────────────────────
        final_scores = np.zeros(len(df))
        used_weights = []

        for col, spec in _ANNOTATION_SPECS.items():
            if col not in column_scores:
                continue
            weight = self.weights.get(col, 0.0)
            if weight <= 0:
                continue

            normalized_scores = column_scores[col]

            # HPO similarity applies as an extra multiplier for hpo column
            if spec[1]:  # is_hpo
                # hpo_similarity already in column_scores
                adjusted = normalized_scores * hpo_similarities
            else:
                # Non-HPO columns: weighted by column weight × HPO similarity of the variant
                # (variants with no relevant HPO get lower scores)
                avg_hpo_sim = np.mean(hpo_similarities) if np.sum(hpo_similarities) > 0 else 1.0
                adjusted = normalized_scores * (0.5 + 0.5 * hpo_similarities)

            df[f"score_{col}"] = adjusted
            final_scores += weight * adjusted
            used_weights.append((col, weight))

        df["seekrare_score"] = final_scores
        df["hpo_similarity"] = hpo_similarities

        logger.info(
            f"Scoring complete. Used columns: {[w[0] for w in used_weights]}. "
            f"Score range: [{final_scores.min():.3f}, {final_scores.max():.3f}]"
        )

        return df
