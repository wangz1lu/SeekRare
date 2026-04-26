"""Rank variants by seekrare_score and export."""

from __future__ import annotations

import pandas as pd


def rank_variants(df: pd.DataFrame, top_k: int = 50) -> pd.DataFrame:
    """
    Sort variants by seekrare_score (descending) and return top_k.

    Parameters
    ----------
    df : pd.DataFrame
        Scored variants DataFrame
    top_k : int
        Number of top candidates to return

    Returns
    -------
    pd.DataFrame
        Sorted, top-k variants
    """
    if "seekrare_score" not in df.columns:
        raise ValueError("seekrare_score column not found. Did you run scoring?")

    ranked = df.sort_values("seekrare_score", ascending=False).head(top_k)
    ranked = ranked.reset_index(drop=True)
    ranked["rank"] = range(1, len(ranked) + 1)

    # Reorder columns: rank first, score second
    cols = ["rank", "seekrare_score", "hpo_similarity"] + [c for c in ranked.columns if c not in ("rank", "seekrare_score", "hpo_similarity")]
    return ranked[[c for c in cols if c in ranked.columns]]
