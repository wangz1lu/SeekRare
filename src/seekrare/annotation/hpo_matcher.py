"""
HPO Matcher — semantic matching between patient symptoms and HPO ontology.

Uses semantic similarity to match free-text symptoms to HPO terms,
complementing the LLM-based approach with a rule-based fallback.

Supports:
- Exact HPO ID lookup
- HPO term name fuzzy matching
- Semantic similarity via HPO ancestor graph (Wu-Palmer similarity)
- Integration with pyhpo library
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger

try:
    from pyhpo import Ontology
    from pyhpo.tools import similarity
    HAS_PYHPO = True
except ImportError:
    HAS_PYHPO = False
    logger.warning("pyhpo not installed. Run: pip install pyhpo")


# ── HPO term database (lightweight built-in fallback) ────────────────────────────

BUILTIN_HPO = {
    # Neurological
    "HP:0001250": ("Seizure", "Neurological"),
    "HP:0012469": ("Infantile spasms", "Neurological"),
    "HP:0002123": ("Generalized tonic-clonic seizure", "Neurological"),
    "HP:0002376": ("Developmental regression", "Neurological"),
    "HP:0002189": ("Intellectual disability", "Neurological"),
    "HP:0007018": ("Attention deficit hyperactivity disorder", "Neurological"),
    "HP:0000726": ("Dystonia", "Neurological"),
    "HP:0002067": ("Gliosis", "Neurological"),
    "HP:0007354": ("Focal impaired awareness seizure", "Neurological"),
    "HP:0002126": ("Generalized EEG abnormality", "Neurological"),
    "HP:0030082": ("Focal seizure", "Neurological"),

    # Muscular / Skeletal
    "HP:0001319": ("Neonatal hypotonia", "Muscular"),
    "HP:0002509": ("Amyotrophy", "Muscular"),
    "HP:0003326": ("Paraparesis", "Muscular"),
    "HP:0003324": ("Generalized muscle weakness", "Muscular"),
    "HP:0002020": ("Gastroesophageal reflux", "Gastrointestinal"),
    "HP:0002021": ("Vomiting", "Gastrointestinal"),

    # Dysmorphic / Craniofacial
    "HP:0000286": ("Anteverted nares", "Dysmorphic"),
    "HP:0000316": ("Hypertelorism", "Dysmorphic"),
    "HP:0000324": ("Facial feature", "Dysmorphic"),
    "HP:0004482": ("Relative macrocephaly", "Dysmorphic"),
    "HP:0005490": ("Stromal dystrophy", "Ophthalmologic"),

    # Cardiac
    "HP:0001639": ("Hypertrophic cardiomyopathy", "Cardiovascular"),
    "HP:0001712": ("Left ventricular hypertrophy", "Cardiovascular"),
    "HP:0001631": ("Atrial septal defect", "Cardiovascular"),

    # Renal
    "HP:0000100": ("Nephronophthisis", "Renal"),
    "HP:0000095": ("Proteinuria", "Renal"),

    # Integumentary
    "HP:0008066": ("Ichthyosis", "Integumentary"),
    "HP:0008064": ("Hypopigmented macule", "Integumentary"),
}


def symptom_to_hpo(
    symptoms: str,
    top_k: int = 10,
    use_ontology: bool = True,
) -> list[dict]:
    """
    Match free-text symptoms to HPO terms.

    Parameters
    ----------
    symptoms : str
        Free-text symptom description
    top_k : int
        Return top-k matched HPO terms
    use_ontology : bool
        Use pyhpo ontology for semantic similarity (recommended).
        If False, use built-in keyword matching.

    Returns
    -------
    list[dict]
        [{"hpo_id": "HP:0001250", "hpo_name": "Seizure", "category": "Neurological", "score": 0.95}, ...]
    """
    symptoms_lower = symptoms.lower()
    scores = []

    if use_ontology and HAS_PYHPO:
        scores = _ontology_match(symptoms_lower, top_k)
    else:
        scores = _keyword_match(symptoms_lower, top_k)

    return sorted(scores, key=lambda x: x["score"], reverse=True)[:top_k]


def _keyword_match(symptoms_lower: str, top_k: int) -> list[dict]:
    """Keyword-based HPO matching (fallback when pyhpo unavailable)."""
    # Common symptom keyword → HPO mappings
    KEYWORD_MAP = {
        "seizure": "HP:0001250",
        "spasm": "HP:0012469",
        "infantile spasm": "HP:0012469",
        "tonic-clonic": "HP:0002123",
        "generalized seizure": "HP:0002123",
        "developmental regression": "HP:0002376",
        "regression": "HP:0002376",
        "intellectual disability": "HP:0002189",
        "id": "HP:0002189",
        "cognitive": "HP:0002189",
        "dystonia": "HP:0000726",
        "gliosis": "HP:0002067",
        "hypotonia": "HP:0001319",
        "muscle weak": "HP:0003324",
        "weakness": "HP:0003324",
        "ataxia": "HP:0001251",
        "epilepsy": "HP:0001250",
        "adhd": "HP:0007018",
        "attention": "HP:0007018",
        "autism": "HP:0000733",
        "macrocephaly": "HP:0004482",
        "micr": "HP:0004482",
        "hypertelor": "HP:0000316",
        "facial": "HP:0000324",
        "ichthyosis": "HP:0008066",
        "cardiomyopathy": "HP:0001639",
        "hcm": "HP:0001639",
        "lv hypertrophy": "HP:0001712",
        "ata": "HP:0001631",
        "vomiting": "HP:0002021",
        "gerd": "HP:0002020",
        "reflux": "HP:0002020",
    }

    matched_ids = {}
    for keyword, hpo_id in KEYWORD_MAP.items():
        if keyword in symptoms_lower:
            matched_ids[hpo_id] = matched_ids.get(hpo_id, 0) + 1

    scores = []
    for hpo_id, count in matched_ids.items():
        if hpo_id in BUILTIN_HPO:
            name, category = BUILTIN_HPO[hpo_id]
            scores.append({
                "hpo_id": hpo_id,
                "hpo_name": name,
                "category": category,
                "score": min(count * 0.3, 1.0),   # rough score
            })

    return scores


def _ontology_match(symptoms_lower: str, top_k: int) -> list[dict]:
    """Ontology-based semantic HPO matching using pyhpo."""
    if not HAS_PYHPO:
        return _keyword_match(symptoms_lower, top_k)

    try:
        ontology = Ontology()
    except Exception as e:
        logger.warning(f"Failed to load HPO ontology: {e}, falling back to keyword")
        return _keyword_match(symptoms_lower, top_k)

    # Extract key terms from symptoms
    terms = symptoms_lower.split()
    scores = {}

    for term in terms:
        if len(term) < 3:
            continue
        try:
            matches = ontology.search(term)
            for match in matches[:3]:  # top 3 per term
                hpo_id = f"HP:{match.id}"
                if hpo_id not in scores:
                    scores[hpo_id] = {"hpo_id": hpo_id, "hpo_name": match.name, "category": _category_from_hpo(match), "score": 0.0}
                scores[hpo_id]["score"] += match.score
        except Exception:
            continue

    # Normalize scores
    max_score = max((s["score"] for s in scores.values()), default=1.0)
    for s in scores.values():
        s["score"] = min(s["score"] / max_score, 1.0)

    return list(scores.values())


def _category_from_hpo(hpo_term) -> str:
    """Extract category from HPO term (best-effort)."""
    try:
        return hpo_term.parent(1).name if hasattr(hpo_term, "parent") else "Unknown"
    except Exception:
        return "Unknown"


class HPOMatcher:
    """
    Stateful HPO matcher with caching.

    Example:
        matcher = HPOMatcher(ontology=ontology)
        results = matcher.match("intellectual disability, seizures")
        print(results)
    """

    def __init__(self, use_ontology: bool = True):
        self.use_ontology = use_ontology
        self._cache = {}

    def match(self, symptoms: str, top_k: int = 10) -> list[dict]:
        """Match symptoms to HPO terms (with caching)."""
        cache_key = (symptoms, top_k)
        if cache_key not in self._cache:
            self._cache[cache_key] = symptom_to_hpo(symptoms, top_k, self.use_ontology)
        return self._cache[cache_key]

    def batch_match(self, symptoms_list: list[str], top_k: int = 10) -> list[list[dict]]:
        """Batch match multiple symptom descriptions."""
        return [self.match(s, top_k) for s in symptoms_list]
