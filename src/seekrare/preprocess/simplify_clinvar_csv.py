"""Simplify ClinVar-annotated CSV: clean columns, map stars, split CLNDISDB."""

import re
from typing import Optional

import numpy as np
import pandas as pd


STAR_MAP = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 0,
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
    "no_classification_provided": 0,
    "no_classification_for_the_single_variant": 0,
}


def _map_star(val) -> int:
    if pd.isna(val) or str(val).strip() == "":
        return 0
    return STAR_MAP.get(str(val).strip(), 0)


def _extract_hpo(val) -> Optional[str]:
    if pd.isna(val):
        return None
    m = re.search(r"HP:\d+", str(val))
    return m.group(0) if m else None


def _extract_omim(val) -> Optional[str]:
    if pd.isna(val):
        return None
    matches = re.findall(r"OMIM:(\d+)", str(val))
    if not matches:
        return None
    return "|".join([f"OMIM:{x}" for x in dict.fromkeys(matches)])


def _extract_orphanet(val) -> Optional[str]:
    if pd.isna(val):
        return None
    matches = re.findall(r"Orphanet:\d+", str(val))
    if not matches:
        return None
    return "|".join(dict.fromkeys(matches))


def _simplify_mc(val) -> Optional[str]:
    if pd.isna(val):
        return None
    text = str(val).strip()
    return text.split("|")[-1] if "|" in text else text


def simplify_clinvar_csv(df: pd.DataFrame) -> pd.DataFrame:
    """Apply all simplification rules to a ClinVar-annotated DataFrame."""
    df = df.copy()

    # 1. 删除 ORIGIN
    if "ORIGIN" in df.columns:
        df.drop(columns=["ORIGIN"], inplace=True)

    # 2. CLNREVSTAT → clinvarstar
    df["clinvarstar"] = df["CLNREVSTAT"].apply(_map_star)

    # 3. 删除 CLNHGVS
    if "CLNHGVS" in df.columns:
        df.drop(columns=["CLNHGVS"], inplace=True)

    # 4. CLNDN → diseasename
    if "CLNDN" in df.columns:
        df.rename(columns={"CLNDN": "diseasename"}, inplace=True)

    # 5. CLNSIG → significance，删除 benign / likely_benign
    if "CLNSIG" in df.columns:
        df.rename(columns={"CLNSIG": "significance"}, inplace=True)
        sl = df["significance"].astype(str).str.lower()
        df = df[~sl.isin(["benign", "likely_benign"])]

    # 6. 删除 CLNVC
    if "CLNVC" in df.columns:
        df.drop(columns=["CLNVC"], inplace=True)

    # 7. MC 列保留 | 右侧
    if "MC" in df.columns:
        df["MC"] = df["MC"].apply(_simplify_mc)

    # 8. CLNDISDB → HPO / OMIM / Orphanet
    if "CLNDISDB" in df.columns:
        df["HPO"]      = df["CLNDISDB"].apply(_extract_hpo)
        df["OMIM"]     = df["CLNDISDB"].apply(_extract_omim)
        df["Orphanet"] = df["CLNDISDB"].apply(_extract_orphanet)

    return df
