"""
ClinVar annotation module.

Merges ClinVar annotations and calculates the distance to the nearest
ClinVar variant in the same gene.

Usage:
    python -m seekrare.preprocess.clinvar_annotation input_annotated.csv clinvar.csv output.csv
"""

from __future__ import annotations

import bisect
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Union

import pandas as pd
from loguru import logger


CLNSIG_PRIORITY = {
    "Pathogenic": 1,
    "Likely_pathogenic": 2,
    "Pathogenic/Likely_pathogenic": 3,
    "Conflicting_classifications_of_pathogenicity": 4,
    "Uncertain_significance": 5,
    "Likely_benign": 6,
    "Benign/Likely_benign": 7,
    "Benign": 8,
}


def norm_chrom(x: str) -> str:
    """Normalize chromosome string."""
    if pd.isna(x) or str(x).strip() == "":
        return ""
    x = str(x).strip()
    if x.startswith("chr"):
        return x
    if x in ("X", "Y", "M", "MT"):
        return "chrM" if x in ("M", "MT") else "chr" + x
    return "chr" + x


def extract_gene(geneinfo: str) -> str:
    """Extract gene names from GENEINFO field (e.g. 'OR4F5:79501|A:1' → 'A;OR4F5')."""
    if pd.isna(geneinfo) or str(geneinfo).strip() == "":
        return ""
    genes = []
    for item in str(geneinfo).split("|"):
        g = item.split(":")[0].strip()
        if g:
            genes.append(g)
    return ";".join(sorted(set(genes)))


def join_unique(series: pd.Series) -> str:
    """Join unique non-empty values from a Series."""
    vals = [str(v) for v in series if v and str(v) not in ("nan", ".", "")]
    return ";".join(sorted(set(vals)))


def get_clinvar_rank(sig: str) -> int:
    """Get ClinVar significance ranking (lower = more pathogenic)."""
    if pd.isna(sig) or str(sig).strip() == "":
        return 99
    items = re.split(r"[;,]", str(sig))
    ranks = [CLNSIG_PRIORITY[i.strip()] for i in items if i.strip() in CLNSIG_PRIORITY]
    return min(ranks) if ranks else 98


def nearest_clinvar_distance(gene_name: str, chrom: str, pos: int, clin_pos_dict: dict) -> int:
    """Find nearest ClinVar variant distance in the same gene/chromosome."""
    if pd.isna(gene_name) or str(gene_name).strip() == "" or pd.isna(pos):
        return -1

    pos = int(pos)
    genes = [g.strip() for g in str(gene_name).split(";") if g.strip()]
    best = None

    for g in genes:
        key = (g, chrom)
        if key not in clin_pos_dict:
            continue
        pos_list = clin_pos_dict[key]
        i = bisect.bisect_left(pos_list, pos)
        candidates = []
        if i < len(pos_list):
            candidates.append(abs(pos_list[i] - pos))
        if i > 0:
            candidates.append(abs(pos_list[i - 1] - pos))
        if candidates:
            d = min(candidates)
            if best is None or d < best:
                best = d

    return best if best is not None else -1


def merge_filter_clinvar(
    annotated_csv: Union[str, Path],
    clinvar_csv: Union[str, Path],
    output_csv: Union[str, Path],
) -> pd.DataFrame:
    """
    Annotate with ClinVar and filter to clinically relevant variants.

    Parameters
    ----------
    annotated_csv : str or Path
        CSV from annotate_by_gtf (with gene_name column)
    clinvar_csv : str or Path
        ClinVar CSV (needs: GENEINFO, CHROM, POS, CLNSIG, MC, CLNDISDB)
    output_csv : str or Path
        Output CSV path

    Returns
    -------
    pd.DataFrame
        ClinVar-annotated and filtered DataFrame with columns:
        clinvar_sig, clinvar_mc, clinvar_hp, clinvar_min_distance
    """
    annotated_csv = str(annotated_csv)
    clinvar_csv = str(clinvar_csv)
    output_csv = str(output_csv)
    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)

    t0 = time.time()
    logger.info(f"Loading: annotated CSV ({annotated_csv}), ClinVar ({clinvar_csv})")

    df = pd.read_csv(annotated_csv, low_memory=False)
    clin = pd.read_csv(clinvar_csv, low_memory=False)

    # Normalize chromosomes
    df["CHROM_norm"] = df["CHROM"].apply(norm_chrom)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    clin["CHROM_norm"] = clin["CHROM"].apply(norm_chrom)
    clin["POS"] = pd.to_numeric(clin["POS"], errors="coerce")

    # Extract gene from ClinVar
    clin["gene"] = clin["GENEINFO"].apply(extract_gene)
    clin = clin.assign(gene=clin["gene"].str.split(";")).explode("gene")
    clin["gene"] = clin["gene"].astype(str).str.strip()

    # Build gene → ClinVar info dict
    clin_group = (
        clin.groupby("gene")
        .agg({"CLNSIG": join_unique, "MC": join_unique, "CLNDISDB": join_unique})
        .reset_index()
    )
    clin_dict = clin_group.set_index("gene").to_dict(orient="index")

    # Build gene+chrom → sorted POS list
    clin_pos_dict = defaultdict(list)
    for _, row in clin.dropna(subset=["gene", "POS"]).iterrows():
        clin_pos_dict[(row["gene"], row["CHROM_norm"])].append(int(row["POS"]))
    for key in clin_pos_dict:
        clin_pos_dict[key] = sorted(set(clin_pos_dict[key]))

    logger.info(f"  ClinVar entries: {len(clin)}, genes: {len(clin_dict)}, gene+chrom keys: {len(clin_pos_dict)}")

    # Match ClinVar info by gene_name
    def match_clinvar_info(gene_name: str):
        if pd.isna(gene_name) or str(gene_name).strip() == "":
            return "", "", ""
        genes = [g.strip() for g in str(gene_name).split(";") if g.strip()]
        sigs, mcs, hps = [], [], []
        for g in genes:
            if g in clin_dict:
                sigs.append(clin_dict[g]["CLNSIG"])
                mcs.append(clin_dict[g]["MC"])
                hps.append(clin_dict[g]["CLNDISDB"])
        return (
            ";".join(sorted(set([x for x in sigs if x]))),
            ";".join(sorted(set([x for x in mcs if x]))),
            ";".join(sorted(set([x for x in hps if x]))),
        )

    logger.info("  Matching ClinVar info...")
    clin_info = df["gene_name"].apply(match_clinvar_info)
    df["clinvar_sig"] = clin_info.apply(lambda x: x[0])
    df["clinvar_mc"] = clin_info.apply(lambda x: x[1])
    df["clinvar_hp"] = clin_info.apply(lambda x: x[2])

    logger.info("  Computing nearest ClinVar distance...")
    df["clinvar_min_distance"] = df.apply(
        lambda row: nearest_clinvar_distance(row.get("gene_name", ""), row.get("CHROM_norm", ""), row.get("POS", 0), clin_pos_dict),
        axis=1,
    )

    df["clinvar_rank"] = df["clinvar_sig"].apply(get_clinvar_rank)

    # Filter: in gene + has ClinVar significance
    mask = (
        df["in_gene"].astype(str).str.lower().isin(["true", "1"]) &
        df["clinvar_sig"].notna() &
        (df["clinvar_sig"].astype(str).str.strip() != "")
    )
    df_filtered = df[mask].copy()

    # Sort
    sort_cols = ["clinvar_rank", "clinvar_min_distance", "gene_name", "CHROM", "POS"]
    df_filtered = df_filtered.sort_values(
        by=[c for c in sort_cols if c in df_filtered.columns],
        ascending=[True] * 5,
    ).reset_index(drop=True)

    # Reorder: clinvar_min_distance first
    cols = list(df_filtered.columns)
    for c in ["clinvar_min_distance", "clinvar_rank", "CHROM_norm"]:
        if c in cols:
            cols.remove(c)
    cols.insert(0, "clinvar_min_distance")
    df_filtered = df_filtered[[c for c in cols if c in df_filtered.columns]]

    df_filtered.to_csv(output_csv, index=False)

    n_path = df_filtered["clinvar_sig"].str.contains("Pathogenic", na=False).sum()
    logger.info(
        f"  Done! {len(df):,} → {len(df_filtered):,} variants "
        f"(Pathogenic: {n_path:,}). Saved to {output_csv} ({time.time()-t0:.1f}s)"
    )

    return df_filtered


# Needed for time import
import time

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    merge_filter_clinvar(sys.argv[1], sys.argv[2], sys.argv[3])
