"""
SpliceVARDB annotation for Stage 2 — simplified.

匹配方式：使用 SpliceVARDB 的 hg38 列（格式 "chr-pos-ref-alt"）
只注释 classification 列，输出单列 splicevardb。

Usage:
    from seekrare.preprocess.splicevardb import stage2_splicevardb_annotation
    stage2_splicevardb_annotation("stage1.csv", "splicevardb.tsv", "output.csv")
"""

import pandas as pd
from loguru import logger
from pathlib import Path
from typing import Union


def stage2_splicevardb_annotation(
    input_csv: Union[str, Path],
    splicevardb_tsv: Union[str, Path],
    output_csv: Union[str, Path],
) -> pd.DataFrame:
    """
    Annotate Stage 1 CSV with SpliceVARDB classification.

    Parameters
    ----------
    input_csv : str   Stage 1 output CSV (with CHROM, POS, REF, ALT)
    splicevardb_tsv : str   SpliceVARDB TSV (header: variant_id, hg19, hg38, ...)
    output_csv : str   Annotated output CSV

    Output column added:
        splicevardb — SpliceVARDB classification (e.g. "Splice-altering", "Low-frequency")
                       Empty string if not found.
    """
    input_csv = str(input_csv)
    splicevardb_tsv = str(splicevardb_tsv)
    output_csv = str(output_csv)

    logger.info(f"[SpliceVARDB] Loading Stage 1: {input_csv}")
    df = pd.read_csv(input_csv, dtype=str)
    logger.info(f"  Variants: {len(df)}")

    logger.info(f"[SpliceVARDB] Loading SpliceVARDB: {splicevardb_tsv}")
    sv = pd.read_csv(splicevardb_tsv, sep="\t", dtype=str)
    logger.info(f"  Records: {len(sv)}")

    # Build hg38 → classification dict
    # hg38 format: "1-100107682-T-C"
    sv_dict = {}
    for _, row in sv.iterrows():
        hg38_val = str(row.get("hg38", "")).strip()
        classification = str(row.get("classification", "")).strip()
        if hg38_val and hg38_val not in ("nan", ""):
            sv_dict[hg38_val] = classification

    # Build CPRA from Stage 1: chr-pos-ref-alt (no chr prefix)
    df["_cpr"] = (
        df["CHROM"].astype(str).str.replace("chr", "", regex=False)
        + "-" + df["POS"].astype(str)
        + "-" + df["REF"].astype(str)
        + "-" + df["ALT"].astype(str)
    )

    df["splicevardb"] = df["_cpr"].map(sv_dict).fillna("")
    df.drop(columns=["_cpr"], inplace=True)

    n_hit = (df["splicevardb"] != "").sum()
    logger.info(f"[SpliceVARDB] {n_hit}/{len(df)} variants annotated")

    df.to_csv(output_csv, index=False)
    logger.info(f"[SpliceVARDB] Saved: {output_csv}")

    return df