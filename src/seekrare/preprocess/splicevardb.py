"""
SpliceVARDB annotation for Stage 2.

Usage:
    from seekrare.preprocess.splicevardb import stage2_splicevardb_annotation
    stage2_splicevardb_annotation("stage1.csv", "splicevardb.tsv", "output.csv")
"""

from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


def _parse_variant_id(vid: str) -> tuple[str, int, str, str]:
    """Parse variant_id like '1-100573238-T-C' → (chr, pos, ref, alt)."""
    parts = vid.strip().split("-")
    if len(parts) != 4:
        raise ValueError(f"Cannot parse variant_id: {vid!r} (expected 'chr-pos-ref-alt')")
    return parts[0], int(parts[1]), parts[2], parts[3]


def _build_cpra(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Build CPRA key for matching."""
    return f"{chrom}:{pos}:{ref}:{alt}"


def stage2_splicevardb_annotation(
    input_csv: Union[str, Path],
    splicevardb_tsv: Union[str, Path],
    output_csv: Union[str, Path],
    chrom_format: str = "no_prefix",  # "chr" or "no_prefix"
) -> pd.DataFrame:
    """
    Annotate Stage 1 CSV with SpliceVARDB data.

    Parameters
    ----------
    input_csv : str   Stage 1 output CSV (with CHROM, POS, REF, ALT columns)
    splicevardb_tsv : str   SpliceVARDB TSV (header: variant_id, hg19, hg38, gene,
                           hgvs, method, classification, location, doi)
    output_csv : str   Annotated output CSV
    chrom_format : str   "chr" or "no_prefix" — how CHROM appears in input CSV

    Output columns added:
        splice_gene, splice_hgvs, splice_method, splice_classification,
        splice_location, splice_doi, splice_variant_id
    """
    input_csv = str(input_csv)
    splicevardb_tsv = str(splicevardb_tsv)
    output_csv = str(output_csv)

    logger.info(f"[SpliceVARDB] Loading Stage 1 CSV: {input_csv}")
    df = pd.read_csv(input_csv, dtype=str)
    logger.info(f"  Variants: {len(df)}")

    logger.info(f"[SpliceVARDB] Loading SpliceVARDB: {splicevardb_tsv}")
    sv = pd.read_csv(splicevardb_tsv, sep="\t", dtype=str)
    logger.info(f"  SpliceVARDB entries: {len(sv)}")

    # Build variant_id → row mapping from SpliceVARDB hg38 column
    # variant_id format in SpliceVARDB: "1-100573238-T-C"
    # hg38 column: "1-100107682-T-C"
    # We use hg38 column for matching (hg38 is GRCh38)
    sv_keys = {}
    for _, row in sv.iterrows():
        vid = str(row.get("variant_id", "")).strip()
        hg38 = str(row.get("hg38", "")).strip()
        key = hg38 if hg38 and hg38 != "nan" else vid
        if key and key != "nan":
            sv_keys[key] = row

    # Build CPRA keys for input CSV
    if chrom_format == "chr":
        df["_cpr"] = df.apply(
            lambda r: _build_cpra(str(r["CHROM"]).replace("chr", ""), int(r["POS"]), r["REF"], r["ALT"]),
            axis=1,
        )
    else:
        df["_cpr"] = df.apply(
            lambda r: _build_cpra(str(r["CHROM"]), int(r["POS"]), r["REF"], r["ALT"]),
            axis=1,
        )

    # Annotate
    records = []
    for _, row in df.iterrows():
        cpr = row["_cpr"]
        # Try exact CPRA match first
        hit = sv_keys.get(cpr)
        # If no hit, try to match via pos-ref-alt in hg38 format "chr-pos-ref-alt"
        if not hit:
            parts = cpr.split(":")
            if len(parts) == 4:
                # Build hg38-style key: chr-pos-ref-alt
                hg38_key = f"{parts[0]}-{parts[1]}-{parts[2]}-{parts[3]}"
                hit = sv_keys.get(hg38_key)
        records.append(hit)

    # Fill annotation columns
    df["splice_variant_id"] = [r["variant_id"] if r is not None else "" for r in records]
    df["splice_gene"] = [r["gene"] if r is not None else "" for r in records]
    df["splice_hgvs"] = [r["hgvs"] if r is not None else "" for r in records]
    df["splice_method"] = [r["method"] if r is not None else "" for r in records]
    df["splice_classification"] = [r["classification"] if r is not None else "" for r in records]
    df["splice_location"] = [r["location"] if r is not None else "" for r in records]
    df["splice_doi"] = [r["doi"] if r is not None else "" for r in records]

    df.drop(columns=["_cpr"], inplace=True)

    n_annotated = (df["splice_gene"] != "").sum()
    logger.info(f"[SpliceVARDB] Annotated: {n_annotated}/{len(df)} variants")

    df.to_csv(output_csv, index=False)
    logger.info(f"[SpliceVARDB] Saved: {output_csv}")

    return df


def stage2_splicevardb_annotation_by_cpra(
    input_csv: Union[str, Path],
    splicevardb_tsv: Union[str, Path],
    output_csv: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Alias for stage2_splicevardb_annotation with CPRA-based matching.
    Kept for backward compatibility.
    """
    return stage2_splicevardb_annotation(input_csv, splicevardb_tsv, output_csv)