"""
Gene annotation from NCBI GTF file.

Uses a sweep-line algorithm to annotate each variant with:
- gene_name
- feature_type (CDS, exon, etc.)
- in_gene flag

Usage:
    python -m seekrare.preprocess.gene_annotation input.csv genomic.gtf output.csv
    # or programmatically:
    from seekrare.preprocess import annotate_by_gtf
    df = annotate_by_gtf("input.csv", "genomic.gtf", "output.csv")
"""

from __future__ import annotations

import os
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Union

import pandas as pd
from loguru import logger


# NCBI accession → chromosome
NCBI_TO_CHR = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9":  "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9":  "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
    "NC_012920.1":  "chrM",
}

FEATURE_PRIORITY = {
    "stop_codon": 0,
    "start_codon": 1,
    "CDS": 2,
    "exon": 3,
    "transcript": 4,
    "gene": 5,
}

MAIN_CHROMOSOMES = {
    f"chr{i}" for i in range(1, 23)
} | {"chrX", "chrY", "chrM"}


def normalize_csv_chrom(chrom: str) -> str:
    """Normalize CSV chromosome: chr11_KI270721v1_random → chr11."""
    chrom = str(chrom).strip()
    m = re.match(r"^(chr[0-9XYM]+)", chrom)
    return m.group(1) if m else chrom


def normalize_gtf_chrom(chrom: str) -> str:
    """Normalize GTF chromosome: NC_000001.11 → chr1."""
    chrom = str(chrom).strip()
    return NCBI_TO_CHR.get(chrom, chrom)


def parse_gene_name(attr: str) -> str:
    """Extract gene_name from GTF attribute field."""
    for pattern in [
        r'gene_name\s+"([^"]+)"',
        r'gene\s+"([^"]+)"',
        r'Name\s+"([^"]+)"',
        r'gene_id\s+"([^"]+)"',
    ]:
        m = re.search(pattern, attr)
        if m:
            return m.group(1)
    return ""


def feature_rank(feature: str) -> int:
    return FEATURE_PRIORITY.get(feature, 99)


def load_gtf(gtf_path: str) -> dict[str, list[dict]]:
    """Load GTF into chromosome-indexed interval lists."""
    t0 = time.time()
    gtf_by_chrom = defaultdict(list)

    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            chrom_raw, feature, start_str, end_str = parts[0], parts[2], parts[3], parts[4]
            attr = parts[8]

            if feature not in FEATURE_PRIORITY:
                continue

            chrom = normalize_gtf_chrom(chrom_raw)
            if chrom not in MAIN_CHROMOSOMES:
                continue

            try:
                start, end = int(start_str), int(end_str)
            except ValueError:
                continue

            gene_name = parse_gene_name(attr)
            gtf_by_chrom[chrom].append({
                "start": start,
                "end": end,
                "gene_name": gene_name,
                "feature_type": feature,
                "rank": feature_rank(feature),
            })

    for chrom in gtf_by_chrom:
        gtf_by_chrom[chrom].sort(key=lambda x: x["start"])

    logger.info(f"  GTF loaded: {len(gtf_by_chrom)} chromosomes, {time.time()-t0:.1f}s")
    return dict(gtf_by_chrom)


def annotate_by_gtf(
    input_csv: Union[str, Path],
    gtf_file: Union[str, Path],
    output_csv: Union[str, Path],
) -> pd.DataFrame:
    """
    Annotate CSV variants with gene information from GTF.

    Parameters
    ----------
    input_csv : str or Path
        CSV from vcf_to_gt_csv (needs CHROM, POS columns)
    gtf_file : str or Path
        NCBI genomic.gtf file
    output_csv : str or Path
        Output CSV path

    Returns
    -------
    pd.DataFrame
        Input CSV with added columns: in_gene, gene_name, feature_type
    """
    input_csv = str(input_csv)
    gtf_file = str(gtf_file)
    output_csv = str(output_csv)
    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)

    t0 = time.time()
    logger.info(f"Annotating: {input_csv} with GTF: {gtf_file}")

    # Load CSV
    df = pd.read_csv(input_csv)
    df.columns = [c.lstrip("#") for c in df.columns]

    if "CHROM" not in df.columns or "POS" not in df.columns:
        raise ValueError("Input CSV must contain CHROM and POS columns")

    df["CHROM"] = df["CHROM"].astype(str).str.strip()
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["_CHROM_KEY"] = df["CHROM"].apply(normalize_csv_chrom)

    logger.info(f"  Variants: {len(df):,}")

    # Load GTF
    gtf_by_chrom = load_gtf(gtf_file)

    # Annotate
    df["in_gene"] = False
    df["gene_name"] = ""
    df["feature_type"] = ""

    for chrom in sorted(df["_CHROM_KEY"].dropna().unique()):
        snp_idx = df["_CHROM_KEY"] == chrom
        intervals = gtf_by_chrom.get(chrom, [])

        if not intervals:
            logger.info(f"  {chrom}: no GTF records, skipping")
            continue

        snp_sub = df.loc[snp_idx].sort_values("POS")
        active = []
        j = 0
        n_int = len(intervals)

        hit_gene = {}
        hit_feat = {}

        for idx, row in snp_sub.iterrows():
            pos = int(row["POS"]) if not pd.isna(row["POS"]) else None
            if pos is None:
                hit_gene[idx] = ""
                hit_feat[idx] = ""
                continue

            # Add intervals starting at or before pos
            while j < n_int and intervals[j]["start"] <= pos:
                active.append(intervals[j])
                j += 1

            # Remove intervals ended before pos
            active = [x for x in active if x["end"] >= pos]

            if not active:
                hit_gene[idx] = ""
                hit_feat[idx] = ""
                continue

            best = sorted(active, key=lambda x: x["rank"])[0]
            genes = sorted(set(x["gene_name"] for x in active if x["gene_name"]))
            hit_gene[idx] = ";".join(genes)
            hit_feat[idx] = best["feature_type"]

        df.loc[list(hit_gene.keys()), "gene_name"] = list(hit_gene.values())
        df.loc[list(hit_feat.keys()), "feature_type"] = list(hit_feat.values())
        df.loc[list(hit_gene.keys()), "in_gene"] = True

        n_hits = sum(bool(v) for v in hit_gene.values())
        logger.info(f"  {chrom}: {snp_idx.sum():,} variants, {n_hits:,} annotated")

    # Clean up and save
    df = df.drop(columns=["_CHROM_KEY"])

    front = ["in_gene", "gene_name", "feature_type"]
    cols = front + [c for c in df.columns if c not in front]
    df = df[cols]

    df = df.sort_values(["in_gene", "CHROM", "POS"], ascending=[False, True, True]).reset_index(drop=True)
    df.to_csv(output_csv, index=False)

    n_in_gene = df["in_gene"].sum()
    logger.info(f"  Done! {n_in_gene:,}/{len(df):,} variants annotated. Saved to {output_csv} ({time.time()-t0:.1f}s)")

    return df


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    annotate_by_gtf(sys.argv[1], sys.argv[2], sys.argv[3])
