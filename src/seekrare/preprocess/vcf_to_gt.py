"""
VCF → GT CSV converter.

Extracts CHROM, POS, REF, ALT, and per-sample GT from VCF files.
Supports .vcf and .vcf.gz.

Usage:
    python -m seekrare.preprocess.vcf_to_gt input.vcf output.csv
    # or programmatically:
    from seekrare.preprocess import vcf_to_gt_csv
    vcf_to_gt_csv("input.vcf", "output.csv")
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path
from typing import Union

import pandas as pd
from loguru import logger


def open_vcf(path: str):
    """Open VCF file (handles .gz)."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def extract_gt(sample_value: str, format_fields: list[str]) -> str:
    """Extract GT field from a sample VCF column."""
    if sample_value in (".", "./.", ".|."):
        return sample_value

    values = sample_value.split(":")
    if "GT" not in format_fields or "GT" not in values:
        return ""

    gt_idx = format_fields.index("GT")
    if gt_idx >= len(values):
        return ""

    return values[gt_idx]


def vcf_to_gt_csv(
    vcf_path: Union[str, Path],
    output_csv: Union[str, Path],
    sample_names: list[str] | None = None,
) -> pd.DataFrame:
    """
    Convert VCF to GT CSV.

    Parameters
    ----------
    vcf_path : str or Path
        Input VCF or VCF.gz file
    output_csv : str or Path
        Output CSV path
    sample_names : list[str], optional
        Specific sample names to extract. If None, extracts all.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: CHROM, POS, REF, ALT, <sample1>, <sample2>, ...
    """
    vcf_path = str(vcf_path)
    output_csv = str(output_csv)

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)

    logger.info(f"Converting VCF: {vcf_path} → {output_csv}")

    rows = []
    file_sample_names = []

    with open_vcf(vcf_path) as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.lstrip("#").split("\t")
                file_sample_names = header[9:]
                continue

            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]
            fmt = parts[8]
            format_fields = fmt.split(":") if fmt else []

            row = {"CHROM": str(chrom), "POS": int(float(pos)), "REF": str(ref), "ALT": str(alt)}

            sample_values = parts[9:]
            names_to_use = sample_names if sample_names else file_sample_names

            for i, sample_name in enumerate(names_to_use):
                if i < len(sample_values):
                    row[sample_name] = extract_gt(sample_values[i], format_fields)

            rows.append(row)

    df = pd.DataFrame(rows)
    df["CHROM"] = df["CHROM"].astype(str)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df = df.sort_values(["CHROM", "POS"]).reset_index(drop=True)

    df.to_csv(output_csv, index=False)

    logger.info(
        f"  Variants: {len(df):,}, Samples: {len(df.columns) - 4}, "
        f"Output columns: {list(df.columns)}"
    )

    return df


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    vcf_to_gt_csv(sys.argv[1], sys.argv[2])
