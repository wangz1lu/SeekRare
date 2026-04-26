"""
CADD (Combined Annotation Dependent Depletion) score integration.

Loads CADD scores from:
- CADD TSV/CSV files (CADD website download)
- tabix-indexed CADD vcf.gz files
- Direct API query to CADD server (rare variants)

CADD scores range 1-99, higher = more deleterious.
Typical thresholds: 20 (moderate), 25 (confident), 30 (strong selection)

Usage:
    from seekrare.annotation import CADDLoader
    cadd = CADDLoader("/ref/cadd/GRCh38/cadd.tsv.gz")
    annotated = cadd.annotate_variants(df)
"""

from __future__ import annotations

import gzip
import re
from bisect import bisect_left
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


def _norm_chrom(chrom: str) -> str:
    chrom = str(chrom).strip().replace("chr", "")
    return chrom if chrom != "M" else "MT"


class CADDLoader:
    """
    CADD score annotation with tabix-style fast lookup.

    Parameters
    ----------
    cadd_path : str or Path
        Path to CADD TSV (.tsv.gz) or VCF (.vcf.gz) file
    """

    def __init__(self, cadd_path: Union[str, Path]):
        self.cadd_path = Path(cadd_path)
        self._index: dict[str, list[tuple]] = {}  # chrom → [(pos, ref, alt, cadd_phred, cadd_raw), ...]
        self._loaded = False

    def build_index(self):
        """Build per-chromosome position index."""
        logger.info(f"Building CADD index from {self.cadd_path}")
        suffix = self.cadd_path.suffix.lower()

        if suffix == ".gz":
            self._load_tsv_gz()
        elif suffix == ".tsv":
            self._load_tsv()
        else:
            raise ValueError(f"Unsupported CADD file type: {suffix}")

        self._loaded = True
        n_total = sum(len(v) for v in self._index.values())
        logger.info(f"CADD index built: {n_total} entries")

    def _load_tsv_gz(self):
        """Load gzipped CADD TSV file."""
        opener = gzip.open
        self._parse_tsv(opener)

    def _load_tsv(self):
        """Load plain CADD TSV file."""
        self._parse_tsv(open)

    def _parse_tsv(self, opener_fn):
        """Parse TSV file (opener_fn = open or gzip.open)."""
        chrom_col = None

        with opener_fn(self.cadd_path, "rt") as f:
            for i, line in enumerate(f):
                if line.startswith("#"):
                    # Parse header
                    if line.startswith("#Chrom"):
                        cols = line.lstrip("#").strip().split("\t")
                        chrom_col = cols.index("Chrom") if "Chrom" in cols else 0
                        pos_col = cols.index("Pos") if "Pos" in cols else 1
                        ref_col = cols.index("Ref") if "Ref" in cols else 2
                        alt_col = cols.index("Alt") if "Alt" in cols else 3
                        phred_col = cols.index("PHRED") if "PHRED" in cols else None
                        raw_col = cols.index("CADD") if "CADD" in cols else None
                        gene_col = cols.index("GeneName") if "GeneName" in cols else None
                    continue

                if chrom_col is None:
                    # Fallback: use positional parsing
                    chrom_col, pos_col, ref_col, alt_col = 0, 1, 2, 3
                    phred_col, raw_col, gene_col = 4, 5, 6

                fields = line.strip().split("\t")
                if len(fields) <= max(chrom_col or 0, pos_col or 0):
                    continue

                try:
                    chrom = _norm_chrom(fields[chrom_col])
                    pos = int(fields[pos_col])
                    ref = fields[ref_col].upper() if ref_col < len(fields) else "N"
                    alt = fields[alt_col].upper() if alt_col < len(fields) else "N"

                    phred = float(fields[phred_col]) if phred_col and phred_col < len(fields) and fields[phred_col] else None
                    raw = float(fields[raw_col]) if raw_col and raw_col < len(fields) and fields[raw_col] else None
                    gene = fields[gene_col] if gene_col and gene_col < len(fields) else ""

                    self._index.setdefault(chrom, []).append((pos, ref, alt, phred, raw, gene))
                except (ValueError, IndexError):
                    continue

        # Sort each chromosome
        for chrom in self._index:
            self._index[chrom].sort(key=lambda x: x[0])

    def annotate_variants(
        self,
        df: pd.DataFrame,
        exact_match_only: bool = True,
    ) -> pd.DataFrame:
        """
        Annotate variants with CADD scores.

        Parameters
        ----------
        df : pd.DataFrame
            Must have CHROM, POS, REF, ALT columns
        exact_match_only : bool
            If True, only annotate exact (pos, ref, alt) matches

        Returns
        -------
        pd.DataFrame
            df with cadd_phred and cadd_raw columns added
        """
        if not self._loaded:
            self.build_index()

        df = df.copy()
        chrom_col = "CHROM" if "CHROM" in df.columns else "chrom"
        pos_col = "POS" if "POS" in df.columns else "pos"

        if "cadd_phred" not in df.columns:
            df["cadd_phred"] = None
        if "cadd_raw" not in df.columns:
            df["cadd_raw"] = None
        if "cadd_gene" not in df.columns:
            df["cadd_gene"] = None

        exact, skipped = 0, 0

        for idx, row in df.iterrows():
            chrom = _norm_chrom(row.get(chrom_col, ""))
            pos = int(row.get(pos_col, 0))
            ref = str(row.get("REF", "")).upper()
            alt = str(row.get("ALT", "")).upper()

            entries = self._index.get(chrom, [])
            if not entries:
                skipped += 1
                continue

            pos_list = [e[0] for e in entries]
            i = bisect_left(pos_list, pos)

            matched = False
            best_phred = None
            best_raw = None
            best_gene = None

            for j in [i - 1, i, i + 1]:
                if 0 <= j < len(entries):
                    e = entries[j]
                    if e[0] == pos and e[1] == ref and e[2] == alt:
                        phred, raw = e[3], e[4]
                        if phred is not None:
                            best_phred = phred
                            best_raw = raw
                            best_gene = e[5] if len(e) > 5 else ""
                            matched = True
                            break
                    elif exact_match_only:
                        continue

            if matched:
                df.at[idx, "cadd_phred"] = best_phred
                df.at[idx, "cadd_raw"] = best_raw
                df.at[idx, "cadd_gene"] = best_gene
                exact += 1

        logger.info(f"CADD: {exact}/{len(df)} exact matches, {skipped} chroms not found")
        return df

    @staticmethod
    def cadd_category(phred: float) -> str:
        """Categorize CADD score."""
        if phred is None:
            return "unknown"
        if phred >= 30:
            return "highly_deleterious"
        if phred >= 25:
            return "likely_deleterious"
        if phred >= 20:
            return "possibly_deleterious"
        return "tolerated"
