"""
SpliceAI / dbscSNV — splicing impact annotation.

SpliceAI (delta score 0-1, higher = stronger splice effect):
- SpliceAI is an deep learning model from Illumina
- delta_score > 0.5 (acceptor gain/loss) suggests strong splicing effect
- delta_score > 0.2 (acceptor gain/loss) is a moderate threshold

dbscSNV (AdaBoost and Random Forest scores):
- Ada score: 0-1, higher = more likely to affect splicing
- RF score: 0-1, higher = more likely to affect splicing
- Both > 0.6 considered pathogenic

Usage:
    from seekrare.annotation import SpliceAILoader
    splice = SpliceAILoader("/ref/spliceai/spliceai_scores.snv.hg38.vcf.gz")
    annotated = splice.annotate_variants(df)
"""

from __future__ import annotations

import gzip
from bisect import bisect_left
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


def _norm_chrom(chrom: str) -> str:
    chrom = str(chrom).strip().replace("chr", "")
    return chrom if chrom != "M" else "MT"


class SpliceAILoader:
    """
    SpliceAI VCF annotation loader.

    Parameters
    ----------
    vcf_path : str or Path
        Path to SpliceAI VCF (.vcf.gz)
    """

    def __init__(self, vcf_path: Union[str, Path]):
        self.vcf_path = Path(vcf_path)
        self._index: dict[str, list[tuple]] = {}
        self._loaded = False

    def build_index(self):
        """Build index from SpliceAI VCF."""
        logger.info(f"Building SpliceAI index from {self.vcf_path}")
        opener = gzip.open if str(self.vcf_path).endswith(".gz") else open

        with opener(self.vcf_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue

                chrom = _norm_chrom(fields[0])
                pos = int(fields[1])
                ref = fields[3].upper()
                alt = fields[4].upper()

                info = {}
                for item in fields[7].split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info[k] = v

                # SpliceAI INFO format: SPLICEAI=SNP|AA|DA|AG|DG...
                # delta_score fields: acceptor_gain, acceptor_loss, donor_gain, donor_loss
                splice_str = info.get("SPLICEAI", "")
                if splice_str:
                    parts = splice_str.split("|")
                    if len(parts) >= 5:
                        self._index.setdefault(chrom, []).append((
                            pos, ref, alt,
                            float(parts[1]) if parts[1] != "." else None,  # acceptor_gain
                            float(parts[2]) if parts[2] != "." else None,  # acceptor_loss
                            float(parts[3]) if parts[3] != "." else None,  # donor_gain
                            float(parts[4]) if parts[4] != "." else None,  # donor_loss
                        ))

        for chrom in self._index:
            self._index[chrom].sort(key=lambda x: x[0])

        n_total = sum(len(v) for v in self._index.values())
        logger.info(f"SpliceAI index built: {n_total} entries")

    def annotate_variants(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotate variants with SpliceAI delta scores.

        Returns df with columns: spliceai_*
        """
        if not self._loaded:
            self.build_index()
            self._loaded = True

        df = df.copy()
        chrom_col = "CHROM" if "CHROM" in df.columns else "chrom"
        pos_col = "POS" if "POS" in df.columns else "pos"

        for col in ["spliceai_acceptor_gain", "spliceai_acceptor_loss",
                    "spliceai_donor_gain", "spliceai_donor_loss",
                    "spliceai_max_delta", "spliceai_category"]:
            if col not in df.columns:
                df[col] = None

        matched = 0
        for idx, row in df.iterrows():
            chrom = _norm_chrom(row.get(chrom_col, ""))
            pos = int(row.get(pos_col, 0))
            ref = str(row.get("REF", "")).upper()
            alt = str(row.get("ALT", "")).upper()

            entries = self._index.get(chrom, [])
            if not entries:
                continue

            pos_list = [e[0] for e in entries]
            i = bisect_left(pos_list, pos)

            for j in [i - 1, i, i + 1]:
                if 0 <= j < len(entries):
                    e = entries[j]
                    if e[0] == pos and e[1] == ref and e[2] == alt:
                        ag, al, dg, dl = e[3], e[4], e[5], e[6]
                        max_delta = max(x for x in [ag, al, dg, dl] if x is not None)

                        df.at[idx, "spliceai_acceptor_gain"] = ag
                        df.at[idx, "spliceai_acceptor_loss"] = al
                        df.at[idx, "spliceai_donor_gain"] = dg
                        df.at[idx, "spliceai_donor_loss"] = dl
                        df.at[idx, "spliceai_max_delta"] = max_delta
                        df.at[idx, "spliceai_category"] = self._categorize(max_delta)
                        matched += 1
                        break

        logger.info(f"SpliceAI: {matched}/{len(df)} variants annotated")
        return df

    @staticmethod
    def _categorize(delta: float) -> str:
        if delta is None:
            return "unknown"
        if delta >= 0.5:
            return "strong_splicing_effect"
        if delta >= 0.2:
            return "moderate_splicing_effect"
        if delta >= 0.1:
            return "weak_splicing_effect"
        return "minimal_splicing_effect"


class dbscSNVLoader:
    """
    dbscSNV (database of splice consequence) loader.

    Two scores: AdaBoost (ADA) and Random Forest (RF), both 0-1.
    Pathogenic threshold: > 0.6 for both scores.

    Parameters
    ----------
    dbsc_path : str or Path
        Path to dbscSNV ANNOVAR file (dbscSNV.txt.gz)
    """

    def __init__(self, dbsc_path: Union[str, Path]):
        self.dbsc_path = Path(dbsc_path)
        self._index: dict[str, list[tuple]] = {}
        self._loaded = False

    def build_index(self):
        """Build index from dbscSNV ANNOVAR file."""
        logger.info(f"Building dbscSNV index from {self.dbsc_path}")
        opener = gzip.open if str(self.dbsc_path).endswith(".gz") else open

        with opener(self.dbsc_path, "rt") as f:
            header = f.readline()  # skip header
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 11:
                    continue

                chrom = _norm_chrom(fields[0])
                pos = int(fields[1])
                ref = fields[2].upper()
                alt = fields[3].upper()
                ada = float(fields[9]) if fields[9] != "." else None
                rf = float(fields[10]) if fields[10] != "." else None

                self._index.setdefault(chrom, []).append((pos, ref, alt, ada, rf))

        for chrom in self._index:
            self._index[chrom].sort(key=lambda x: x[0])

        n_total = sum(len(v) for v in self._index.values())
        logger.info(f"dbscSNV index built: {n_total} entries")

    def annotate_variants(self, df: pd.DataFrame) -> pd.DataFrame:
        """Annotate with dbscSNV ADA and RF scores."""
        if not self._loaded:
            self.build_index()
            self._loaded = True

        df = df.copy()
        chrom_col = "CHROM" if "CHROM" in df.columns else "chrom"
        pos_col = "POS" if "POS" in df.columns else "pos"

        for col in ["dbscsnv_ada", "dbscsnv_rf", "dbscsnv_pathogenic"]:
            if col not in df.columns:
                df[col] = None

        matched = 0
        for idx, row in df.iterrows():
            chrom = _norm_chrom(row.get(chrom_col, ""))
            pos = int(row.get(pos_col, 0))
            ref = str(row.get("REF", "")).upper()
            alt = str(row.get("ALT", "")).upper()

            entries = self._index.get(chrom, [])
            if not entries:
                continue

            pos_list = [e[0] for e in entries]
            i = bisect_left(pos_list, pos)

            for j in [i - 1, i, i + 1]:
                if 0 <= j < len(entries):
                    e = entries[j]
                    if e[0] == pos and e[1] == ref and e[2] == alt:
                        ada, rf = e[3], e[4]
                        df.at[idx, "dbscsnv_ada"] = ada
                        df.at[idx, "dbscsnv_rf"] = rf
                        df.at[idx, "dbscsnv_pathogenic"] = (
                            "pathogenic" if (ada and ada > 0.6 and rf and rf > 0.6) else "uncertain"
                        )
                        matched += 1
                        break

        logger.info(f"dbscSNV: {matched}/{len(df)} variants annotated")
        return df
