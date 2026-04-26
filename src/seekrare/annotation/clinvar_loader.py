"""
ClinVar annotation loader.

Loads ClinVar VCF or CSV, performs chromosomal distance-based
annotation of the nearest ClinVar variant for each query variant.

Supports:
- ClinVar VCF (GRCh38)
- ClinVar CSV dump
- Lookup by (chrom, pos, ref, alt) exact match
- Distance-based fallback for intronic/flanking variants

Usage:
    from seekrare.annotation import ClinVarLoader
    clinvar = ClinVarLoader("clinvar_2024.vcf.gz")
    annotated = clinvar.annotate_variants(df)   # df needs CHROM, POS, REF, ALT
"""

from __future__ import annotations

import gzip
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


# ClinVar Clinical Significance → priority (lower = more significant)
CLNSIG_PRIORITY = {
    "Pathogenic": 1,
    "Likely pathogenic": 2,
    "Uncertain significance": 5,
    "Likely benign": 7,
    "Benign": 8,
    "Pathogenic/Likely pathogenic": 1,
    "Likely pathogenic/Pathogenic": 1,
    "Conflicting interpretations of pathogenicity": 6,
    "not provided": 4,
    "other": 5,
}


def _norm_chrom(chrom: str) -> str:
    """Normalize chromosome name to simple form."""
    chrom = str(chrom).strip().replace("chr", "")
    return chrom if chrom not in ("M", "MT") else "M"


class ClinVarLoader:
    """
    ClinVar annotation loader with distance-based lookup.

    Parameters
    ----------
    clinvar_path : str or Path
        Path to ClinVar VCF.gz or CSV file
    genome : str
        "GRCh38" or "GRCh37"
    """

    def __init__(self, clinvar_path: Union[str, Path], genome: str = "GRCh38"):
        self.clinvar_path = Path(clinvar_path)
        self.genome = genome
        self._index: dict[str, dict] = {}  # chrom → {pos_list: [(pos, ann_dict), ...]}
        self._loaded = False

    def build_index(self):
        """Build per-chromosome sorted position index."""
        logger.info(f"Building ClinVar index from {self.clinvar_path}")
        suffix = self.clinvar_path.suffix.lower()

        if suffix == ".gz" and self._is_vcf():
            self._load_vcf_index()
        elif suffix == ".csv":
            self._load_csv_index()
        elif suffix in (".vcf", ".vcf.gz"):
            self._load_vcf_index()
        else:
            raise ValueError(f"Unsupported ClinVar file format: {suffix}")

        self._loaded = True
        n_total = sum(len(v) for v in self._index.values())
        logger.info(f"ClinVar index built: {n_total} variants across {len(self._index)} chromosomes")

    def _is_vcf(self) -> bool:
        with gzip.open(self.clinvar_path, "rt") as f:
            first = next((l for l in f if l.strip() and not l.startswith("#")), "")
            return "##" in first or "CHROM" in first

    def _load_vcf_index(self):
        """Load from ClinVar VCF.gz (gzipped text or compressed)."""
        prefix = "chr" if self.genome == "GRCh38" else ""

        with gzip.open(self.clinvar_path, "rt") as f:
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

                info_str = fields[7]
                info = dict(x.split("=", 1) if "=" in x else (x, True)
                            for x in info_str.split(";"))

                # Extract key annotations from INFO field
                clnsig = info.get("CLNSIG", "")
                clndisdb = info.get("CLNDISDB", "")
                clnrevstat = info.get("CLNREVSTAT", "")
                gene = info.get("GENEINFO", "").split(":")[0] if "GENEINFO" in info else ""
                mc = info.get("MC", "")  # Molecular consequence
                rsid = info.get("RS", info.get("dbSNPBuildID", ""))

                # Multiple ALTs → multiple entries
                alts = alt.split(",") if "," in alt else [alt]
                refs = ref.split(",") if "," in ref else [ref]

                for i, (a, r) in enumerate(zip(alts, refs)):
                    key = (chrom, pos, r, a)
                    self._index.setdefault(chrom, []).append({
                        "pos": pos,
                        "ref": r,
                        "alt": a,
                        "clinvar_sig": clnsig.split("|")[i] if "|" in clnsig else clnsig,
                        "clinvar_mc": mc.split("|")[i] if "|" in mc else mc,
                        "clinvar_gene": gene.split("|")[i] if "|" in gene else gene,
                        "clinvar_rsid": rsid.split("|")[i] if "|" in str(rsid) else str(rsid),
                        "clinvar_disdb": clndisdb.split("|")[i] if "|" in clndisdb else clndisdb,
                        "clinvar_revstat": clnrevstat.split("|")[i] if "|" in clnrevstat else clnrevstat,
                    })

        # Sort each chromosome by position
        for chrom in self._index:
            self._index[chrom] = sorted(self._index[chrom], key=lambda x: x["pos"])

    def _load_csv_index(self):
        """Load from ClinVar CSV dump."""
        df = pd.read_csv(self.clinvar_path, low_memory=False, dtype={"#CHROM": str})
        df.rename(columns={"#CHROM": "CHROM"}, inplace=True)

        for _, row in df.iterrows():
            chrom = _norm_chrom(row.get("CHROM", ""))
            pos = int(row.get("POS", 0))
            ref = str(row.get("REF", "")).upper()
            alt = str(row.get("ALT", "")).upper()

            ann = {
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "clinvar_sig": str(row.get("CLNSIG", "")),
                "clinvar_gene": str(row.get("GENEINFO", "")).split(":")[0],
                "clinvar_rsid": str(row.get("RS", "")),
            }
            self._index.setdefault(chrom, []).append(ann)

        for chrom in self._index:
            self._index[chrom] = sorted(self._index[chrom], key=lambda x: x["pos"])

    def annotate_variants(
        self,
        df: pd.DataFrame,
        max_distance: int = 5000,
        exact_match_only: bool = False,
    ) -> pd.DataFrame:
        """
        Annotate variants with ClinVar information.

        Parameters
        ----------
        df : pd.DataFrame
            Must have CHROM, POS, REF, ALT columns
        max_distance : int
            Max distance (bp) for nearest-variant fallback lookup
        exact_match_only : bool
            If True, only annotate exact matches (no distance lookup)

        Returns
        -------
        pd.DataFrame
            Original df with new ClinVar columns
        """
        if not self._loaded:
            self.build_index()

        df = df.copy()
        chrom_col = "CHROM" if "CHROM" in df.columns else "chrom"
        pos_col = "POS" if "POS" in df.columns else "pos"

        # Initialize output columns
        for col in ["clinvar_sig", "clinvar_mc", "clinvar_gene",
                    "clinvar_rsid", "clinvar_disdb", "clinvar_revstat",
                    "clinvar_min_distance", "clinvar_match_type"]:
            if col not in df.columns:
                df[col] = None

        exact, distant = 0, 0

        for idx, row in df.iterrows():
            chrom = _norm_chrom(row.get(chrom_col, ""))
            pos = int(row.get(pos_col, 0))
            ref = str(row.get("REF", "")).upper()
            alt = str(row.get("ALT", "")).upper()

            chrom_entries = self._index.get(chrom, [])
            if not chrom_entries:
                continue

            # Exact match lookup
            match = None
            key = (pos, ref, alt)
            pos_list = [e["pos"] for e in chrom_entries]
            i = bisect_left(pos_list, pos)

            # Check exact positions around i
            for j in [i - 1, i, i + 1]:
                if 0 <= j < len(chrom_entries):
                    e = chrom_entries[j]
                    if e["pos"] == pos and e["ref"] == ref and e["alt"] == alt:
                        match = (e, 0, "exact")
                        break

            if match:
                e, dist, mtype = match
                exact += 1
            elif not exact_match_only and max_distance > 0:
                # Nearest variant within max_distance
                nearest = self._find_nearest(chrom_entries, pos, max_distance)
                if nearest:
                    e, dist = nearest
                    match = (e, dist, "nearest")
                    distant += 1

            if match:
                e, dist, mtype = match
                df.at[idx, "clinvar_sig"] = e.get("clinvar_sig")
                df.at[idx, "clinvar_mc"] = e.get("clinvar_mc")
                df.at[idx, "clinvar_gene"] = e.get("clinvar_gene")
                df.at[idx, "clinvar_rsid"] = e.get("clinvar_rsid")
                df.at[idx, "clinvar_disdb"] = e.get("clinvar_disdb")
                df.at[idx, "clinvar_revstat"] = e.get("clinvar_revstat")
                df.at[idx, "clinvar_min_distance"] = dist
                df.at[idx, "clinvar_match_type"] = mtype

        logger.info(f"ClinVar: {exact} exact + {distant} nearest / {len(df)} variants annotated")
        return df

    def _find_nearest(
        self, entries: list[dict], pos: int, max_dist: int
    ) -> Optional[tuple[dict, int]]:
        """Find nearest ClinVar variant within max_distance."""
        pos_list = [e["pos"] for e in entries]
        i = bisect_left(pos_list, pos)

        best = None
        best_dist = max_dist + 1

        for j in [i - 1, i]:
            if 0 <= j < len(entries):
                d = abs(entries[j]["pos"] - pos)
                if d <= max_dist and d < best_dist:
                    best_dist = d
                    best = entries[j]

        return (best, best_dist) if best else None
