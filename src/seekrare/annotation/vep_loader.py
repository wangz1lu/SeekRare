"""
VEP (Variant Effect Predictor) annotation loader.

Loads VEP annotations from cache files (--cache --everything --json).
Supports: consequence, gene symbol, protein change, SIFT, PolyPhen,
CADD, gnomAD,ClinVar, OMIM, UniProt.

Usage:
    from seekrare.annotation import VEPAnnotationLoader
    vep = VEPAnnotationLoader(cache_dir="/ref/vep_cache")
    annotations = vep.annotate_variants(variants_df)   # variants_df needs CHROM, POS, REF, ALT
    # or
    vep.build_index()   # pre-index for faster lookup
"""

from __future__ import annotations

import gzip
import json
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


class VEPAnnotationLoader:
    """
    VEP annotation loader from cache or VCF output.

    Parameters
    ----------
    cache_dir : str or Path
        Directory containing VEP cache files (per chromosome VEP_{chrom}.json.gz)
    genome_version : str
        "GRCh37" or "GRCh38"
    """

    # VEP consequence term → simplified category
    CONSEQUENCE_CATEGORIES = {
        "transcript_ablation": "HIGH",
        "splice_acceptor_variant": "HIGH",
        "splice_donor_variant": "HIGH",
        "stop_gained": "HIGH",
        "frameshift_variant": "HIGH",
        "stop_lost": "HIGH",
        "start_lost": "HIGH",
        "transcript_amplification": "HIGH",
        "regulatory_region_variant": "MODERATE",
        "TF_binding_site_variant": "MODERATE",
        "protein_altering_variant": "MODERATE",
        "splice_region_variant": "MODERATE",
        "inframe_insertion": "MODERATE",
        "inframe_deletion": "MODERATE",
        "missense_variant": "MODERATE",
        "protein_coding_transcript": "MODERATE",
        "3_prime_UTR_variant": "LOW",
        "5_prime_UTR_variant": "LOW",
        "intron_variant": "LOW",
        "non_coding_transcript_variant": "LOW",
        "upstream_gene_variant": "LOW",
        "downstream_gene_variant": "LOW",
        "intergenic_variant": "LOW",
        "synonymous_variant": "LOW",
        "stop_retained_variant": "LOW",
    }

    def __init__(self, cache_dir: Union[str, Path], genome_version: str = "GRCh38"):
        self.cache_dir = Path(cache_dir)
        self.genome_version = genome_version
        self._index: dict[tuple, dict] = {}
        self._loaded_chroms: set[str] = set()

    def build_index(self, chroms: Optional[list[str]] = None):
        """
        Pre-load VEP annotations into memory index for fast lookup.

        Parameters
        ----------
        chroms : list[str], optional
            Specific chromosomes to load (e.g., ["1", "2", "X"]). 
            If None, loads all found in cache_dir.
        """
        if chroms is None:
            chrom_files = sorted(self.cache_dir.glob("VEP_*.json.gz"))
        else:
            chrom_files = [self.cache_dir / f"VEP_{c}.json.gz" for c in chroms]

        logger.info(f"Building VEP index from {len(chrom_files)} chromosome files")

        for fpath in chrom_files:
            self._load_chrom_file(fpath)

        logger.info(f"VEP index built: {len(self._index)} variants loaded")

    def _load_chrom_file(self, fpath: Path):
        """Load a single chromosome VEP JSON.gz file."""
        chrom = fpath.stem.removeprefix("VEP_").removesuffix(".json")
        self._loaded_chroms.add(chrom)

        try:
            opener = gzip.open if fpath.suffix == ".gz" else open
            with opener(fpath, "rt", encoding="utf-8") as f:
                for line in f:
                    if not line.strip():
                        continue
                    try:
                        entry = json.loads(line)
                    except json.JSONDecodeError:
                        continue

                    # Key: (chrom, pos, ref, alt)
                    key = self._make_key(entry)
                    if key:
                        self._index[key] = self._parse_entry(entry)
        except Exception as e:
            logger.error(f"Failed to load {fpath}: {e}")

    def _make_key(self, entry: dict) -> Optional[tuple]:
        chrom = str(entry.get("seq_region_name", "")).replace("chr", "")
        pos = entry.get("start")
        ref = entry.get("allele_string", "").split("/")[0]
        alt = entry.get("variant")
        if not all([chrom, pos, ref, alt]):
            return None
        return (chrom, int(pos), ref.upper(), alt.upper())

    def _parse_entry(self, entry: dict) -> dict:
        """Parse a VEP JSON entry into standardized fields."""
        transcripts = entry.get("transcript_consequences", [])
        gene_info = {}

        if transcripts:
            # Pick the canonical/best transcript
            best = transcripts[0]
            gene_info = {
                "gene_symbol": best.get("gene_symbol", ""),
                "gene_id": best.get("gene_id", ""),
                "transcript_id": best.get("transcript_id", ""),
                "consequence": "|".join(best.get("consequences", [])),
                "consequence_category": self._best_category(
                    best.get("consequences", [])
                ),
                "protein_change": best.get("amino_acids", ""),  # e.g. "R/W"
                "codon_change": best.get("codons", ""),
                "sift_score": best.get("sift_score"),
                "sift_prediction": best.get("sift_pred"),
                "polyphen_score": best.get("polyphen_score"),
                "polyphen_prediction": best.get("polyphen_pred"),
                "canonical": best.get("canonical"),
                "biotype": best.get("biotype", ""),
                "exon": best.get("exon", ""),
                "intron": best.get("intron", ""),
                "uniprot_isoform": best.get("uniparc", ""),
            }

        # Variant-level info
        intergenic = entry.get("intergenic_consequences", [])
        if intergenic and not gene_info:
            gene_info = {
                "consequence": "intergenic_variant",
                "consequence_category": "LOW",
            }

        # Colocated variants
        colocalated = entry.get("colocated_variant", {})
        gene_info.update({
            "gnomad_af": colocalated.get("gnomad_AF"),
            "gnomad_afr_af": colocalated.get("gnomad_AF_AFR"),
            "gnomad_amr_af": colocalated.get("gnomad_AF_AMR"),
            "gnomad_eas_af": colocalated.get("gnomad_AF_EAS"),
            "gnomad_nfe_af": colocalated.get("gnomad_AF_NFE"),
            "gnomad_sas_af": colocalated.get("gnomad_AF_SAS"),
            "clinvar_id": colocalated.get("clinvar"),
            "dbSNP_id": colocalated.get("id"),
            "vep_id": entry.get("id", ""),
        })

        # Population frequencies
        frequencies = entry.get("frequencies", {})
        for freq_key, freq_val in frequencies.items():
            gene_info[f"freq_{freq_key}"] = freq_val

        return gene_info

    def _best_category(self, consequences: list[str]) -> str:
        """Return highest impact category from consequence list."""
        for csq in consequences:
            cat = self.CONSEQUENCE_CATEGORIES.get(csq)
            if cat == "HIGH":
                return "HIGH"
        for csq in consequences:
            cat = self.CONSEQUENCE_CATEGORIES.get(csq)
            if cat == "MODERATE":
                return "MODERATE"
        for csq in consequences:
            cat = self.CONSEQUENCE_CATEGORIES.get(csq)
            if cat == "LOW":
                return "LOW"
        return "MODIFIER"

    def annotate_variants(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Add VEP annotations to a variants DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Must have columns: CHROM (or chrom), POS (or pos), REF, ALT

        Returns
        -------
        pd.DataFrame
            Original df with new VEP annotation columns
        """
        if self._index:
            return self._annotate_from_index(df)
        else:
            logger.warning("VEP index not built. Call build_index() first for efficient annotation.")
            return self._annotate_lazy(df)

    def _annotate_from_index(self, df: pd.DataFrame) -> pd.DataFrame:
        """Annotate using pre-built index."""
        df = df.copy()
        chrom_col = "CHROM" if "CHROM" in df.columns else "chrom"
        pos_col = "POS" if "POS" in df.columns else "pos"

        # Initialize annotation columns
        for col in ["gene_symbol", "consequence", "consequence_category",
                    "sift_prediction", "polyphen_prediction", "gnomad_af",
                    "clinvar_id", "dbSNP_id"]:
            if col not in df.columns:
                df[col] = None

        matched = 0
        for idx, row in df.iterrows():
            key = (str(row[chrom_col]).replace("chr", ""),
                   int(row[pos_col]),
                   str(row.get("REF", "")).upper(),
                   str(row.get("ALT", "")).upper())
            ann = self._index.get(key)
            if ann:
                for k, v in ann.items():
                    if v is not None and k in df.columns:
                        df.at[idx, k] = v
                    elif v is not None:
                        df[k] = v
                matched += 1

        logger.info(f"VEP annotation: {matched}/{len(df)} variants matched")
        return df

    def _annotate_lazy(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Lazy mode: load chromosome file on demand for each variant.
        Much slower — use for one-off queries only.
        """
        chrom_col = "CHROM" if "CHROM" in df.columns else "chrom"
        pos_col = "POS" if "POS" in df.columns else "pos"

        # Group variants by chromosome
        df["_chrom"] = df[chrom_col].astype(str).str.replace("chr", "")

        chrom_groups = df.groupby("_chrom")
        results = []

        for chrom, group in chrom_groups:
            fpath = self.cache_dir / f"VEP_{chrom}.json.gz"
            if not fpath.exists():
                results.append(group)
                continue

            # Build temporary index for this chromosome
            chrom_index = {}
            opener = gzip.open if fpath.suffix == ".gz" else open
            with opener(fpath, "rt", encoding="utf-8") as f:
                for line in f:
                    if not line.strip():
                        continue
                    try:
                        entry = json.loads(line)
                        key = self._make_key(entry)
                        if key:
                            chrom_index[key] = self._parse_entry(entry)
                    except Exception:
                        continue

            group = group.copy()
            for idx, row in group.iterrows():
                key = (chrom, int(row[pos_col]),
                       str(row.get("REF", "")).upper(),
                       str(row.get("ALT", "")).upper())
                ann = chrom_index.get(key, {})
                for k, v in ann.items():
                    if v is not None:
                        group.at[idx, k] = v
            results.append(group)

        df = pd.concat(results, ignore_index=True)
        df.drop(columns=["_chrom"], inplace=True, errors="ignore")
        return df

    @staticmethod
    def from_vep_vcf(vep_vcf_path: Union[str, Path]) -> pd.DataFrame:
        """
        Parse VEP output in VCF format (--vcf flag).

        Returns DataFrame with parsed CSQ fields.
        """
        vep_vcf_path = Path(vep_vcf_path)
        csq_data = []
        csq_fields = []

        with gzip.open(vep_vcf_path, "rt") as f:
            for line in f:
                if line.startswith("##INFO=<ID=CSQ"):
                    # Extract CSQ field names from header
                    m = re.search(r'Format: "([^"]+)"', line)
                    if m:
                        csq_fields = m.group(1).split("|")
                elif line.startswith("#CHROM"):
                    continue
                elif "CSQ=" in line:
                    # Parse INFO field
                    info = dict(p.split("=", 1) for p in line.split("\t")[7].split(";") if "=" in p)
                    csq_str = info.get("CSQ", "")
                    for alt_allele in csq_str.split(","):
                        values = alt_allele.split("|")
                        if len(values) == len(csq_fields):
                            csq_data.append(dict(zip(csq_fields, values)))

        df = pd.DataFrame(csq_data)
        logger.info(f"Parsed {len(df)} VEP annotations from VCF")
        return df
