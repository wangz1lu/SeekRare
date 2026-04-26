"""
VCF Parser — supports trio (proband + father + mother) VCF files.

Merges family VCFs and infers inheritance mode annotations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import pandas as pd
import pysam
from loguru import logger


def load_trio_vcf(
    proband_vcf: Union[str, Path],
    father_vcf: Optional[Union[str, Path]] = None,
    mother_vcf: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Load and merge VCF files for a family trio.

    Adds columns:
        - sample: sample name
        - father_gt: father's genotype (if provided)
        - mother_gt: mother's genotype (if provided)
        - inheritance: inferred inheritance pattern (AD/AR/XR/XD/denovo/unknown)

    Parameters
    ----------
    proband_vcf : str or Path
        Path to proband VCF
    father_vcf : str or Path, optional
        Path to father's VCF
    mother_vcf : str or Path, optional
        Path to mother's VCF

    Returns
    -------
    pd.DataFrame
        One row per variant, with genotype and inheritance annotations.
    """
    records = []

    # ── Parse proband VCF ──────────────────────────────────────
    proband_records = _parse_vcf(proband_vcf, sample_name="proband")
    records.extend(proband_records)
    logger.info(f"  Proband: {len(proband_records)} variants from {proband_vcf}")

    # ── Parse father VCF ──────────────────────────────────────
    father_records = {}
    if father_vcf:
        father_records_list = _parse_vcf(father_vcf, sample_name="father")
        father_records = {f"{r['chrom']}-{r['pos']}-{r['ref']}-{r['alt']}": r
                          for r in father_records_list}
        logger.info(f"  Father: {len(father_records)} variants from {father_vcf}")

    # ── Parse mother VCF ─────────────────────────────────────
    mother_records = {}
    if mother_vcf:
        mother_records_list = _parse_vcf(mother_vcf, sample_name="mother")
        mother_records = {f"{r['chrom']}-{r['pos']}-{r['ref']}-{r['alt']}": r
                         for r in mother_records_list}
        logger.info(f"  Mother: {len(mother_records)} variants from {mother_vcf}")

    # ── Merge and annotate inheritance ──────────────────────
    df = pd.DataFrame(records)

    if father_vcf:
        df["father_gt"] = df.apply(
            lambda r: father_records.get(f"{r['chrom']}-{r['pos']}-{r['ref']}-{r['alt']}", {}).get("gt", "./."),
            axis=1
        )
    if mother_vcf:
        df["mother_gt"] = df.apply(
            lambda r: mother_records.get(f"{r['chrom']}-{r['pos']}-{r['ref']}-{r['alt']}", {}).get("gt", "./."),
            axis=1
        )

    # ── Infer inheritance ────────────────────────────────────
    df = _annotate_inheritance(df)

    return df.reset_index(drop=True)


def _parse_vcf(vcf_path: Union[str, Path], sample_name: str) -> list[dict]:
    """Parse a single VCF file into a list of record dicts."""
    records = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        for i, rec in enumerate(vf):
            for sample in rec.samples:
                gt = rec.samples[sample].get("GT", None)
                ad = rec.samples[sample].get("AD", None)
                dp = rec.samples[sample].get("DP", None)
                gq = rec.samples[sample].get("GQ", None)

                records.append({
                    "chrom": str(rec.chrom),
                    "pos": int(rec.pos),
                    "rsid": rec.id if rec.id else ".",
                    "ref": str(rec.ref),
                    "alt": ",".join(str(a) for a in rec.alts) if rec.alts else ".",
                    "qual": float(rec.qual) if rec.qual else ".",
                    "filter": ";".join(rec.filter) if rec.filter else "PASS",
                    "gene": _extract_gene_ann(rec),
                    "sample": sample_name,
                    "gt": "/".join(str(g) if g is not None else "." for g in gt) if gt else "./.",
                    "ad": ",".join(str(a) for a in ad) if ad else ".",
                    "dp": dp if dp is not None else ".",
                    "gq": gq if gq is not None else ".",
                })
    return records


def _extract_gene_ann(rec) -> str:
    """Extract gene symbol from VCF INFO field."""
    # Common gene annotations in INFO field
    for field in ["GENE", "ANN", "CSQ"]:
        if field in rec.info:
            val = rec.info[field]
            if isinstance(val, (list, tuple)):
                val = val[0]
            # Try to extract gene from ANN/CSQ field
            if field in ["ANN", "CSQ"]:
                parts = str(val).split("|")
                if len(parts) > 3:
                    return parts[3] if parts[3] != "" else "."
            return str(val)
    return "."


def _annotate_inheritance(df: pd.DataFrame) -> pd.DataFrame:
    """
    Infer inheritance pattern from trio genotypes.

    Genotype format: "0/1", "1/1", "0/0", "./."
    - 0 = ref allele
    - 1 = alt allele
    """

    def _gt_to_tuple(gt_str: str) -> tuple:
        if pd.isna(gt_str) or gt_str == "./.":
            return (None, None)
        try:
            return tuple(int(x) for x in gt_str.split("/"))
        except Exception:
            return (None, None)

    def _infer(row):
        p_gt = _gt_to_tuple(row.get("gt", "./."))
        f_gt = _gt_to_tuple(row.get("father_gt", "./.")) if "father_gt" in row else (None, None)
        m_gt = _gt_to_tuple(row.get("mother_gt", "./.")) if "mother_gt" in row else (None, None)

        # De novo: proband has alt, both parents are ref (or missing)
        if p_gt[0] is not None and p_gt[0] > 0:
            p_has_alt = any(a > 0 for a in p_gt if a is not None)
            if p_has_alt:
                f_has_alt = any(a > 0 for a in f_gt if a is not None)
                m_has_alt = any(a > 0 for a in m_gt if a is not None)
                if not f_has_alt and not m_has_alt:
                    return "denovo"

        # Autosomal recessive: both parents heterozygous, proband homozygous alt
        if (f_gt[0] is not None and m_gt[0] is not None and
            p_gt[0] is not None and p_gt[1] is not None):
            f_het = (f_gt[0] != f_gt[1]) and (f_gt[0] > 0 or f_gt[1] > 0)
            m_het = (m_gt[0] != m_gt[1]) and (m_gt[0] > 0 or m_gt[1] > 0)
            if f_het and m_het and p_gt[0] > 0 and p_gt[1] > 0:
                return "AR_comp"
            # AR homozygous
            if not f_has_alt and not m_has_alt if 'f_has_alt' in dir() else False:
                pass

        return "unknown"

    if "father_gt" not in df.columns:
        df["inheritance"] = "unknown"
    else:
        df["inheritance"] = df.apply(_infer, axis=1)

    return df
