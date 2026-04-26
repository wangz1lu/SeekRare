"""
Annotation Combiner — merges all annotation sources into a single DataFrame.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


def annotate_variants(
    variants_df: pd.DataFrame,
    clinvar_vcf: Optional[Union[str, Path]] = None,
    vep_cache_dir: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Add annotation columns to a variants DataFrame.

    Parameters
    ----------
    variants_df : pd.DataFrame
        Output from load_trio_vcf()
    clinvar_vcf : str or Path, optional
        Path to ClinVar VCF
    vep_cache_dir : str or Path, optional
        Path to VEP cache directory

    Returns
    -------
    pd.DataFrame
        variants_df with new annotation columns added
    """
    df = variants_df.copy()
    n_before = len(df)

    if clinvar_vcf:
        df = _annotate_clinvar(df, clinvar_vcf)

    if vep_cache_dir:
        df = _annotate_vep(df, vep_cache_dir)

    # ── Add placeholder columns if not present ─────────────────
    default_cols = {
        "cadd_score": 0.0,
        "gnomad_af": 0.0,
        "sift_score": 0.5,
        "polyphen_score": 0.5,
        "impact": "UNKNOWN",
        "hpo_terms": ".",
    }
    for col, default in default_cols.items():
        if col not in df.columns:
            df[col] = default

    n_after = len(df)
    logger.info(f"  Annotation complete: {n_before} variants, columns: {list(df.columns)}")
    return df.reset_index(drop=True)


def _annotate_clinvar(df: pd.DataFrame, clinvar_vcf: Union[str, Path]) -> pd.DataFrame:
    """Annotate variants with ClinVar data."""
    import pysam

    try:
        clinvar = pysam.VariantFile(str(clinvar_vcf))
    except Exception as e:
        logger.warning(f"Could not open ClinVar VCF: {e}")
        df["clinvar_significance"] = "not provided"
        df["clinvar_stars"] = 0
        return df

    # Build ClinVar lookup: (chrom, pos, ref, alt) → (significance, stars)
    clinvar_map = {}
    for rec in clinvar:
        chrom = str(rec.chrom).replace("chr", "")
        for alt in rec.alts:
            key = (chrom, int(rec.pos), str(rec.ref), str(alt))
            sig = rec.info.get("CLNSIG", ["not provided"])[0]
            try:
                stars = int(rec.info.get("CLNREVSTAT", ["none"])[0])
            except Exception:
                stars = 0
            clinvar_map[key] = (str(sig), stars)

    clinvar.close()

    def lookup(row):
        key = (str(row["chrom"]).replace("chr", ""), int(row["pos"]), str(row["ref"]), str(row["alt"]))
        return clinvar_map.get(key, ("not provided", 0))

    results = df.apply(lookup, axis=1, result_type="expand")
    df["clinvar_significance"] = results[0]
    df["clinvar_stars"] = results[1]

    n_annotated = (df["clinvar_significance"] != "not provided").sum()
    logger.info(f"  ClinVar: {n_annotated}/{len(df)} variants annotated")
    return df


def _annotate_vep(df: pd.DataFrame, vep_cache_dir: Union[str, Path]) -> pd.DataFrame:
    """Annotate with VEP cache."""
    # Placeholder: would use pyensys or ANNOVAR-style lookup
    logger.info(f"  VEP annotation not yet implemented (cache: {vep_cache_dir})")
    return df
