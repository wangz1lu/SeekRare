"""
clinvar_annotation.py — ClinVar 注释

输入:
    - 2_gtf_annotated.csv: GTF 注释后的变异列表（含 gene_name）
    - clinvar CSV: CHROM, POS, REF, ALT, CLNDISDB, CLNDN, CLNREVSTAT, CLNSIG, CLNVC, ORIGIN

输出:
    - 第一列: CPRA = CHROM:POS:REF:ALT
    - CHROM, POS, REF, ALT, in_gene, gene_name, feature_type, D2007018873, ...
    - CLNDISDB, CLNDN, CLNREVSTAT, CLNSIG, CLNVC, ORIGIN（匹配不上留空）

匹配: CHROM + POS + REF + ALT 四键精确匹配（全字符串，无 dtype 问题）
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


# CLNSIG 数值→文字映射
CLNSIG_MAP = {
    "0": "not provided",
    "1": "Pathogenic",
    "2": "Likely pathogenic",
    "3": "Likely benign",
    "4": "Benign",
    "5": "Uncertain significance",
    "6": "Conflicting interpretations",
    "255": "other",
    "247": "Pathogenic/Likely pathogenic",
    "248": "Likely pathogenic/Pathogenic",
    "249": "Likely benign/Benign",
    "250": "not provided",
    "251": "not provided",
}


def _build_cpra(df: pd.DataFrame) -> pd.Series:
    """
    向量化构建 CPRA 列: CHROM:POS:REF:ALT
    CHROM 统一去 chr 前缀，REF/ALT 转大写
    """
    return (
        df["CHROM"].str.lstrip("chr").astype(str)
        + ":"
        + df["POS"].astype(str)
        + ":"
        + df["REF"].str.upper()
        + ":"
        + df["ALT"].str.upper()
    )


def merge_filter_clinvar(
    input_csv: Union[str, Path],
    clinvar_csv: Union[str, Path],
    output_csv: Union[str, Path],
) -> pd.DataFrame:
    """
    将 ClinVar 注释合入 GTF 注释后的变异 CSV。

    参数
    ----
    input_csv : str   2_gtf_annotated.csv 路径
    clinvar_csv : str  ClinVar CSV 路径
    output_csv : str   输出路径

    返回
    ----
    pd.DataFrame
        第一列 CPRA，其余列按顺序:
        CHROM, POS, REF, ALT, in_gene, gene_name, feature_type, ...,
        D2007018873, CLNDISDB, CLNDN, CLNREVSTAT, CLNSIG, CLNVC, ORIGIN
    """
    input_csv = str(input_csv)
    clinvar_csv = str(clinvar_csv)
    output_csv = str(output_csv)

    # ── 1. 加载 GTF CSV ─────────────────────────────────────────────────────
    logger.info(f"Loading GTF CSV: {input_csv}")
    df = pd.read_csv(input_csv, dtype=str)
    logger.info(f"  Variants: {len(df)}")

    # ── 2. 向量化构建 CPRA ──────────────────────────────────────────────────
    logger.info("Building CPRA keys (vectorized)...")
    df["CPRA"] = _build_cpra(df)
    logger.info(f"  CPRA built, sample: {df['CPRA'].iloc[0]}")

    # ── 3. 加载 ClinVar ─────────────────────────────────────────────────────
    logger.info(f"Loading ClinVar: {clinvar_csv}")

    cv = pd.read_csv(clinvar_csv, dtype=str, low_memory=False)
    cv.columns = cv.columns.str.strip()

    if "#CHROM" in cv.columns:
        cv.rename(columns={"#CHROM": "CHROM"}, inplace=True)

    clinvar_cols = ["CLNDISDB", "CLNDN", "CLNREVSTAT", "CLNSIG", "CLNVC", "ORIGIN"]
    for col in clinvar_cols:
        if col not in cv.columns:
            cv[col] = ""

    cv = cv[["CHROM", "POS", "REF", "ALT"] + clinvar_cols].copy()

    # 归一化 CHROM（去 chr 前缀，chrM→MT）
    cv["CHROM"] = cv["CHROM"].str.lstrip("chr").str.replace("M", "MT", regex=False)

    # 构建 ClinVar CPRA
    cv["CPRA"] = (
        cv["CHROM"].astype(str)
        + ":"
        + cv["POS"].astype(str)
        + ":"
        + cv["REF"].str.upper()
        + ":"
        + cv["ALT"].str.upper()
    )

    # CLNSIG 数字→文字
    cv["CLNSIG"] = cv["CLNSIG"].map(CLNSIG_MAP).fillna(cv["CLNSIG"])

    # 去重：同 CPRA 保留第一个
    cv = cv.drop_duplicates(subset=["CPRA"], keep="first")

    # ── 4. pandas merge（CPRA 对 CPRA，字符串对字符串）──────────────────────────
    logger.info("Merging on CPRA...")
    n_before = len(df)
    df = df.merge(
        cv[["CPRA"] + clinvar_cols],
        on="CPRA",
        how="left",
    )
    n_matched = df["CLNSIG"].notna().sum()
    logger.info(f"  Matched: {n_matched}/{n_before} variants with ClinVar records")

    # ── 5. 整理列顺序 ────────────────────────────────────────────────────────
    # 优先级顺序
    first_cols = ["CPRA", "CHROM", "POS", "REF", "ALT"]
    # 原 GTF 的核心注释列
    gtf_cols = ["in_gene", "gene_name", "feature_type"]
    # ClinVar 列放最后
    last_cols = clinvar_cols

    # 收集已有列
    ordered = []
    for c in first_cols:
        if c in df.columns:
            ordered.append(c)
    for c in df.columns:
        if c not in first_cols and c not in last_cols:
            ordered.append(c)
    for c in last_cols:
        if c in df.columns:
            ordered.append(c)

    df = df[ordered]

    # ── 6. 保存 ──────────────────────────────────────────────────────────────
    df.to_csv(output_csv, index=False)
    logger.info(f"ClinVar annotation saved → {output_csv}")
    logger.info(f"  Columns: {df.columns.tolist()}")

    return df
