"""
clinvar_annotation.py — ClinVar 注释（简化版）

输入:
    - 2_gtf_annotated.csv: GTF 注释后的变异列表（含 gene_name）
      列: CHROM, POS, REF, ALT, gene_name, in_gene, feature_type, ...
    - clinvar CSV: ClinVar 注释库
      列: CHROM, POS, ID, REF, ALT, CLNDISDB, CLNDN, CLNREVSTAT, CLNSIG, CLNVC, ORIGIN

输出:
    - 继承 2_gtf_annotated.csv 的所有列，
    - 新增 CLNDISDB, CLNDN, CLNREVSTAT, CLNSIG, CLNVC, ORIGIN 列
    - 匹配不上的位点留空

匹配方式: CHROM + POS + REF + ALT 四键精确匹配（区分大小写）
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


# ClinVar CLNSIG 数值→文字映射（常见值）
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


def _normalize_chrom(chrom: str) -> str:
    """统一染色体格式: chr1→1, chrM→MT, 1→1"""
    chrom = str(chrom).strip().lstrip("chr")
    return chrom if chrom != "M" else "MT"


def _load_clinvar(clinvar_path: str) -> pd.DataFrame:
    """
    加载 ClinVar CSV，只取需要的列。

    Parameters
    ----------
    clinvar_path : str
        ClinVar CSV 路径

    Returns
    -------
    pd.DataFrame
        列: CHROM, POS, REF, ALT, CLNDISDB, CLNDN, CLNREVSTAT, CLNSIG, CLNVC, ORIGIN
        CHROM 已归一化，index 为 (CHROM, POS, REF, ALT) 四元组
    """
    needed = [
        "CHROM", "POS", "REF", "ALT",
        "CLNDISDB", "CLNDN", "CLNREVSTAT", "CLNSIG", "CLNVC", "ORIGIN",
    ]

    logger.info(f"Loading ClinVar: {clinvar_path}")

    # 检测是否有 CHROM 列（有的 CSV 用 #CHROM）
    df = pd.read_csv(clinvar_path, dtype=str, low_memory=False)
    df.columns = df.columns.str.strip()

    # 处理 #CHROM 列名
    if "#CHROM" in df.columns:
        df.rename(columns={"#CHROM": "CHROM"}, inplace=True)

    # 只保留需要的列（缺失的列填空）
    for col in needed:
        if col not in df.columns:
            df[col] = ""

    df = df[needed].copy()

    # 归一化 CHROM
    df["CHROM"] = df["CHROM"].apply(_normalize_chrom)

    # POS 转为整数
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce").fillna(0).astype("Int64")

    # 归一化 REF/ALT（大写，去除空白）
    df["REF"] = df["REF"].str.strip().str.upper()
    df["ALT"] = df["ALT"].str.strip().str.upper()

    # 建索引: (CHROM, POS, REF, ALT)
    df = df.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"])
    df = df.set_index(["CHROM", "POS", "REF", "ALT"])

    n = len(df)
    logger.info(f"ClinVar loaded: {n} entries, {df.index.nunique()} unique variants")
    return df


def merge_filter_clinvar(
    input_csv: Union[str, Path],
    clinvar_csv: Union[str, Path],
    output_csv: Union[str, Path],
) -> pd.DataFrame:
    """
    将 ClinVar 注释合入 GTF 注释后的变异 CSV。

    直接 CHROM + POS + REF + ALT 四键精确匹配，
    匹配不上则 CLNDISDB/CLNDN/CLNREVSTAT/CLNSIG/CLNVC/ORIGIN 全部留空。

    Parameters
    ----------
    input_csv : str
        2_gtf_annotated.csv 路径
    clinvar_csv : str
        ClinVar CSV 路径
    output_csv : str
        输出路径

    Returns
    -------
    pd.DataFrame
        合并后的 DataFrame
    """
    input_csv = str(input_csv)
    clinvar_csv = str(clinvar_csv)
    output_csv = str(output_csv)

    # ── 1. 加载 GTF 注释后的变异 CSV ──────────────────────────────────────────
    logger.info(f"Loading GTF-annotated CSV: {input_csv}")
    df = pd.read_csv(input_csv, dtype=str)
    logger.info(f"  Variants: {len(df)}")

    # 归一化 CHROM/POS/REF/ALT
    if "CHROM" not in df.columns:
        raise ValueError(f"input_csv 缺少 CHROM 列，列名: {df.columns.tolist()}")
    df["CHROM"] = df["CHROM"].apply(_normalize_chrom)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
    df["REF"] = df["REF"].str.strip().str.upper()
    df["ALT"] = df["ALT"].str.strip().str.upper()

    # ── 2. 加载 ClinVar ────────────────────────────────────────────────────────
    clin = _load_clinvar(clinvar_csv)

    # ── 3. 四键精确匹配 ────────────────────────────────────────────────────────
    logger.info("Matching variants to ClinVar (exact CHROM+POS+REF+ALT)...")

    n_before = len(df)
    matched = df.merge(
        clin,
        left_on=["CHROM", "POS", "REF", "ALT"],
        right_index=True,
        how="left",
    )
    df = matched

    # 统计
    n_clinvar = df["CLNSIG"].notna().sum()
    logger.info(f"  Matched: {n_clinvar}/{n_before} variants with ClinVar records")

    # ── 4. 清理 CLNSIG 数值→文字（如果需要）───────────────────────────────────
    # CLNSIG 如果是数字，转为文字
    def _map_clinsig(val):
        if pd.isna(val) or str(val).strip() == "":
            return val
        val_str = str(val).strip()
        return CLNSIG_MAP.get(val_str, val_str)

    if "CLNSIG" in df.columns:
        df["CLNSIG"] = df["CLNSIG"].apply(_map_clinsig)

    # ── 5. 保存 ────────────────────────────────────────────────────────────────
    df.to_csv(output_csv, index=False)
    logger.info(f"ClinVar annotation saved → {output_csv}")

    return df


def clinvar_filter_by_significance(
    input_csv: Union[str, Path],
    output_csv: Union[str, Path],
    sig_filter: Optional[list[str]] = None,
    exclude_benign: bool = True,
) -> pd.DataFrame:
    """
    按 ClinVar 致病性过滤变异（可选的后处理步骤）。

    Parameters
    ----------
    input_csv : str
        clinvar_annotation 后的 CSV
    output_csv : str
        输出路径
    sig_filter : list[str], optional
        只保留指定 CLNSIG 值，如 ["Pathogenic", "Likely pathogenic"]
    exclude_benign : bool
        若 True，排除 Benign / Likely benign

    Returns
    -------
    pd.DataFrame
        过滤后的 DataFrame
    """
    df = pd.read_csv(input_csv, dtype=str)
    n_before = len(df)

    if "CLNSIG" not in df.columns:
        logger.warning("CLNSIG column not found, skipping filter")
        return df

    if sig_filter:
        df = df[df["CLNSIG"].isin(sig_filter)]

    if exclude_benign:
        benign_terms = {"Benign", "Likely benign", "Likely benign/Benign"}
        df = df[~df["CLNSIG"].isin(benign_terms)]

    n_after = len(df)
    logger.info(f"ClinVar filter: {n_before} → {n_after} variants ({n_before - n_after} removed)")
    df.to_csv(output_csv, index=False)
    return df
