"""
clinvar_annotation.py — ClinVar 注释（简化版）

输入:
    - 2_gtf_annotated.csv: GTF 注释后的变异列表（含 gene_name）
    - clinvar CSV: CHROM, POS, ID, REF, ALT, CLNDISDB, CLNDN, CLNREVSTAT, CLNSIG, CLNVC, ORIGIN

输出:
    - 第一列: CPRA = CHROM:POS:REF:ALT（合并主键）
    - 其余列: CHROM, POS, REF, ALT, in_gene, gene_name, feature_type,
              [原 GTF 列], D2007018873, CLNDISDB, CLNDN, CLNREVSTAT, CLNSIG, CLNVC, ORIGIN
    - 匹配不上的 ClinVar 列留空

匹配: CHROM + POS + REF + ALT 四键精确匹配（全部转字符串后匹配）
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


def _build_cpra(row_or_dict) -> str:
    """Build CPRA key: CHROM:POS:REF:ALT (all as strings)."""
    chrom = str(row_or_dict.get("CHROM", "")).strip().lstrip("chr")
    pos = str(int(float(row_or_dict.get("POS", 0))))
    ref = str(row_or_dict.get("REF", "")).strip().upper()
    alt = str(row_or_dict.get("ALT", "")).strip().upper()
    return f"{chrom}:{pos}:{ref}:{alt}"


def _load_clinvar(clinvar_path: str) -> pd.DataFrame:
    """
    加载 ClinVar CSV，构建 CPRA→ClinVar列 的字典。
    返回 dict: {cpra_key → {col: value, ...}}
    """
    logger.info(f"Loading ClinVar: {clinvar_path}")

    df = pd.read_csv(clinvar_path, dtype=str, low_memory=False)
    df.columns = df.columns.str.strip()

    # 处理 #CHROM
    if "#CHROM" in df.columns:
        df.rename(columns={"#CHROM": "CHROM"}, inplace=True)

    clinvar_cols = ["CLNDISDB", "CLNDN", "CLNREVSTAT", "CLNSIG", "CLNVC", "ORIGIN"]
    for col in clinvar_cols:
        if col not in df.columns:
            df[col] = ""

    # 只保留需要的列
    df = df[["CHROM", "POS", "REF", "ALT"] + clinvar_cols].copy()

    # 归一化 CHROM
    df["CHROM"] = df["CHROM"].apply(
        lambda x: str(x).strip().lstrip("chr").replace("M", "MT")
    )

    # 构建 CPRA key
    df["_cpra"] = df.apply(_build_cpra, axis=1)

    # 转 CLNSIG 数字→文字
    df["CLNSIG"] = df["CLNSIG"].map(CLNSIG_MAP).fillna(df["CLNSIG"])

    # 去重：同 CPRA 只保留第一个
    df = df.drop_duplicates(subset=["_cpra"], keep="first")

    # 转 dict: cpra → {CLNDISDB, CLNDN, ...}
    records = df.set_index("_cpra")[clinvar_cols].to_dict("index")

    logger.info(f"ClinVar loaded: {len(records)} unique CPRA keys")
    return records


def merge_filter_clinvar(
    input_csv: Union[str, Path],
    clinvar_csv: Union[str, Path],
    output_csv: Union[str, Path],
) -> pd.DataFrame:
    """
    将 ClinVar 注释合入 GTF 注释后的变异 CSV。

    参数
    ----
    input_csv : str
        2_gtf_annotated.csv 路径
    clinvar_csv : str
        ClinVar CSV 路径
    output_csv : str
        输出路径

    返回
    ----
    pd.DataFrame
        第一列为 CPRA (CHROM:POS:REF:ALT)，其余列依次为：
        CHROM, POS, REF, ALT, [原 GTF 列],
        D2007018873, CLNDISDB, CLNDN, CLNREVSTAT, CLNSIG, CLNVC, ORIGIN
    """
    input_csv = str(input_csv)
    clinvar_csv = str(clinvar_csv)
    output_csv = str(output_csv)

    # ── 1. 加载 GTF CSV ──────────────────────────────────────────────────────
    logger.info(f"Loading GTF CSV: {input_csv}")
    df = pd.read_csv(input_csv, dtype=str)
    logger.info(f"  Variants: {len(df)}")

    # ── 2. 构建 CPRA 列 ─────────────────────────────────────────────────────
    logger.info("Building CPRA keys...")
    df["CPRA"] = df.apply(_build_cpra, axis=1)

    # ── 3. 加载 ClinVar ─────────────────────────────────────────────────────
    clin_dict = _load_clinvar(clinvar_csv)

    # ── 4. 用 CPRA 字典填充 ClinVar 列 ─────────────────────────────────────
    clinvar_cols = ["CLNDISDB", "CLNDN", "CLNREVSTAT", "CLNSIG", "CLNVC", "ORIGIN"]
    for col in clinvar_cols:
        df[col] = ""   # 先填空

    n_matched = 0
    for idx, cpra in enumerate(df["CPRA"]):
        if cpra in clin_dict:
            n_matched += 1
            for col in clinvar_cols:
                df.at[idx, col] = clin_dict[cpra].get(col, "")

        if idx > 0 and idx % 5_000_000 == 0:
            logger.info(f"  Processed {idx:,} / {len(df):,} variants ({n_matched:,} matched)")

    logger.info(f"  Matched: {n_matched}/{len(df)} variants with ClinVar records")

    # ── 5. 调整列顺序 ────────────────────────────────────────────────────────
    # CPRA 作为第一列，其余按: CHROM, POS, REF, ALT, 原GTF列, GT列, ClinVar列
    priority_first = ["CPRA", "CHROM", "POS", "REF", "ALT"]
    priority_last = clinvar_cols   # CLNDISDB CLNDN CLNREVSTAT CLNSIG CLNVC ORIGIN

    cols = [c for c in df.columns if c in priority_first]
    rest = [c for c in df.columns if c not in priority_first and c not in priority_last]
    last = [c for c in df.columns if c in priority_last]

    df = df[cols + rest + last]

    # ── 6. 保存 ──────────────────────────────────────────────────────────────
    df.to_csv(output_csv, index=False)
    logger.info(f"ClinVar annotation saved → {output_csv}")
    logger.info(f"  Columns: {df.columns.tolist()}")

    return df
