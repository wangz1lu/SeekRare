"""
stage3_prep.py — Stage 3 数据预处理

功能:
1. 解析 CLNDISDB / CLNDN / HPO / OMIM / Orphanet 复杂标签
2. 统计各动态列 unique 值分布（gene_name, HPO, OMIM, Orphanet, inheritance_mode, MC）
3. 生成 LLM prompt 摘要
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


# ── 标签解析 ─────────────────────────────────────────────────────────────────

def parse_clndisdb_tags(value: str) -> list[str]:
    """解析 CLNDISDB，只保留 ontology 前缀标签（OMIM: / Orphanet: / HP: / MedGen:）。"""
    if pd.isna(value) or str(value).strip() == "":
        return []
    tags = []
    groups = str(value).strip().split("|")
    for group in groups:
        for part in group.split(","):
            part = part.strip()
            if not part or ":" not in part:
                continue
            tags.append(part)
    seen = set()
    unique = []
    for t in tags:
        if t not in seen:
            seen.add(t)
            unique.append(t)
    return unique


def parse_clndn_tags(value: str) -> list[str]:
    """解析 CLNDN，包含 ontology 标签 + 疾病名称。"""
    if pd.isna(value) or str(value).strip() == "":
        return []
    tags = []
    groups = str(value).strip().split("|")
    for group in groups:
        for part in group.split(","):
            part = part.strip()
            if not part:
                continue
            if ":" in part:
                tags.append(part)
            else:
                tags.append(part.replace("_", " ").lower())
    seen = set()
    unique = []
    for t in tags:
        if t not in seen:
            seen.add(t)
            unique.append(t)
    return unique


def parse_hpo_tags(value: str) -> list[str]:
    """解析 HPO 列（HP:xxxx 格式，可能用 ; 分隔）"""
    if pd.isna(value) or str(value).strip() == "":
        return []
    return [t.strip() for t in str(value).strip().split(";") if t.strip() and t.strip().startswith("HP:")]


def parse_omim_tags(value: str) -> list[str]:
    """解析 OMIM 列（OMIM:xxxxx 格式，可能用 ; 分隔）"""
    if pd.isna(value) or str(value).strip() == "":
        return []
    return [t.strip() for t in str(value).strip().split(";") if t.strip() and t.strip().startswith("OMIM:")]


def parse_orphanet_tags(value: str) -> list[str]:
    """解析 Orphanet 列（Orphanet:xxxxx 格式，可能用 ; 分隔）"""
    if pd.isna(value) or str(value).strip() == "":
        return []
    return [t.strip() for t in str(value).strip().split(";") if t.strip() and t.strip().startswith("Orphanet:")]


# ── 统计 ─────────────────────────────────────────────────────────────────────

def summarize_stage1(csv_path: str, sample_limit: int = 5_000_000) -> dict:
    """
    读取 Stage 1/2 CSV，生成 LLM 摘要（兼容新旧列结构）。

    Returns dict:
        n_total, n_with_clinvar,
        gene_name_counts, HPO_tag_counts, OMIM_tag_counts, Orphanet_tag_counts,
        inheritance_mode_counts, MC_counts,
        feature_type_counts, significance_counts (from old CLNSIG), clinvarstar_counts,
        eqtl_tissue_counts, splicevardb_counts,
        sample_rows
    """
    logger.info(f"Loading Stage CSV: {csv_path}")
    df = pd.read_csv(csv_path, dtype=str, low_memory=False)
    n_total = len(df)
    logger.info(f"  Total rows: {n_total:,}")

    result = {"n_total": int(n_total)}

    # ── gene_name ─────────────────────────────────────────────────────────
    result["gene_name_counts"] = df["gene_name"].value_counts().to_dict()
    result["diseasename_counts"] = df.get("diseasename", pd.Series(dtype=str)).value_counts().to_dict()

    # ── inheritance_mode ──────────────────────────────────────────────────
    result["inheritance_mode_counts"] = df["inheritance_mode"].value_counts().to_dict()

    # ── HPO ───────────────────────────────────────────────────────────────
    hpo_counter = Counter()
    for val in df.get("HPO", pd.Series(dtype=str)).dropna():
        for tag in parse_hpo_tags(val):
            hpo_counter[tag] += 1
    result["HPO_tag_counts"] = dict(hpo_counter.most_common(200))
    logger.info(f"  HPO unique tags: {len(hpo_counter)}")

    # ── OMIM ──────────────────────────────────────────────────────────────
    omim_counter = Counter()
    for val in df.get("OMIM", pd.Series(dtype=str)).dropna():
        for tag in parse_omim_tags(val):
            omim_counter[tag] += 1
    result["OMIM_tag_counts"] = dict(omim_counter.most_common(200))
    logger.info(f"  OMIM unique tags: {len(omim_counter)}")

    # ── Orphanet ───────────────────────────────────────────────────────────
    orphanet_counter = Counter()
    for val in df.get("Orphanet", pd.Series(dtype=str)).dropna():
        for tag in parse_orphanet_tags(val):
            orphanet_counter[tag] += 1
    result["Orphanet_tag_counts"] = dict(orphanet_counter.most_common(200))
    logger.info(f"  Orphanet unique tags: {len(orphanet_counter)}")

    # ── MC ────────────────────────────────────────────────────────────────
    result["MC_counts"] = df.get("MC", pd.Series(dtype=str)).value_counts().to_dict()
    logger.info(f"  MC unique values: {len(result['MC_counts'])}")

    # ── feature_type ──────────────────────────────────────────────────────
    result["feature_type_counts"] = df["feature_type"].value_counts().to_dict()

    # ── significance（从 CLNSIG 继承，简化 CSV 后改名）──────────────────────
    # 简化 CSV 后，significance 列存在，CLNSIG 列已删除
    if "significance" in df.columns:
        result["significance_counts"] = df["significance"].value_counts().to_dict()
    elif "CLNSIG" in df.columns:
        result["significance_counts"] = df["CLNSIG"].value_counts().to_dict()
    else:
        result["significance_counts"] = {}

    # ── clinvarstar ────────────────────────────────────────────────────────
    if "clinvarstar" in df.columns:
        result["clinvarstar_counts"] = df["clinvarstar"].value_counts().to_dict()

    # ── eqtl_tissue ────────────────────────────────────────────────────────
    if "eqtl_tissue" in df.columns:
        result["eqtl_tissue_counts"] = df["eqtl_tissue"].value_counts().to_dict()

    # ── splicevardb ────────────────────────────────────────────────────────
    if "splicevardb" in df.columns:
        result["splicevardb_counts"] = df["splicevardb"].value_counts().to_dict()
        logger.info(f"  splicevardb values: {len(result['splicevardb_counts'])}")

    # ── n_with_clinvar ────────────────────────────────────────────────────
    if "significance" in df.columns:
        result["n_with_clinvar"] = df["significance"].notna().sum()
    elif "CLNSIG" in df.columns:
        result["n_with_clinvar"] = df["CLNSIG"].notna().sum()
    else:
        result["n_with_clinvar"] = 0

    # ── sample rows ────────────────────────────────────────────────────────
    sample_df = df if n_total <= sample_limit else df.sample(n=sample_limit, random_state=42)
    cols_show = ["CPRA", "gene_name", "inheritance_mode", "feature_type",
                 "significance", "clinvarstar", "HPO", "OMIM", "Orphanet",
                 "MC", "eqtl_tissue", "splicevardb"]
    cols_show = [c for c in cols_show if c in df.columns]
    result["sample_rows"] = sample_df.head(5)[cols_show].fillna("").to_dict("records")

    return result


def build_llm_prompt(symptoms: str, summary: dict) -> str:
    """兼容旧接口，指向新版 prompt builder。"""
    from seekrare.scoring.stage3_annotate import build_llm_prompt_new
    return build_llm_prompt_new(symptoms, "unknown", summary)