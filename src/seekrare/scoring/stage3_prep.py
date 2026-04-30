"""
stage3_prep.py — Stage 3 数据预处理

输入: Stage 1 输出的 CSV (3_clinvar_annotated.csv)
输出: LLM 需要的结构化统计摘要 + 解析后的 per-row 打分原料

功能:
1. 解析 CLNDISDB / CLNDN 的复杂标签（用 | 分组，用 , 分隔 ontology）
2. 统计 unique 标签集合
3. 统计每列值分布
4. 生成 LLM prompt 摘要
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


# ── 标签解析 ─────────────────────────────────────────────────────────────────

def parse_ontology_tags(value: str) -> list[str]:
    """
    解析 CLNDISDB / CLNDN 列的复杂标签。

    格式示例:
      CLNDISDB: "MedGen:CN169374|MONDO:MONDO:0014502,MedGen:C4015293,OMIM:616126,Orphanet:319563"
      CLNDN:    "not_specified|Mendelian_susceptibility_to...|not_provided"

    规则:
      - | 分隔不同的疾病/表型组
      - , 在每个组内分隔多个 ontology 标签
      - 标签格式: MedGen:xxx / MONDO:xxx / OMIM:xxx / Orphanet:xxx / MeSH:xxx
      - 也可能直接是疾病名称（无前缀）

    返回: 所有唯一标签的列表（去重）
    """
    if pd.isna(value) or str(value).strip() == "":
        return []

    tags = []
    value = str(value).strip()

    # 先按 | 分组
    groups = value.split("|")

    for group in groups:
        # 组内按 , 分隔
        parts = group.split(",")
        for part in parts:
            part = part.strip()
            if not part:
                continue
            # 提取 ontology 标签（MedGen: / MONDO: / OMIM: / Orphanet: / MeSH:）
            if ":" in part:
                # 保留完整标签，如 "MedGen:CN169374"
                tags.append(part)
            else:
                # 纯疾病名称，如 "Mendelian_susceptibility_to_mycobacterial_diseases..."
                # 转下划线为空格，转小写，作为简化标签
                tags.append(part.replace("_", " ").lower())

    # 去重，保持顺序
    seen = set()
    unique = []
    for t in tags:
        if t not in seen:
            seen.add(t)
            unique.append(t)

    return unique


def parse_clndisdb_tags(value: str) -> list[str]:
    """专门解析 CLNDISDB，只保留 ontology 前缀标签。"""
    if pd.isna(value) or str(value).strip() == "":
        return []

    tags = []
    groups = str(value).strip().split("|")
    for group in groups:
        for part in group.split(","):
            part = part.strip()
            if not part or ":" not in part:
                continue
            tags.append(part)  # MedGen:xxx, MONDO:xxx, OMIM:xxx, Orphanet:xxx, MeSH:xxx

    seen = set()
    unique = []
    for t in tags:
        if t not in seen:
            seen.add(t)
            unique.append(t)
    return unique


def parse_clndn_tags(value: str) -> list[str]:
    """专门解析 CLNDN，包含 ontology 标签 + 疾病名称。"""
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
                # 疾病名称转小写
                tags.append(part.replace("_", " ").lower())

    seen = set()
    unique = []
    for t in tags:
        if t not in seen:
            seen.add(t)
            unique.append(t)
    return unique


# ── 统计 ──────────────────────────────────────────────────────────────────────

def summarize_stage1(csv_path: str, sample_limit: int = 5_000_000) -> dict:
    """
    读取 Stage 1 CSV，生成 LLM 摘要。

    Parameters
    ----------
    csv_path : str
        3_clinvar_annotated.csv 路径
    sample_limit : int
        超过此行数时随机采样（避免太大）

    Returns
    -------
    dict: {
        "n_total": int,
        "n_with_clinvar": int,
        "feature_type_counts": dict,
        "clnsig_counts": dict,
        "clnrevstat_counts": dict,
        "clndisdb_unique_tags": list[str],   # 所有 unique ontology 标签
        "clndn_unique_tags": list[str],       # 所有 unique 标签（含疾病名）
        "clndisdb_tag_counts": dict,         # 标签出现频次
        "clndn_tag_counts": dict,
        "sample_rows": list[dict],            # 随机几行样例
    }
    """
    logger.info(f"Loading Stage 1 CSV: {csv_path}")
    df = pd.read_csv(csv_path, dtype=str, low_memory=False)
    n_total = len(df)
    logger.info(f"  Total rows: {n_total:,}")

    # ── 基本计数 ──────────────────────────────────────────────────────────
    has_clinvar = df["CLNSIG"].notna() & (df["CLNSIG"] != "")
    n_with_clinvar = has_clinvar.sum()

    feature_counts = df["feature_type"].value_counts().to_dict()
    clnsig_counts = df["CLNSIG"].value_counts().to_dict()
    clnrevstat_counts = df["CLNREVSTAT"].value_counts().to_dict()

    # ── CLNDISDB 解析 ─────────────────────────────────────────────────────
    logger.info("  Parsing CLNDISDB tags...")
    all_disdb_tags = []
    disdb_tag_counter = Counter()
    all_disdb_tags_raw = []

    for val in df["CLNDISDB"].dropna():
        tags = parse_clndisdb_tags(val)
        all_disdb_tags_raw.extend(tags)
        disdb_tag_counter.update(tags)

    # 取 top 100 常用标签（避免 prompt 爆炸）
    top_disdb = dict(disdb_tag_counter.most_common(100))
    logger.info(f"  CLNDISDB unique tags: {len(disdb_tag_counter)}, top: {len(top_disdb)}")

    # ── CLNDN 解析 ────────────────────────────────────────────────────────
    logger.info("  Parsing CLNDN tags...")
    clndn_tag_counter = Counter()
    for val in df["CLNDN"].dropna():
        tags = parse_clndn_tags(val)
        clndn_tag_counter.update(tags)

    top_clndn = dict(clndn_tag_counter.most_common(100))
    logger.info(f"  CLNDN unique tags: {len(clndn_tag_counter)}, top: {len(top_clndn)}")

    # ── 样例行 ─────────────────────────────────────────────────────────────
    if n_total > sample_limit:
        sample_df = df.sample(n=sample_limit, random_state=42)
    else:
        sample_df = df

    sample_rows = sample_df.head(5)[
        ["CPRA", "gene_name", "feature_type", "CLNSIG", "CLNREVSTAT",
         "CLNDISDB", "CLNDN"]
    ].fillna("").to_dict("records")

    result = {
        "n_total": int(n_total),
        "n_with_clinvar": int(n_with_clinvar),
        "feature_type_counts": feature_counts,
        "clnsig_counts": clnsig_counts,
        "clnrevstat_counts": clnrevstat_counts,
        "clndisdb_unique_tags": list(top_disdb.keys()),
        "clndn_unique_tags": list(top_clndn.keys()),
        "clndisdb_tag_counts": top_disdb,
        "clndn_tag_counts": top_clndn,
        "sample_rows": sample_rows,
    }

    return result


def build_llm_prompt(symptoms: str, summary: dict) -> str:
    """
    生成发给 LLM 的打分 prompt。

    summary 包含:
      - n_total, n_with_clinvar
      - feature_type_counts, clnsig_counts, clnrevstat_counts
      - clndisdb_unique_tags, clndn_unique_tags, clndisdb_tag_counts, clndn_tag_counts
      - sample_rows
    """
    lines = [
        f"患者症状: {symptoms}",
        f"",
        f"变异数据统计 (共 {summary['n_total']:,} 行, 其中 {summary['n_with_clinvar']:,} 行有 ClinVar 注释):",
        f"",
        f"【feature_type 分布】",
    ]
    for k, v in summary["feature_type_counts"].items():
        lines.append(f"  {k}: {v:,}")

    lines.extend([
        f"",
        f"【CLNSIG 分布】(ClinVar 致病性评级)",
    ])
    for k, v in summary["clnsig_counts"].items():
        lines.append(f"  {k}: {v:,}")

    lines.extend([
        f"",
        f"【CLNREVSTAT 分布】(ClinVar 证据状态)",
    ])
    for k, v in summary["clnrevstat_counts"].items():
        lines.append(f"  {k}: {v:,}")

    lines.extend([
        f"",
        f"【CLNDISDB 标签出现频次】(top 80, 格式: ontology:code)",
        f"（每个标签代表一种疾病/表型关联，按频次排列）",
    ])
    for tag, cnt in list(summary["clndisdb_tag_counts"].items())[:80]:
        lines.append(f"  [{cnt:,}] {tag}")

    lines.extend([
        f"",
        f"【CLNDN 标签出现频次】(top 80, 格式: ontology:code 或疾病名)",
    ])
    for tag, cnt in list(summary["clndn_tag_counts"].items())[:80]:
        lines.append(f"  [{cnt:,}] {tag}")

    lines.extend([
        f"",
        f"【样例行】",
    ])
    for row in summary["sample_rows"]:
        lines.append(f"  {row}")

    lines.extend([
        f"",
        f"请根据患者症状，为以上各列的每个取值给出 0~1 的相关性分数（1=最相关，0=不相关）：",
        f"  1. feature_type 各取值的分数",
        f"  2. CLNSIG 各取值的分数",
        f"  3. CLNREVSTAT 各取值的分数",
        f"  4. CLNDISDB 每个标签在本病例下的相关性分数",
        f"  5. CLNDN 每个标签在本病例下的相关性分数",
        f"  6. 各列的权重比例 (col_weights: feature_type/CLNSIG/CLNREVSTAT/CLNDISDB/CLNDN，合计=1.0)",
        f"",
        f"返回 JSON 格式：",
        f'''{{
  "feature_type_scores": {{"exon": 1.0, "CDS": 1.0, ...}},
  "clnsig_scores": {{"Pathogenic": 1.0, ...}},
  "clnrevstat_scores": {{"criteria_provided,_multiple_submitters": 0.9, ...}},
  "clndisdb_tag_scores": {{"MONDO:0014502": 0.95, "OMIM:616126": 0.8, ...}},
  "clndn_tag_scores": {{"Mendelian disease related": 0.7, ...}},
  "col_weights": {{"feature_type": 0.15, "CLNSIG": 0.30, "CLNREVSTAT": 0.10, "CLNDISDB": 0.30, "CLNDN": 0.15}}
}}''',
    ])

    return "\n".join(lines)
