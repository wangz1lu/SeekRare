"""
eqtl_annotation.py — Stage 2: GTEx eQTL 注释

用户从 GTEx 下载 eQTL parquet 数据（每组织一个 .parquet 文件），
LLM 根据患者症状筛选相关组织，然后对 Stage 1 CSV 进行 eQTL 注释。

GTEx 下载地址: https://gtexportal.org/downloads
文件: GTEx_Analysis_v11_eQTL.tar → *.eQTLs.signif_pairs.parquet

Stage 2 设计:
1. LLM 接收患者症状 + 所有可用组织列表 → 选择相关组织
2. 对选中的组织，运行 eQTL 匹配
3. eQTL 注释以新增列形式追加到 Stage 1 CSV
   - 同一基因有多条 eQTL → 取 p-value 最显著的那条
   - 没有 eQTL 注释的位点 → 列值留空
"""

from __future__ import annotations

import gzip
import json
import os
import glob
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


# ── GTEx 组织列表（56个组织）──────────────────────────────────────────────
GTEX_TISSUES = [
    "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland",
    "Artery_Aorta", "Artery_Coronary", "Artery_Tibial",
    "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
    "Breast_Mammary_Tissue", "Skin_Not_Sun_Exposed_Suprapubic",
    "Skin_Sun_Exposed_Lower_leg", "Cells_Cultured_fibroblasts",
    "Cells_EBV-transformed_lymphocytes",
    "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction",
    "Esophagus_Mucosa", "Esophagus_Muscularis",
    "Heart_Atrial_Appendage", "Heart_Left_Ventricle",
    "Kidney_Cortex", "Kidney_Medulla",
    "Liver", "Lung", "Minor_Salivary_Gland",
    "Muscle_Skeletal", "Nerve_Tibial",
    "Ovary", "Pancreas", "Pituitary", "Prostate",
    "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis",
    "Thyroid", "Uterus", "Vagina",
    "Whole_Blood",
]

# 组织→HPO 关联（帮助 LLM 判断哪些组织相关）
TISSUE_HPO_HINT = {
    "Eye": ["Retina", "Choroid", "Optic nerve"],
    "Brain": ["Brain", "Neuron", "Cortex", "Cerebellum"],
    "Heart": ["Heart", "Cardiac", "Myocardium"],
    "Liver": ["Liver", "Hepatocyte"],
    "Kidney": ["Kidney", "Nephron", "Glomerulus"],
    "Muscle": ["Muscle", "Skeletal muscle"],
    "Nerve": ["Nerve", "Peripheral neuropathy"],
    "Blood": ["Blood", "Lymphocyte", "Neutrophil"],
    "Skin": ["Skin", "Epidermis", "Dermis"],
    "Lung": ["Lung", "Alveolus", "Bronchus"],
}


def get_available_tissues(tissue_dir: str) -> list[str]:
    """扫描 tissue_dir，返回所有可用的 GTEx parquet 文件对应的组织名。"""
    tissue_files = glob.glob(os.path.join(tissue_dir, "*.eQTLs.signif_pairs.parquet"))
    tissues = []
    for f in tissue_files:
        # 文件名格式: Brain_Amygdala.v11.signif_pairs.parquet
        name = os.path.basename(f).split(".v11")[0]
        tissues.append(name)
    logger.info(f"发现 {len(tissues)} 个 GTEx 组织: {tissues}")
    return tissues


def build_llm_tissue_prompt(symptoms: str, tissues: list[str]) -> str:
    """生成 LLM 组织选择 prompt。"""
    tissue_str = "\n".join(f"  - {t}" for t in tissues)
    return f"""患者症状: {symptoms}

GTEx 数据库有以下 {len(tissues)} 个可用组织:
{tissue_str}

请根据患者症状，选择最相关的 3~8 个组织。
只返回与症状可能相关的组织（例如神经系统症状选 Brain 相关组织，眼科症状选 Eye 相关组织）。

返回 JSON 数组格式（只列出组织名）:
{{"selected_tissues": ["Brain_Cortex", "Brain_Hippocampus", "Thyroid"]}}
"""


def parse_variant_id(chrom: str, pos: str, ref: str, alt: str, suffix: str) -> str:
    """构建 GTEx variant_id 格式: chr1_12345_A_G_b38"""
    chrom = str(chrom).lstrip("chr").replace("chr", "")
    return f"chr{chrom}_{pos}_{ref.upper()}_{alt.upper()}{suffix}"


def run_eqtl_annotation(
    stage1_csv: str,
    tissue_dir: str,
    selected_tissues: list[str],
    output_csv: Optional[str] = None,
    min_pval: float = 1e-3,
) -> pd.DataFrame:
    """
    对 Stage 1 CSV 进行 eQTL 注释。

    参数
    ----
    stage1_csv : str
        Stage 1 输出 CSV (3_clinvar_annotated.csv)
    tissue_dir : str
        GTEx eQTL parquet 文件目录（每组织一个 .parquet）
    selected_tissues : list[str]
        LLM 筛选出的相关组织列表
    output_csv : str, optional
        输出路径
    min_pval : float
        eQTL p-value 阈值 (default: 1e-3)

    返回
    ----
    pd.DataFrame
        Stage 1 CSV 基础上追加 eQTL 列:
        eqtl_gene, eqtl_slope, eqtl_pval, eqtl_tissue
        （无 eQTL 注释则留空）
    """
    logger.info(f"Stage 2: GTEx eQTL annotation ({len(selected_tissues)} tissues)")
    logger.info(f"  Tissue dir: {tissue_dir}")

    # ── 1. 检测 suffix ──────────────────────────────────────────────────────
    tissue_files = glob.glob(os.path.join(tissue_dir, "*.eQTLs.signif_pairs.parquet"))
    if not tissue_files:
        raise FileNotFoundError(f"未找到 GTEx parquet 文件: {tissue_dir}")

    # 找第一个在 selected_tissues 里的组织
    sample_file = None
    for tf in tissue_files:
        t_name = os.path.basename(tf).split(".v11")[0]
        if t_name in selected_tissues:
            sample_file = tf
            break

    if sample_file is None:
        # Fallback: 用第一个文件检测 suffix
        sample_file = tissue_files[0]
        logger.warning(f"未找到匹配组织的 parquet，用第一个代替: {sample_file}")

    sample_row = pd.read_parquet(sample_file, columns=["variant_id"]).iloc[0]
    variant_id_sample = sample_row["variant_id"]
    # e.g. "chr1_12345_A_G_b38" → "_b38"
    suffix = "_" + variant_id_sample.rsplit("_", 1)[-1]
    logger.info(f"  GTEx variant_id suffix: '{suffix}'")

    # ── 2. 加载 Stage 1 CSV ────────────────────────────────────────────────
    df = pd.read_csv(stage1_csv, dtype=str)
    logger.info(f"  Stage 1 rows: {len(df)}")

    # 构建 variant_id 列
    chrom_col = "CHROM" if "CHROM" in df.columns else "chrom"
    df["_vid"] = (
        "chr" + df[chrom_col].astype(str).str.lstrip("chr")
        + "_" + df["POS"].astype(str)
        + "_" + df["REF"].str.upper()
        + "_" + df["ALT"].str.upper()
        + suffix
    )
    vcf_set = set(df["_vid"])

    # ── 3. 初始化 eQTL 列 ──────────────────────────────────────────────────
    for col in ["eqtl_gene", "eqtl_slope", "eqtl_pval", "eqtl_tissue", "n_eqtl_tissues"]:
        if col not in df.columns:
            df[col] = None

    # ── 4. 逐组织处理 ──────────────────────────────────────────────────────
    n_matched_total = 0

    for tissue_file in tissue_files:
        t_name = os.path.basename(tissue_file).split(".v11")[0]
        if t_name not in selected_tissues:
            continue

        logger.info(f"  --> {t_name}")
        try:
            chunk = pd.read_parquet(
                tissue_file,
                columns=["variant_id", "phenotype_id", "slope", "pval_nominal"],
            )
            chunk = chunk[chunk["pval_nominal"] < min_pval]

            matched = chunk[chunk["variant_id"].isin(vcf_set)]
            if matched.empty:
                logger.info(f"      无匹配")
                continue

            # 取每个 variant 最显著的 eQTL
            best = (
                matched.sort_values("pval_nominal")
                .groupby("variant_id")
                .first()
                .reset_index()
            )

            # 合并到 df
            merge_cols = ["variant_id", "phenotype_id", "slope", "pval_nominal"]
            merged = df[["_vid"]].merge(
                best[merge_cols].rename(columns={
                    "variant_id": "_vid",
                    "phenotype_id": f"eqtl_gene_{t_name}",
                    "slope": f"eqtl_slope_{t_name}",
                    "pval_nominal": f"eqtl_pval_{t_name}",
                }),
                on="_vid",
                how="left",
            )

            # 更新有 eQTL 的行（取最显著的组织）
            for idx in merged[merged[f"eqtl_gene_{t_name}"].notna()].index:
                row_equiv = df.loc[idx]
                # 如果该行已有 eQTL，取 p-value 更小的
                if pd.isna(df.at[idx, "eqtl_gene"]) or \
                   merged.at[idx, f"eqtl_pval_{t_name}"] < float(df.at[idx, "eqtl_pval"] or 1):
                    df.at[idx, "eqtl_gene"] = merged.at[idx, f"eqtl_gene_{t_name}"]
                    df.at[idx, "eqtl_slope"] = merged.at[idx, f"eqtl_slope_{t_name}"]
                    df.at[idx, "eqtl_pval"] = merged.at[idx, f"eqtl_pval_{t_name}"]
                    df.at[idx, "eqtl_tissue"] = t_name
                    n_matched_total += 1

            logger.info(f"      {len(best)} variants with eQTL")
            del chunk, matched, best, merged

        except Exception as e:
            logger.warning(f"      失败: {e}")
            continue

    # ── 5. 统计匹配数 ──────────────────────────────────────────────────────
    n_eqtl = df["eqtl_gene"].notna().sum()
    logger.info(f"  eQTL annotated: {n_eqtl}/{len(df)} variants")
    df.drop(columns=["_vid"], inplace=True, errors="ignore")

    # ── 6. 保存 ────────────────────────────────────────────────────────────
    if output_csv:
        df.to_csv(output_csv, index=False)
        logger.info(f"  Stage 2 输出: {output_csv}")

    return df
