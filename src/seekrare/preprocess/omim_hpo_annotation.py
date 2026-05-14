"""
Stage 2 后处理：OMIM + HPO 二次注释。

流程：
  1. OMIM 注释：根据 gene_name 从 genemap2.txt 查 MIM 号，更新 OMIM 列
     （若 CSV 已有 OMIM 值，保留不覆盖）
  2. HPO 注释：用 phenotype.hpoa（database_id=OMIM:xxxx 映射 hpo_id=disease_name）
     更新 HPO 列（追加，不覆盖已有 HPO）
  3. disease_name：用 phenotype.hpoa 的 disease_name 补充 diseasename 列
     （若 CSV 已有值，保留不覆盖）

Usage:
    from seekrare.preprocess.omim_hpo_annotation import stage2_omim_hpo_annotation
    stage2_omim_hpo_annotation("stage2_output.csv", ...)
"""

import re
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


# ─────────────────────────────────────────────────────────────────────────────
# 内部解析函数
# ─────────────────────────────────────────────────────────────────────────────

def build_gene_to_mim(genemap2_path: str) -> dict[str, str]:
    """
    genemap2.txt → gene_name → MIM号列表（逗号分隔）
    列：第9列 approved_symbol，第6列 entrez，隔12列 phenotype
    """
    gene_to_mim = {}
    with open(genemap2_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 13:
                continue
            approved = parts[8].strip()
            gene = approved if (approved and approved != "-") else parts[6].split(",")[0].strip()
            if not gene:
                continue
            pheno = parts[12] if len(parts) > 12 else ""
            mims = re.findall(r"(?:,\s*|\s)(\d{6})\s*\(", pheno)
            if mims:
                seen = set()
                unique = []
                for m in mims:
                    if m not in seen:
                        seen.add(m)
                        unique.append(m)
                gene_to_mim[gene] = ",".join(unique)
    return gene_to_mim


def build_mim_to_name(mimtitles_path: str) -> dict[str, str]:
    """mimTitles.txt → MIM号 → 疾病名称"""
    mim_to_name = {}
    with open(mimtitles_path, "r") as f:
        for line in f:
            if line.startswith("#") or line.startswith("Prefix"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            mim_num = parts[1].strip()
            title = parts[2].strip()
            if ";" in title:
                title = title.split(";")[0].strip()
            if mim_num and title:
                mim_to_name[mim_num] = title
    return mim_to_name


def build_hpo_map(phenotype_hpoa_path: str) -> dict[str, list[tuple[str, str]]]:
    """
    phenotype.hpoa → {OMIM:xxxx: [(hpo_id, disease_name), ...]}
    关键列：database_id（如 OMIM:607745），hpo_id（如 HP:0001250），disease_name
    """
    omim_to_hpos = {}
    with open(phenotype_hpoa_path, "r") as f:
        header = None
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if header is None:
                header = cols
                try:
                    db_idx = header.index("database_id")
                    hpo_idx = header.index("hpo_id")
                    dn_idx = header.index("disease_name")
                except ValueError:
                    logger.error(f"[HPO] phenotype.hpoa 缺少必要列: {header}")
                    return {}
                continue
            if len(cols) <= max(db_idx, hpo_idx, dn_idx):
                continue
            db_id = cols[db_idx].strip()
            hpo_id = cols[hpo_idx].strip()
            disease_name = cols[dn_idx].strip()

            if not db_id or not hpo_id or db_id == "NOT":
                continue
            if not db_id.startswith("OMIM:"):
                continue

            mim = db_id.replace("OMIM:", "").strip()
            if mim not in omim_to_hpos:
                omim_to_hpos[mim] = []
            omim_to_hpos[mim].append((hpo_id, disease_name))

    return omim_to_hpos


# ─────────────────────────────────────────────────────────────────────────────
# 主函数
# ─────────────────────────────────────────────────────────────────────────────

def stage2_omim_hpo_annotation(
    input_csv: Union[str, Path],
    genemap2_path: Optional[str] = None,
    mimtitles_path: Optional[str] = None,
    phenotype_hpoa_path: Optional[str] = None,
    output_csv: Optional[Union[str, Path]] = None,
    gene_col: str = "gene_name",
    omim_col: str = "OMIM",
    hpo_col: str = "HPO",
    disease_col: str = "diseasename",
) -> pd.DataFrame:
    """
    Stage 2 后处理：OMIM + HPO 二次注释。

    Parameters
    ----------
    input_csv : str   Stage 2 输出 CSV（含 gene_name, OMIM, HPO, diseasename 列）
    genemap2_path : str   genemap2.txt（ OMIM: gene → MIM 号）
    mimtitles_path : str   mimTitles.txt（MIM 号 → 疾病名）
    phenotype_hpoa_path : str   phenotype.hpoa（OMIM → HPO 列表 + disease_name）
    output_csv : str, optional   输出路径（默认覆盖 input_csv）
    gene_col, omim_col, hpo_col, disease_col : 列名（默认标准列名）

    规则：
      - OMIM 列：已有值则保留，不覆盖；无值则用 genemap2 补充
      - disease_name 列：已有值则保留，不覆盖；无值则用 mimTitles 或 phenotype.hpoa 补充
      - HPO 列：追加新 HPO（用 ; 分隔），不覆盖已有
    """
    input_csv = str(input_csv)
    if output_csv is None:
        output_csv = input_csv
    else:
        output_csv = str(output_csv)

    logger.info(f"[OMIM/HPO] Loading: {input_csv}")
    df = pd.read_csv(input_csv, dtype=str)

    # 确保目标列存在
    for col in [omim_col, hpo_col, disease_col]:
        if col not in df.columns:
            df[col] = ""

    n = len(df)

    # ── Step 1: OMIM 注释（genemap2）────────────────────────────────────────
    if genemap2_path and Path(genemap2_path).exists():
        logger.info(f"[OMIM] Building gene→MIM from: {genemap2_path}")
        gene_to_mim = build_gene_to_mim(genemap2_path)
        omim_to_name = build_mim_to_name(mimtitles_path) if mimtitles_path and Path(mimtitles_path).exists() else {}

        updated_omim = 0
        for idx, row in df.iterrows():
            current_omim = str(row.get(omim_col, "")).strip()
            # 已有 OMIM → 跳过（不覆盖）
            if current_omim and current_omim not in ("nan", ""):
                continue

            gene = str(row.get(gene_col, "")).strip()
            if not gene:
                continue

            new_mim = gene_to_mim.get(gene, "")
            if not new_mim:
                continue

            df.at[idx, omim_col] = new_mim

            # 补充 disease_name（若缺失）
            current_dn = str(row.get(disease_col, "")).strip()
            if current_dn in ("nan", ""):
                names = [omim_to_name.get(m, f"OMIM:{m}") for m in new_mim.split(",")]
                df.at[idx, disease_col] = "|".join(names)

            updated_omim += 1

        logger.info(f"[OMIM] Updated {updated_omim}/{n} rows")
    else:
        logger.info("[OMIM] Skipped (genemap2_path not provided or not found)")

    # ── Step 2: HPO + disease_name 注释（phenotype.hpoa）────────────────────
    if phenotype_hpoa_path and Path(phenotype_hpoa_path).exists():
        logger.info(f"[HPO] Building OMIM→HPO map from: {phenotype_hpoa_path}")
        omim_to_hpos = build_hpo_map(phenotype_hpoa_path)
        logger.info(f"  {len(omim_to_hpos)} OMIM entries with HPO annotations")

        updated_hpo = 0
        updated_dn = 0
        for idx, row in df.iterrows():
            omim_val = str(row.get(omim_col, "")).strip()
            if not omim_val or omim_val in ("nan", ""):
                continue

            # 解析所有 OMIM 号（逗号分隔）
            mim_list = [m.strip() for m in omim_val.split(",") if m.strip()]

            new_hpos = []
            new_dn_candidates = []

            for mim in mim_list:
                entries = omim_to_hpos.get(mim, [])
                for hpo_id, disease_name in entries:
                    if hpo_id not in new_hpos:
                        new_hpos.append(hpo_id)
                    if disease_name:
                        new_dn_candidates.append(disease_name)

            # 更新 HPO 列（追加，不覆盖）
            if new_hpos:
                current_hpo = str(row.get(hpo_col, "")).strip()
                existing = current_hpo.split(";") if current_hpo and current_hpo not in ("nan", "") else []
                existing = [h.strip() for h in existing if h.strip()]
                for hp in new_hpos:
                    if hp not in existing:
                        existing.append(hp)
                df.at[idx, hpo_col] = ";".join(existing)
                updated_hpo += 1

            # 更新 disease_name（取第一个新 disease_name；已有则保留）
            if new_dn_candidates:
                current_dn = str(row.get(disease_col, "")).strip()
                if current_dn in ("nan", ""):
                    df.at[idx, disease_col] = new_dn_candidates[0]
                    updated_dn += 1

        logger.info(f"[HPO] Updated HPO for {updated_hpo}/{n} rows")
        logger.info(f"[HPO] Updated disease_name for {updated_dn}/{n} rows")
    else:
        logger.info("[HPO] Skipped (phenotype_hpoa_path not provided or not found)")

    df.to_csv(output_csv, index=False)
    logger.info(f"[OMIM/HPO] Saved: {output_csv}")

    return df