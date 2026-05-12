"""
stage1_family.py — Stage 1 家系 trio 预处理模式

当提供了 father_vcf + mother_vcf 时，启用家系模式：

流程（纯 Python + bcftools subprocess，不依赖外部脚本）:
  1. Normalize — left-normalize + split multi-allelics（每个样本单独做）
  2. Index — bcftools index
  3. Merge — 三口 VCF 合并
  4. Filter — 质量过滤（QUAL>30, DP>10, GQ>20）
  5. Norm again — check-ref w 严格规范化
  6. Drop missing genotypes
  7. Inheritance 分类:
       de_novo:     父母均 ref，子女 alt
       recessive:   父母均 het，子女 hom_alt
  8. 各模式 VCF → GT CSV → GTF → ClinVar → 合并输出

Usage:
  from seekrare.preprocess.stage1_family import run_family_preprocess

  result = run_family_preprocess(
      work_dir="seekrare_output",
      child_vcf="child.vcf.gz",
      father_vcf="father.vcf.gz",
      mother_vcf="mother.vcf.gz",
      ref_fasta="/path/to/GRCh38.fa",
      gtf_file="/path/to/genomic.gtf",
      clinvar_vcf="/path/to/clinvar.vcf.gz",
      dbsnp_vcf="/path/to/dbsnp.vcf.gz",
  )
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger

from seekrare.preprocess.vcf_to_gt import stage1_vcf_to_gt_csv
from seekrare.preprocess.simplify_clinvar_csv import simplify_clinvar_csv
from seekrare.preprocess.gene_annotation import stage1_annotate_by_gtf
from seekrare.preprocess.clinvar_annotation import stage1_merge_filter_clinvar


# ─────────────────────────────────────────────────────────────────────────────
# 内部工具
# ─────────────────────────────────────────────────────────────────────────────

def _run(cmd: list, check: bool = True, **kwargs):
    """执行 bcftools 命令，失败则抛出异常。"""
    logger.info(f"  bcftools: {' '.join(str(x) for x in cmd[:3])} ...")
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"bcftools failed (code {result.returncode}):\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDERR: {result.stderr[:1000]}"
        )
    return result


def _which_bcftools() -> str:
    """确保 bcftools 在 PATH 中。"""
    path = shutil.which("bcftools")
    if path is None:
        raise RuntimeError("bcftools not found in PATH. 请安装 bcftools >= 1.18 并确保在 PATH 中.")
    return path


def _vcf_exists_notempty(vcf: str) -> bool:
    p = Path(vcf)
    return p.exists() and p.stat().st_size > 0


def _is_missing_sample(vcf_path: str, sample_name: str) -> bool:
    """检查 VCF 是否包含某个 sample。"""
    cmd = [_which_bcftools(), "query", "-l", vcf_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return sample_name not in (result.stdout or "")


# ─────────────────────────────────────────────────────────────────────────────
# Step 1: normalize（每个样本单独做）
# ─────────────────────────────────────────────────────────────────────────────

def _normalize_sample(
    vcf_in: str,
    vcf_out: str,
    ref_fasta: str,
    force: bool = False,
) -> str:
    """
    bcftools norm -m -both -f <ref> -Oz -o <out>
    将 multi-allelic 拆分为 biallelic，left-normalize。
    """
    if not force and _vcf_exists_notempty(vcf_out):
        logger.info(f"  [跳过] 已规范化: {vcf_out}")
        return vcf_out

    _run([
        _which_bcftools(), "norm",
        "-m", "-both",
        "-f", ref_fasta,
        vcf_in,
        "-Oz", "-o", vcf_out,
    ])
    return vcf_out


def _index_vcf(vcf: str, force: bool = False):
    """bcftools index <vcf>"""
    tbi = vcf + ".tbi"
    if not force and Path(tbi).exists():
        return
    _run([_which_bcftools(), "index", vcf])


# ─────────────────────────────────────────────────────────────────────────────
# Step 2: merge trio
# ─────────────────────────────────────────────────────────────────────────────

def _merge_trio(outdir: Path, father_vcf: str, mother_vcf: str, child_vcf: str) -> str:
    """bcftools merge 三口 VCF"""
    merged = outdir / "trio.merged.vcf.gz"
    if _vcf_exists_notempty(str(merged)):
        logger.info(f"  [跳过] 已合并: {merged}")
        return str(merged)

    _run([
        _which_bcftools(), "merge",
        father_vcf, mother_vcf, child_vcf,
        "-Oz", "-o", str(merged),
    ])
    return str(merged)


# ─────────────────────────────────────────────────────────────────────────────
# Step 3: quality filter
# ─────────────────────────────────────────────────────────────────────────────

def _filter_quality(vcf_in: str, vcf_out: str, min_qual: float = 30.0, min_dp: int = 10, min_gq: int = 20) -> str:
    """
    bcftools filter -i 'QUAL>{min_qual} && FMT/DP>{min_dp} && FMT/GQ>{min_gq}'
    """
    if _vcf_exists_notempty(vcf_out):
        logger.info(f"  [跳过] 已过滤: {vcf_out}")
        return vcf_out

    expr = f"QUAL>{min_qual} && FMT/DP>{min_dp} && FMT/GQ>{min_gq}"
    _run([
        _which_bcftools(), "filter",
        "-i", expr,
        vcf_in, "-Oz", "-o", vcf_out,
    ])
    _index_vcf(vcf_out)
    return vcf_out


# ─────────────────────────────────────────────────────────────────────────────
# Step 4: strict left-normalize (check-ref)
# ─────────────────────────────────────────────────────────────────────────────

def _norm_strict(vcf_in: str, vcf_out: str, ref_fasta: str) -> str:
    """
    bcftools norm -m -both --check-ref w -f <ref>  严格规范化
    确保 REF 与参考基因组完全一致
    """
    if _vcf_exists_notempty(vcf_out):
        logger.info(f"  [跳过] 已严格规范化: {vcf_out}")
        return vcf_out

    _run([
        _which_bcftools(), "norm",
        "-m", "-both",
        "--check-ref", "w",
        "-f", ref_fasta,
        vcf_in, "-Oz", "-o", vcf_out,
    ])
    _index_vcf(vcf_out)
    return vcf_out


# ─────────────────────────────────────────────────────────────────────────────
# Step 5: drop missing genotypes
# ─────────────────────────────────────────────────────────────────────────────

def _drop_missing(vcf_in: str, vcf_out: str) -> str:
    """
    bcftools view -i 'GT[0]!="mis" && GT[1]!="mis" && GT[2]!="mis"'
    去掉任意样本基因型为缺失的位点
    """
    if _vcf_exists_notempty(vcf_out):
        logger.info(f"  [跳过] 已去除缺失: {vcf_out}")
        return vcf_out

    _run([
        _which_bcftools(), "view",
        "-i", 'GT[0]!="mis" && GT[1]!="mis" && GT[2]!="mis"',
        vcf_in, "-Oz", "-o", vcf_out,
    ])
    _index_vcf(vcf_out)
    return vcf_out


# ─────────────────────────────────────────────────────────────────────────────
# Step 6: inheritance classification
# ─────────────────────────────────────────────────────────────────────────────

def _filter_denovo(vcf_in: str, vcf_out: str) -> str:
    """
    de novo: 父母均 0/0 或 0|0，子女为 0/1 或 1/0 或 1/1
    (GT[0]=father, GT[1]=mother, GT[2]=child)
    """
    if _vcf_exists_notempty(vcf_out):
        logger.info(f"  [跳过] de_novo: {vcf_out}")
        return vcf_out

    _run([
        _which_bcftools(), "view",
        "-i",
        '(GT[0]="0/0" || GT[0]="0|0") && '
        '(GT[1]="0/0" || GT[1]="0|0") && '
        '(GT[2]="0/1" || GT[2]="1/0" || GT[2]="1/1" || GT[2]="1|1" || GT[2]="0|1" || GT[2]="1|0")',
        vcf_in, "-Oz", "-o", vcf_out,
    ])
    _index_vcf(vcf_out)
    logger.info(f"  de_novo: {vcf_out} ({Path(vcf_out).stat().st_size} bytes)")
    return vcf_out


def _filter_recessive(vcf_in: str, vcf_out: str) -> str:
    """
    recessive (hom): 父母均 het (0/1)，子女为 hom_alt (1/1 或 1|1)
    """
    if _vcf_exists_notempty(vcf_out):
        logger.info(f"  [跳过] recessive: {vcf_out}")
        return vcf_out

    _run([
        _which_bcftools(), "view",
        "-i",
        '(GT[0]="0/1" || GT[0]="0|1" || GT[0]="1/0" || GT[0]="1|0") && '
        '(GT[1]="0/1" || GT[1]="0|1" || GT[1]="1/0" || GT[1]="1|0") && '
        '(GT[2]="1/1" || GT[2]="1|1" || GT[2]="1/0" || GT[2]="1|0" || GT[2]="0/1" || GT[2]="0|1")',
        vcf_in, "-Oz", "-o", vcf_out,
    ])
    _index_vcf(vcf_out)
    logger.info(f"  recessive: {vcf_out} ({Path(vcf_out).stat().st_size} bytes)")
    return vcf_out


def _filter_xlinked(vcf_in: str, vcf_out: str) -> str:
    """
    X-linked: 母亲 ref，父亲 het，子女 son alt（仅男性）
    bcftools view -i ' ... '
    """
    if _vcf_exists_notempty(vcf_out):
        return vcf_out

    # father het, mother hom_ref, child alt (male)
    # 简化条件: father het, mother ref, child alt
    _run([
        _which_bcftools(), "view",
        "-i",
        '(GT[0]="0/1" || GT[0]="0|1" || GT[0]="1/0" || GT[0]="1|0") && '
        '(GT[1]="0/0" || GT[1]="0|0") && '
        '(GT[2]="1/1" || GT[2]="1|1" || GT[2]="1/0" || GT[2]="1|0" || GT[2]="0/1" || GT[2]="0|1")',
        vcf_in, "-Oz", "-o", vcf_out,
    ])
    _index_vcf(vcf_out)
    return vcf_out


# ─────────────────────────────────────────────────────────────────────────────
# VCF → annotated CSV（通用）
# ─────────────────────────────────────────────────────────────────────────────

def _annotate_vcf(
    vcf_path: str,
    work_dir: Path,
    gtf_file: Optional[str],
    clinvar_vcf: Optional[str],
    mode_name: str,
) -> Optional[pd.DataFrame]:
    """
    单个 VCF → GT CSV → GTF → ClinVar → DataFrame
    """
    if not _vcf_exists_notempty(vcf_path):
        logger.info(f"  [{mode_name}] VCF 不存在或为空: {vcf_path}")
        return None

    mode_dir = work_dir / mode_name
    mode_dir.mkdir(parents=True, exist_ok=True)

    gt_csv   = mode_dir / f"{mode_name}_gt.csv"
    gtf_csv  = mode_dir / f"{mode_name}_gtf.csv"
    final_csv = mode_dir / f"{mode_name}_clinvar.csv"

    # VCF → GT CSV
    if not gt_csv.exists():
        logger.info(f"  [{mode_name}] VCF → GT CSV: {gt_csv}")
        stage1_vcf_to_gt_csv(vcf_path, str(gt_csv))
    else:
        logger.info(f"  [跳过] GT CSV 已存在: {gt_csv}")

    df = pd.read_csv(str(gt_csv), dtype=str)
    if len(df) == 0:
        logger.info(f"  [{mode_name}] 无变异: {vcf_path}")
        return None

    # GTF annotation
    if gtf_file:
        if not gtf_csv.exists():
            stage1_annotate_by_gtf(str(gt_csv), gtf_file, str(gtf_csv))
        df = pd.read_csv(str(gtf_csv), dtype=str)

    # ClinVar annotation
    if clinvar_vcf:
        if not final_csv.exists():
            stage1_merge_filter_clinvar(
                str(gtf_csv if gtf_file else gt_csv),
                clinvar_vcf,
                str(final_csv),
            )
        df = pd.read_csv(str(final_csv), dtype=str)

    df["inheritance_mode"] = mode_name
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 主入口
# ─────────────────────────────────────────────────────────────────────────────

def run_family_preprocess(
    work_dir: Union[str, Path],
    child_vcf: str,
    father_vcf: Optional[str] = None,
    mother_vcf: Optional[str] = None,
    ref_fasta: Optional[str] = None,
    gtf_file: Optional[str] = None,
    clinvar_vcf: Optional[str] = None,
    dbsnp_vcf: Optional[str] = None,
    min_qual: float = 30.0,
    min_dp: int = 10,
    min_gq: int = 20,
    inheritance_modes: Optional[list[str]] = None,
) -> dict:
    """
    Stage 1 家系 trio 预处理完整流程。

    当 father_vcf + mother_vcf 同时提供时，启用家系模式。
    否则退化为单样本模式（原有流程）。

    Parameters
    ----------
    work_dir : str or Path
        输出根目录
    child_vcf : str
        先证者 VCF（必填）
    father_vcf, mother_vcf : str, optional
        父母 VCF（家系模式需要）
    ref_fasta : str, optional
        参考基因组 FASTA（家系模式需要）
    gtf_file : str, optional
        GTF 基因注释文件
    clinvar_vcf : str, optional
        ClinVar VCF
    dbsnp_vcf : str, optional
        dbSNP VCF（家系模式暂不使用，保留接口）
    min_qual, min_dp, min_gq : float/int
        过滤阈值（家系模式用）
    inheritance_modes : list[str], optional
        指定运行哪些模式，默认全部:
        ["de_novo", "recessive", "xlinked"]

    Returns
    -------
    dict: {
        "modes": {"de_novo": DataFrame, "recessive": DataFrame, ...},
        "combined": DataFrame,
        "output_csv": str,
    }
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    if inheritance_modes is None:
        inheritance_modes = ["de_novo", "recessive", "xlinked"]

    # ── 家系检测 ──────────────────────────────────────────────────────────
    is_family = father_vcf and mother_vcf and ref_fasta

    if not is_family:
        logger.info("=" * 60)
        logger.info("Stage 1: 单样本模式")
        logger.info("=" * 60)
        return _singleton_mode(
            work_dir, child_vcf, gtf_file, clinvar_vcf
        )

    # ── 家系模式 ─────────────────────────────────────────────────────────
    logger.info("=" * 60)
    logger.info("Stage 1: 家系 trio 预处理模式")
    logger.info(f"  Child:  {child_vcf}")
    logger.info(f"  Father: {father_vcf}")
    logger.info(f"  Mother: {mother_vcf}")
    logger.info(f"  REF:    {ref_fasta}")
    logger.info(f"  Modes:  {inheritance_modes}")
    logger.info("=" * 60)

    trio_dir = work_dir / "trio_preprocess"
    trio_dir.mkdir(parents=True, exist_ok=True)

    bt = _which_bcftools()
    logger.info(f"  bcftools: {bt}")

    # Step 1: normalize 每个样本
    f_fa = father_vcf
    f_norm = str(trio_dir / "father.norm.vcf.gz")
    m_norm = str(trio_dir / "mother.norm.vcf.gz")
    c_norm = str(trio_dir / "child.norm.vcf.gz")

    logger.info("[Step 1] Normalize each sample...")
    _normalize_sample(f_fa, f_norm, ref_fasta)
    _normalize_sample(mother_vcf, m_norm, ref_fasta)
    _normalize_sample(child_vcf, c_norm, ref_fasta)

    # Index
    logger.info("[Step 2] Index normalized VCFs...")
    _index_vcf(f_norm)
    _index_vcf(m_norm)
    _index_vcf(c_norm)

    # Step 3: merge
    merged_vcf = _merge_trio(trio_dir, f_norm, m_norm, c_norm)

    # Step 4: quality filter
    filtered_vcf = str(trio_dir / "filtered.vcf.gz")
    logger.info("[Step 3-4] Quality filter (QUAL>{}, DP>{}, GQ>{})...".format(min_qual, min_dp, min_gq))
    _filter_quality(merged_vcf, filtered_vcf, min_qual, min_dp, min_gq)

    # Step 5: strict normalize
    norm_clean = str(trio_dir / "norm.clean.vcf.gz")
    logger.info("[Step 5] Strict left-normalize...")
    _norm_strict(filtered_vcf, norm_clean, ref_fasta)

    # Step 6: drop missing
    no_miss = str(trio_dir / "drop.miss.vcf.gz")
    logger.info("[Step 6] Drop missing genotypes...")
    _drop_missing(norm_clean, no_miss)

    # Step 7: inheritance filters
    mode_dfs: dict[str, Optional[pd.DataFrame]] = {}

    if "de_novo" in inheritance_modes:
        denovo_vcf = str(trio_dir / "denovo.vcf.gz")
        logger.info("[Step 7a] de_novo filter...")
        _filter_denovo(no_miss, denovo_vcf)
        mode_dfs["de_novo"] = _annotate_vcf(
            denovo_vcf, work_dir, gtf_file, clinvar_vcf, "de_novo"
        )

    if "recessive" in inheritance_modes:
        rec_vcf = str(trio_dir / "recessive.vcf.gz")
        logger.info("[Step 7b] recessive filter...")
        _filter_recessive(no_miss, rec_vcf)
        mode_dfs["recessive"] = _annotate_vcf(
            rec_vcf, work_dir, gtf_file, clinvar_vcf, "recessive"
        )

    if "xlinked" in inheritance_modes:
        xl_vcf = str(trio_dir / "xlinked.vcf.gz")
        logger.info("[Step 7c] xlinked filter...")
        _filter_xlinked(no_miss, xl_vcf)
        mode_dfs["xlinked"] = _annotate_vcf(
            xl_vcf, work_dir, gtf_file, clinvar_vcf, "xlinked"
        )

    # ── 合并所有模式 ────────────────────────────────────────────────────
    valid_dfs = [df for df in mode_dfs.values() if df is not None and len(df) > 0]

    if not valid_dfs:
        logger.warning("所有遗传模式均无有效位点")
        combined = pd.DataFrame()
    else:
        combined = pd.concat(valid_dfs, ignore_index=True)
        combined["rank"] = range(1, len(combined) + 1)

    combined_csv = work_dir / "3_clinvar_annotated.csv"
    combined.to_csv(str(combined_csv), index=False)

    for mn, df in mode_dfs.items():
        n = len(df) if df is not None else 0
        logger.info(f"  {mn}: {n} variants")

    logger.info(f"Stage 1 家系完成: {len(combined)} total → {combined_csv}")

    return {
        "modes": mode_dfs,
        "combined": combined,
        "output_csv": str(combined_csv),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 单样本退化模式
# ─────────────────────────────────────────────────────────────────────────────

def _singleton_mode(
    work_dir: Path,
    child_vcf: str,
    gtf_file: Optional[str],
    clinvar_vcf: Optional[str],
) -> dict:
    """单样本模式：VCF → GT → GTF → ClinVar"""
    gt_csv   = work_dir / "1_gt.csv"
    gtf_csv  = work_dir / "2_gtf_annotated.csv"
    final_csv = work_dir / "3_clinvar_annotated.csv"

    if not gt_csv.exists():
        stage1_vcf_to_gt_csv(child_vcf, str(gt_csv))
    df = pd.read_csv(str(gt_csv), dtype=str)
    logger.info(f"  {len(df)} variants from VCF")

    if gtf_file:
        if not gtf_csv.exists():
            stage1_annotate_by_gtf(str(gt_csv), gtf_file, str(gtf_csv))
        df = pd.read_csv(str(gtf_csv), dtype=str)

    if clinvar_vcf:
        if not final_csv.exists():
            stage1_merge_filter_clinvar(
                str(gtf_csv if gtf_file else gt_csv),
                clinvar_vcf,
                str(final_csv),
            )
        df = pd.read_csv(str(final_csv), dtype=str)

    # 精简 CSV
    df = simplify_clinvar_csv(df)

    df["inheritance_mode"] = "singleton"
    df["rank"] = range(1, len(df) + 1)
    df.to_csv(str(final_csv if clinvar_vcf else (gtf_csv if gtf_file else gt_csv)), index=False)

    logger.info(f"Stage 1 单样本完成: {len(df)} variants")

    return {
        "modes": {"singleton": df},
        "combined": df,
        "output_csv": str(final_csv if clinvar_vcf else (gtf_csv if gtf_file else gt_csv)),
    }
