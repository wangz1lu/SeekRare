"""
stage4_genos.py — Stage 4: Genos 模型分析

从 Stage 3 排序结果中选择感兴趣的位点，运行 Genos 全套分析：

Pipeline（4步）:
  Step 1: 1_build_variant_ref_csv.py
    输入: CHROM, POS, REF, ALT 列表 + 参考基因组 FASTA
    输出: variant_ref_input.csv（variant + reference 双序列）

  Step 2: 2_export_attention_by_variant_no_vcf.py
    输入: variant_ref_input.csv + Genos 模型路径
    输出: attention_result/ (每variant一个子目录)

  Step 3: 3_plot_variant_reference_logfc.py
    输入: attention_result/
    输出: 每variant的 log2FC 图 + CSV

  Step 4: 4_check_truth_peak.py
    输入: attention_result/
    输出: is_reasonable_truth (是否在正确位置有富集)

Usage:
    from seekrare.preprocess.stage4_genos import run_genos_analysis
    result = run_genos_analysis(
        sites=[("chr1", 123456, "A", "G"), ...],
        genome_fa="/path/to/GRCh38.fa",
        model_path="/path/to/Genos-1.2B",
        output_dir="/path/to/attention_result",
    )
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


# ─────────────────────────────────────────────────────────────────────────────
# 内部工具
# ─────────────────────────────────────────────────────────────────────────────

def _genos_scripts_dir() -> Path:
    """
    返回 scripts/genos/ 目录（包内自带的 Genos 工具脚本）。
    """
    import seekrare
    pkg_root = Path(seekrare.__file__).parent.parent
    return pkg_root / "scripts" / "genos"


def _script_path(name: str) -> str:
    """
    返回脚本的绝对路径。
    优先用包内脚本；必要时可扩展为支持部署机独立路径。
    """
    scripts_dir = _genos_scripts_dir()
    p = scripts_dir / name
    if p.exists():
        return str(p)
    raise FileNotFoundError(
        f"Cannot find Genos script '{name}' in {scripts_dir}"
    )


# ─────────────────────────────────────────────────────────────────────────────
# Step 1: build_variant_ref_csv
# ─────────────────────────────────────────────────────────────────────────────

def step1_build_variant_ref(
    sites: list[tuple],
    genome_fa: str,
    output_csv: str,
    flank: int = 2000,
) -> str:
    """
    Step 1: 从位点列表生成 variant reference CSV。

    Parameters
    ----------
    sites : list of (chrom, pos, ref, alt)
    genome_fa : 参考基因组 FASTA 路径
    output_csv : 输出 CSV 路径
    flank : 上下游扩展长度，默认 2000

    Returns
    -------
    str: 输出 CSV 路径
    """
    sites_csv = Path(output_csv).parent / "sites_input.csv"
    with open(sites_csv, "w", newline="") as f:
        for chrom, pos, ref, alt in sites:
            f.write(f"{chrom},{pos},{ref},{alt}\n")

    logger.info(f"Step 1: 构建 variant reference CSV ({len(sites)} sites)")

    cmd = [
        sys.executable,
        _script_path("1_build_variant_ref_csv.py"),
        "--sites_csv", str(sites_csv),
        "--fa", str(genome_fa),
        "--output_csv", str(output_csv),
        "--flank", str(flank),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Step 1 failed:\n{result.stderr}")

    logger.info(f"  → {output_csv}")
    return output_csv


# ─────────────────────────────────────────────────────────────────────────────
# Step 2: export_attention_by_variant
# ─────────────────────────────────────────────────────────────────────────────

def step2_export_attention(
    input_csv: str,
    model_path: str,
    output_dir: str,
    seq_chunk_size: int = 4096,
    seq_overlap: int = 1000,
    gpu: int = 0,
    max_variants: Optional[int] = None,
) -> str:
    """
    Step 2: 运行 Genos 模型，计算每 variant 的 attention matrix。

    Parameters
    ----------
    input_csv : step1 输出的 variant_ref_input.csv
    model_path : Genos 模型路径
    output_dir : attention 结果输出目录
    seq_chunk_size : 序列分块大小
    seq_overlap : 块重叠长度
    gpu : GPU 编号
    max_variants : 最大处理位点数（用于测试）

    Returns
    -------
    str: output_dir
    """
    logger.info(f"Step 2: Genos attention 分析 → {output_dir}")

    cmd = [
        sys.executable,
        _script_path("2_export_attention_by_variant_no_vcf.py"),
        "--input_csv", str(input_csv),
        "--model_path", str(model_path),
        "--output_dir", str(output_dir),
        "--seq_chunk_size", str(seq_chunk_size),
        "--seq_overlap", str(seq_overlap),
        "--gpu", str(gpu),
    ]
    if max_variants is not None:
        cmd.extend(["--max_variants", str(max_variants)])

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Step 2 failed:\n{result.stderr}")

    logger.info(f"  → {output_dir}")
    return output_dir


# ─────────────────────────────────────────────────────────────────────────────
# Step 3: plot log2FC
# ─────────────────────────────────────────────────────────────────────────────

def step3_plot_logfc(
    input_dir: str,
    max_variants: Optional[int] = None,
) -> str:
    """
    Step 3: 绘制 variant vs reference log2FC 图。

    Parameters
    ----------
    input_dir : step2 输出目录（attention_result/）
    max_variants : 最大处理位数（用于测试）

    Returns
    -------
    str: input_dir（原地修改）
    """
    logger.info(f"Step 3: 绘制 log2FC 图")

    cmd = [
        sys.executable,
        _script_path("3_plot_variant_reference_logfc.py"),
        "--input_dir", str(input_dir),
    ]
    if max_variants is not None:
        cmd.extend(["--max_variants", str(max_variants)])

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Step 3 failed:\n{result.stderr}")

    logger.info(f"  → logfc_plot_result/ 子目录已生成")
    return input_dir


# ─────────────────────────────────────────────────────────────────────────────
# Step 4: check truth peak
# ─────────────────────────────────────────────────────────────────────────────

def step4_check_peak(
    input_dir: str,
    output_csv: Optional[str] = None,
    window: int = 50,
    background_exclude: int = 500,
    top_quantile: float = 0.95,
    fold_change: float = 2.0,
    max_variants: Optional[int] = None,
) -> pd.DataFrame:
    """
    Step 4: 检验 variant 正确位置是否有 peak 富集。

    Parameters
    ----------
    input_dir : step2 输出目录
    output_csv : 总结果 CSV 路径（默认 input_dir/truth_peak_validation_all_summary.csv）
    window : 检验窗口大小（bp）
    background_exclude : 背景排除窗口
    top_quantile : global peak 阈值分位数
    fold_change : 富集倍数阈值
    max_variants : 最大处理位数

    Returns
    -------
    pd.DataFrame: 所有 variant 的验证结果
    """
    if output_csv is None:
        output_csv = str(Path(input_dir) / "truth_peak_validation_all_summary.csv")

    logger.info(f"Step 4: Peak 验证 → {output_csv}")

    cmd = [
        sys.executable,
        _script_path("4_check_truth_peak.py"),
        "--input_dir", str(input_dir),
        "--output_csv", str(output_csv),
        "--window", str(window),
        "--background_exclude", str(background_exclude),
        "--top_quantile", str(top_quantile),
        "--fold_change", str(fold_change),
    ]
    if max_variants is not None:
        cmd.extend(["--max_variants", str(max_variants)])

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Step 4 failed:\n{result.stderr}")

    df = pd.read_csv(output_csv)
    n_reasonable = int(df["is_reasonable_truth"].sum()) if "is_reasonable_truth" in df.columns else 0
    logger.info(f"  → {len(df)} variants validated, {n_reasonable}/{len(df)} passed")

    return df


# ─────────────────────────────────────────────────────────────────────────────
# 一键运行
# ─────────────────────────────────────────────────────────────────────────────

def run_genos_analysis(
    sites: list[tuple],
    genome_fa: str,
    model_path: str,
    output_dir: str,
    flank: int = 2000,
    seq_chunk_size: int = 4096,
    seq_overlap: int = 1000,
    gpu: int = 0,
    max_variants: Optional[int] = None,
    window: int = 50,
    background_exclude: int = 500,
    top_quantile: float = 0.95,
    fold_change: float = 2.0,
) -> pd.DataFrame:
    """
    Stage 4: Genos 全套分析（一键）。

    从 Stage 3 排序结果中选取感兴趣的位点，运行完整 Genos pipeline。

    Parameters
    ----------
    sites : list of (chrom, pos, ref, alt)
    genome_fa : str
        参考基因组 FASTA（GRCh38）
    model_path : str
        Genos 模型路径（部署机器上的绝对路径）
    output_dir : str
        输出根目录
    flank, seq_chunk_size, seq_overlap, gpu
        各 step 参数
    max_variants : int, optional
        最大处理位数（测试用）
    window, background_exclude, top_quantile, fold_change
        Step 4 检验参数

    Returns
    -------
    pd.DataFrame
        Step 4 汇总结果（含 is_reasonable_truth）
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Stage 4: Genos Analysis")
    logger.info("=" * 60)
    logger.info(f"  Sites: {len(sites)}")
    logger.info(f"  Genome: {genome_fa}")
    logger.info(f"  Model: {model_path}")
    logger.info(f"  Output: {output_dir}")

    # Step 1
    variant_ref_csv = output_dir / "variant_ref_input.csv"
    step1_build_variant_ref(
        sites=sites,
        genome_fa=genome_fa,
        output_csv=str(variant_ref_csv),
        flank=flank,
    )

    # Step 2
    attention_dir = output_dir / "attention_result"
    step2_export_attention(
        input_csv=str(variant_ref_csv),
        model_path=model_path,
        output_dir=str(attention_dir),
        seq_chunk_size=seq_chunk_size,
        seq_overlap=seq_overlap,
        gpu=gpu,
        max_variants=max_variants,
    )

    # Step 3
    step3_plot_logfc(
        input_dir=str(attention_dir),
        max_variants=max_variants,
    )

    # Step 4
    summary_csv = output_dir / "truth_peak_validation_all_summary.csv"
    result_df = step4_check_peak(
        input_dir=str(attention_dir),
        output_csv=str(summary_csv),
        window=window,
        background_exclude=background_exclude,
        top_quantile=top_quantile,
        fold_change=fold_change,
        max_variants=max_variants,
    )

    logger.info(f"Stage 4 完成: {output_dir}")
    return result_df


def read_sites_from_csv(
    csv_path: str,
    top_n: Optional[int] = None,
) -> list[tuple]:
    """
    从 Stage 3 输出的 CSV 读取位点。

    Parameters
    ----------
    csv_path : str
        Stage 3 排序后的 CSV 路径
    top_n : int, optional
        只取前 N 行

    Returns
    -------
    list of (chrom, pos, ref, alt)
    """
    df = pd.read_csv(csv_path)
    if top_n:
        df = df.head(top_n)

    sites = []
    for _, row in df.iterrows():
        try:
            chrom = str(row.get("CHROM", row.get("chrom", ""))).strip()
            pos = int(float(row.get("POS", row.get("pos", 0))))
            ref = str(row.get("REF", row.get("ref", ""))).strip()
            alt = str(row.get("ALT", row.get("alt", ""))).strip()
            if chrom and pos and ref and alt:
                sites.append((chrom, pos, ref, alt))
        except Exception:
            continue

    logger.info(f"从 {csv_path} 读取 {len(sites)} 个位点")
    return sites


# ─────────────────────────────────────────────────────────────────────────────
# pipeline.py 调用的 wrapper（与旧接口兼容）
# ─────────────────────────────────────────────────────────────────────────────

def stage4_genos_analysis(
    stage3_csv: str,
    genome_fa: str,
    model_path: str,
    output_dir: str,
    top_n: int = 10,
    gpu: int = 0,
    flank: int = 2000,
    seq_chunk_size: int = 4096,
    seq_overlap: int = 1000,
    window: int = 50,
    background_exclude: int = 500,
    top_quantile: float = 0.95,
    fold_change: float = 2.0,
) -> pd.DataFrame:
    """
    Stage 4 Genos wrapper（供 pipeline.py 调用）。

    从 Stage 3 CSV 取 top-N 位点，运行 Genos 全套分析。
    """
    sites = read_sites_from_csv(stage3_csv, top_n=top_n)
    if not sites:
        raise ValueError(f"Stage 3 CSV 中未找到有效位点: {stage3_csv}")

    return run_genos_analysis(
        sites=sites,
        genome_fa=genome_fa,
        model_path=model_path,
        output_dir=output_dir,
        flank=flank,
        seq_chunk_size=seq_chunk_size,
        seq_overlap=seq_overlap,
        gpu=gpu,
        max_variants=None,
        window=window,
        background_exclude=background_exclude,
        top_quantile=top_quantile,
        fold_change=fold_change,
    )