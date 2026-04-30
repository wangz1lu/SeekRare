"""
dbsnp_filter.py — dbSNP common 变异过滤

策略：
  1. 生成 chr 重命名映射: 1→chr1, 2→chr2, ..., X→chrX, Y→chrY
  2. bcftools annotate --rename-chrs rename.txt dbSNP → dbSNP_chr.vcf.gz
  3. bcftools view -T ^dbSNP_chr.vcf.gz input.vcf.gz → output.vcf.gz

输入 VCF 和 dbSNP 都不修改，只操作 dbSNP 的副本。
"""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path
from typing import Union


STANDARD_CHROMS = [
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22",
    "X","Y","MT","M"
]


def _opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def count_vcf_variants(vcf_path: str) -> int:
    count = 0
    with _opener(vcf_path) as f:
        for line in f:
            if not line.startswith("#"):
                count += 1
    return count


def run_dbsnp_filter(
    input_vcf: Union[str, Path],
    dbsnp_vcf: Union[str, Path],
    output_vcf: Union[str, Path],
) -> dict:
    """
    执行 dbSNP common 过滤。

    Parameters
    ----------
    input_vcf : str   原始 VCF.gz
    dbsnp_vcf : str   dbSNP common VCF.gz（原始文件，不修改）
    output_vcf : str  输出 VCF.gz

    Returns dict with n_before, n_after, n_removed, input_format, dbsnp_format
    """
    input_vcf, dbsnp_vcf, output_vcf = str(input_vcf), str(dbsnp_vcf), str(output_vcf)
    work_dir = Path(output_vcf).parent

    # ── Step 1: 生成 chr 重命名映射文件 ─────────────────────────────────
    # 格式: "1\tchr1" （和用户给的脚本完全一致）
    rename_txt = work_dir / "chr_rename.txt"
    with open(rename_txt, "w") as f:
        for chrom in STANDARD_CHROMS:
            f.write(f"{chrom}\tchr{chrom}\n")

    # ── Step 2: 给 dbSNP 的 CHROM 列加上 chr 前缀 ─────────────────────
    # 输出到工作目录，不修改原始 dbSNP
    dbsnp_chr = work_dir / "dbsnp_chr.vcf.gz"

    if not dbsnp_chr.exists():
        print(f"  Step 2: bcftools annotate --rename-chrs dbSNP...")
        r = subprocess.run(
            ["bcftools", "annotate",
             "--rename-chrs", str(rename_txt),
             "-Oz", "-o", str(dbsnp_chr),
             dbsnp_vcf],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            # annotate 失败（可能是 header 问题），打印警告但继续
            stderr = r.stderr.strip()
            print(f"  [WARN] annotate --rename-chrs: {stderr[:200]}", file=sys.stderr)
            print(f"  Trying with --force flag...", file=sys.stderr)
            r = subprocess.run(
                ["bcftools", "annotate",
                 "--rename-chrs", str(rename_txt),
                 "--force",
                 "-Oz", "-o", str(dbsnp_chr),
                 dbsnp_vcf],
                capture_output=True, text=True,
            )
            if r.returncode != 0:
                raise RuntimeError(f"bcftools annotate --rename-chrs --force failed: {r.stderr}")

        subprocess.run(["bcftools", "index", "-f", str(dbsnp_chr)], capture_output=True)
        print(f"  dbSNP chr version: {dbsnp_chr}")
    else:
        print(f"  [跳过] dbSNP chr 版本已存在: {dbsnp_chr}")

    # ── Step 3: 排除 dbSNP common 变异 ──────────────────────────────────
    print(f"  Step 3: bcftools view -T ^{dbsnp_vcf} input.vcf.gz...")
    n_before = count_vcf_variants(input_vcf)

    r = subprocess.run(
        ["bcftools", "view",
         "-T", f"^{dbsnp_chr}",
         "-Oz", "-o", output_vcf,
         input_vcf],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"bcftools view -T failed: {r.stderr}")

    n_after = count_vcf_variants(output_vcf)

    return {
        "n_before": n_before,
        "n_after": n_after,
        "n_removed": n_before - n_after,
        "input_format": "chr",
        "dbsnp_format": "nochr",
    }
