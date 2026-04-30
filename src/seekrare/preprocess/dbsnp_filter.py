"""
dbsnp_filter.py — dbSNP common 变异过滤（稳定版）

策略：以 dbSNP 的染色体格式为基准。
  - 若 input VCF 有 chr 前缀但 dbSNP 没有 → 用 bcftools reheader -c
    把 input 的 CHROM 列批量去掉 chr 前缀
  - 若 input VCF 无 chr 但 dbSNP 有 → 用 bcftools reheader -c 给 input 加上 chr 前缀
  - 然后 bcftools view -T ^dbSNP 排除

不需要修改 dbSNP 文件本身。
"""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path
from typing import Union


def _opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def detect_chrom_format(vcf_path: str) -> str:
    """检测 'chr' 前缀 或 'nochr'。"""
    chroms = set()
    with _opener(vcf_path) as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            chroms.add(line.split("\t")[0])
            if len(chroms) >= 50 or i > 200:
                break
    if not chroms:
        return "nochr"
    chr_count = sum(1 for c in chroms if str(c).startswith("chr"))
    return "chr" if chr_count >= len(chroms) * 0.5 else "nochr"


def get_chromosomes(vcf_path: str, max_lines: int = 500) -> list[str]:
    """从 VCF 数据行提取所有唯一染色体名称（有序）。"""
    chroms = []
    seen = set()
    with _opener(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")[0]
            if c not in seen:
                seen.add(c)
                chroms.append(c)
            if len(seen) >= 60 or len(chroms) >= max_lines:
                break
    return chroms


def count_vcf_variants(vcf_path: str) -> int:
    count = 0
    with _opener(vcf_path) as f:
        for line in f:
            if not line.startswith("#"):
                count += 1
    return count


def _write_chr_map(chroms: list[str], from_fmt: str, to_fmt: str, out_path: Path):
    """写 chr 映射文件 (old → new)。"""
    with open(out_path, "w") as f:
        for chrom in chroms:
            if from_fmt == "chr" and chrom.startswith("chr") and to_fmt == "nochr":
                # chr1 → 1
                f.write(f"{chrom}\t{chrom[3:]}\n")
            elif from_fmt == "nochr" and not chrom.startswith("chr") and to_fmt == "chr":
                # 1 → chr1
                f.write(f"{chrom}\tchr{chrom}\n")


def run_dbsnp_filter(
    input_vcf: Union[str, Path],
    dbsnp_vcf: Union[str, Path],
    output_vcf: Union[str, Path],
) -> dict:
    """
    执行 dbSNP common 过滤，自动统一染色体格式。

    Parameters
    ----------
    input_vcf, dbsnp_vcf, output_vcf : str

    Returns dict with n_before, n_after, n_removed, input_format, dbsnp_format
    """
    input_vcf, dbsnp_vcf, output_vcf = str(input_vcf), str(dbsnp_vcf), str(output_vcf)
    work_dir = Path(output_vcf).parent

    # ── 1. 检测格式 ───────────────────────────────────────────────────────
    input_fmt = detect_chrom_format(input_vcf)
    dbsnp_fmt = detect_chrom_format(dbsnp_vcf)

    # ── 2. 统一 input VCF 到 dbSNP 的格式 ──────────────────────────────
    input_for_filter = input_vcf

    if input_fmt != dbsnp_fmt:
        chroms = get_chromosomes(input_vcf)
        chr_map = work_dir / "input_chr_map.txt"
        _write_chr_map(chroms, input_fmt, dbsnp_fmt, chr_map)

        renamed = work_dir / "input_renamed.vcf.gz"

        # bcftools reheader -c: old→new 两列文件，重命名 CHROM 列
        r = subprocess.run(
            ["bcftools", "reheader",
             "-- chromosomes", str(chr_map),   # -c FILE
             "-o", str(renamed),
             input_vcf],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            print(f"  [WARN] bcftools reheader failed: {r.stderr[:300]}", file=sys.stderr)
            raise RuntimeError(f"bcftools reheader -c failed: {r.stderr}")

        subprocess.run(["bcftools", "index", "-f", str(renamed)], capture_output=True)
        input_for_filter = str(renamed)
        print(f"  Converted input {input_fmt}→{dbsnp_fmt}: {renamed}")

    # ── 3. 统计过滤前 ─────────────────────────────────────────────────────
    n_before = count_vcf_variants(input_for_filter)

    # ── 4. 排除 dbSNP common ───────────────────────────────────────────────
    r = subprocess.run(
        ["bcftools", "view",
         "-T", f"^{dbsnp_vcf}",
         "-Oz", "-o", output_vcf,
         input_for_filter],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"bcftools view -T failed: {r.stderr}")

    n_after = count_vcf_variants(output_vcf)
    n_removed = n_before - n_after

    return {
        "n_before": n_before,
        "n_after": n_after,
        "n_removed": n_removed,
        "input_format": input_fmt,
        "dbsnp_format": dbsnp_fmt,
    }
