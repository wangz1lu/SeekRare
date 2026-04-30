"""
dbsnp_filter.py — dbSNP common 变异过滤

策略：以 dbSNP 的染色体格式为基准。
  - input VCF 有 chr 前缀，dbSNP 没有 → 给 input VCF 的 CHROM 列去掉 chr 前缀
  - input VCF 无 chr 前缀，dbSNP 有 → 给 input VCF 的 CHROM 列加上 chr 前缀
  - 然后 bcftools view -T ^dbSNP 过滤

使用 bcftools reheader -c 重命名 CHROM 列（不修改 header）。
bcftools -T 过滤时只需要数据行 CHROM 匹配，不强制要求 header contig 定义一致。
"""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path
from typing import Union


STANDARD_CHROMS = {
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22",
    "X","Y","MT","M"
}


def _opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def detect_chrom_format(vcf_path: str) -> str:
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


def get_chromosomes_from_lines(vcf_path: str, max_lines: int = 500) -> list[str]:
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
            if len(seen) >= 60:
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
    """写 reheader -c 映射文件 (from → to)。"""
    with open(out_path, "w") as f:
        for chrom in chroms:
            base = chrom.lstrip("chr")
            if base not in STANDARD_CHROMS:
                continue
            if from_fmt == "chr" and to_fmt == "nochr":
                f.write(f"{chrom}\t{base}\n")
            elif from_fmt == "nochr" and to_fmt == "chr":
                f.write(f"{chrom}\tchr{base}\n")


def _reheader_chroms(vcf_path: str, work_dir: Path, from_fmt: str, to_fmt: str) -> str:
    """
    用 bcftools reheader -c 重命名 CHROM 列。
    返回重命名后的 VCF 路径。
    """
    chroms = get_chromosomes_from_lines(vcf_path)
    chr_map = work_dir / "input_chr_map.txt"
    _write_chr_map(chroms, from_fmt, to_fmt, chr_map)

    output = work_dir / f"input_renamed.{Path(vcf_path).name}"

    r = subprocess.run(
        ["bcftools", "reheader",
         "-- chromosomes", str(chr_map),   # -c FILE: rename CHROM column
         "-o", str(output),
         vcf_path],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"bcftools reheader -c failed: {r.stderr}")

    subprocess.run(["bcftools", "index", "-f", str(output)], capture_output=True)
    return str(output)


def run_dbsnp_filter(
    input_vcf: Union[str, Path],
    dbsnp_vcf: Union[str, Path],
    output_vcf: Union[str, Path],
) -> dict:
    """
    执行 dbSNP common 过滤。

    Returns dict with n_before, n_after, n_removed, input_format, dbsnp_format
    """
    input_vcf, dbsnp_vcf, output_vcf = str(input_vcf), str(dbsnp_vcf), str(output_vcf)
    work_dir = Path(output_vcf).parent

    input_fmt = detect_chrom_format(input_vcf)
    dbsnp_fmt = detect_chrom_format(dbsnp_vcf)

    input_for_filter = input_vcf

    if input_fmt != dbsnp_fmt:
        print(f"  Format mismatch: input={input_fmt}, dbsnp={dbsnp_fmt}")
        print(f"  Converting input VCF CHROM column to '{dbsnp_fmt}' format via bcftools reheader -c...")
        try:
            input_for_filter = _reheader_chroms(input_vcf, work_dir, input_fmt, dbsnp_fmt)
            print(f"  Input VCF converted: {input_for_filter}")
        except RuntimeError as e:
            print(f"  [ERROR] {e}", file=sys.stderr)
            raise

    n_before = count_vcf_variants(input_for_filter)

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

    return {
        "n_before": n_before,
        "n_after": n_after,
        "n_removed": n_before - n_after,
        "input_format": input_fmt,
        "dbsnp_format": dbsnp_fmt,
    }
