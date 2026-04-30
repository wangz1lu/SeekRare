"""
dbsnp_filter.py — dbSNP common 变异过滤（稳定版）

策略：以 input VCF 格式为基准。
  - 若 input 有 chr 前缀，dbSNP 没有 → 给 dbSNP 加上 chr 前缀（优先用 bcftools annotate --rename-chrs）
  - 若 input 没有 chr 前缀，dbSNP 有 → 给 dbSNP 去掉 chr 前缀（bcftools reheader -c）
  - 然后 bcftools view -T ^dbSNP 过滤

不需要修改 input VCF，只改 dbSNP。
"""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path
from typing import Union


STANDARD_CHROMS = ["1","2","3","4","5","6","7","8","9","10",
                   "11","12","13","14","15","16","17","18","19","20","21","22",
                   "X","Y","MT","M"]


def _opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def detect_chrom_format(vcf_path: str) -> str:
    """返回 'chr' 或 'nochr'。"""
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


def _write_rename_txt(chroms: list[str], from_fmt: str, to_fmt: str, out_path: Path):
    """写 chr 映射文件 (old → new)。"""
    with open(out_path, "w") as f:
        for chrom in chroms:
            # 只处理标准染色体
            base = chrom.lstrip("chr")
            if base not in STANDARD_CHROMS:
                continue
            if from_fmt == "nochr" and to_fmt == "chr":
                # 1 → chr1
                f.write(f"{chrom}\tchr{base}\n")
            elif from_fmt == "chr" and to_fmt == "nochr":
                # chr1 → 1
                f.write(f"{chrom}\t{base}\n")


def _rename_vcf_chroms(vcf_path: str, work_dir: Path, from_fmt: str, to_fmt: str) -> str:
    """
    用 bcftools annotate --rename-chrs 给 VCF 的染色体加上/去掉前缀。
    返回重命名后的 VCF 路径。
    """
    chroms = get_chromosomes_from_lines(vcf_path)
    rename_txt = work_dir / "chr_rename.txt"
    _write_rename_txt(chroms, from_fmt, to_fmt, rename_txt)

    output = work_dir / f"renamed_{Path(vcf_path).name}"
    r = subprocess.run(
        ["bcftools", "annotate",
         "--rename-chrs", str(rename_txt),
         "-Oz", "-o", str(output),
         vcf_path],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"bcftools annotate --rename-chrs failed: {r.stderr}")

    subprocess.run(["bcftools", "index", "-f", str(output)], capture_output=True)
    return str(output)


def run_dbsnp_filter(
    input_vcf: Union[str, Path],
    dbsnp_vcf: Union[str, Path],
    output_vcf: Union[str, Path],
) -> dict:
    """
    执行 dbSNP common 过滤。

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

    # ── 2. 统一 dbSNP 到 input 的格式 ────────────────────────────────────
    # 永远不动 input VCF，只改 dbSNP
    dbsnp_for_filter = dbsnp_vcf

    if input_fmt != dbsnp_fmt:
        print(f"  Chromosome format: input={input_fmt}, dbsnp={dbsnp_fmt}")
        print(f"  Converting dbSNP to '{input_fmt}' format...")

        try:
            dbsnp_for_filter = _rename_vcf_chroms(dbsnp_vcf, work_dir, dbsnp_fmt, input_fmt)
            print(f"  dbSNP converted → {dbsnp_for_filter}")
        except RuntimeError as e:
            print(f"  [WARN] annotate --rename-chrs failed: {e}", file=sys.stderr)
            raise

    # ── 3. 统计过滤前 ─────────────────────────────────────────────────────
    n_before = count_vcf_variants(input_vcf)

    # ── 4. bcftools view -T 排除 dbSNP ────────────────────────────────────
    r = subprocess.run(
        ["bcftools", "view",
         "-T", f"^{dbsnp_for_filter}",
         "-Oz", "-o", output_vcf,
         input_vcf],
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
