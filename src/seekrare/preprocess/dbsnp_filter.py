"""
dbsnp_filter.py — dbSNP common 变异过滤

自动检测 VCF 和 dbSNP 的染色体命名格式（chr1 vs 1），
自动统一后执行 bcftools view -T 排除。
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from pathlib import Path
from typing import Union


def _opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def detect_chrom_format(vcf_path: str) -> str:
    """检测 VCF 染色体格式: 'chr' (带前缀) 或 'nochr' (不带前缀)。"""
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


def get_chromosomes(vcf_path: str, max_lines: int = 200) -> set[str]:
    """提取 VCF 中所有唯一染色体名称。"""
    chroms = set()
    with _opener(vcf_path) as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            chroms.add(line.split("\t")[0])
            if i > max_lines:
                break
    return chroms


def count_vcf_variants(vcf_path: str) -> int:
    """统计 VCF 非头部行数。"""
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
    执行 dbSNP common 变异过滤（自动处理染色体命名不一致）。

    Parameters
    ----------
    input_vcf : str
        输入 VCF.gz
    dbsnp_vcf : str
        dbSNP common VCF.gz
    output_vcf : str
        输出 VCF.gz（剔除 common 后）

    Returns
    -------
    dict: {"n_before", "n_after", "n_removed", "input_format", "dbsnp_format"}
    """
    input_vcf, dbsnp_vcf, output_vcf = str(input_vcf), str(dbsnp_vcf), str(output_vcf)
    work_dir = Path(output_vcf).parent

    # ── 1. 检测染色体命名格式 ──────────────────────────────────────────────
    input_fmt = detect_chrom_format(input_vcf)
    dbsnp_fmt = detect_chrom_format(dbsnp_vcf)

    # ── 2. 统一 dbSNP 格式以匹配输入 VCF ──────────────────────────────────
    if input_fmt != dbsnp_fmt:
        chroms = get_chromosomes(dbsnp_vcf)
        rename_txt = work_dir / "dbsnp_chr_rename.txt"
        with open(rename_txt, "w") as f:
            for chrom in sorted(chroms):
                if input_fmt == "chr" and not chrom.startswith("chr"):
                    f.write(f"{chrom}\tchr{chrom}\n")
                elif input_fmt == "nochr" and chrom.startswith("chr"):
                    f.write(f"{chrom}\t{chrom[3:]}\n")

        dbsnp_renamed = work_dir / "dbsnp_renamed.vcf.gz"
        r = subprocess.run(
            ["bcftools", "annotate", "--rename-chrs", str(rename_txt),
             "-Oz", "-o", str(dbsnp_renamed), dbsnp_vcf],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            print(f"[WARN] bcftools annotate rename-chrs failed: {r.stderr[:200]}", file=sys.stderr)
            dbsnp_for_filter = dbsnp_vcf
        else:
            subprocess.run(["bcftools", "index", "-f", str(dbsnp_renamed)], capture_output=True)
            dbsnp_for_filter = str(dbsnp_renamed)
    else:
        dbsnp_for_filter = dbsnp_vcf

    # ── 3. 统计过滤前 ───────────────────────────────────────────────────────
    n_before = count_vcf_variants(input_vcf)

    # ── 4. bcftools view -T 排除 dbSNP ─────────────────────────────────────
    r = subprocess.run(
        ["bcftools", "view", "-T", f"^{dbsnp_for_filter}",
         "-Oz", "-o", output_vcf, input_vcf],
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


def main():
    parser = argparse.ArgumentParser(description="dbSNP common 变异过滤")
    parser.add_argument("--input", required=True)
    parser.add_argument("--dbsnp", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = run_dbsnp_filter(args.input, args.dbsnp, args.output)
    print(
        f"dbSNP filter: {result['n_removed']}/{result['n_before']} removed "
        f"(input={result['input_format']}, dbsnp={result['dbsnp_format']})"
    )
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
