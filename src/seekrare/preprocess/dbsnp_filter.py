"""
dbsnp_filter.py — dbSNP common 变异过滤

直接执行 bcftools view -T 过滤。

过滤前检查 input VCF 和 dbSNP 的染色体格式是否兼容：
  - 预期格式：input VCF 和 dbSNP 均有 chr 前缀（chr1, chr2, ...）
  - 若格式不一致（一方有 chr 前缀，另一方没有），直接报错退出，
    提示用户自行用 bcftools annotate --rename-chrs 统一格式后再来。
"""

from __future__ import annotations

import gzip
import subprocess
from pathlib import Path
from typing import Union


def _opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def detect_chrom_format(vcf_path: str) -> str:
    """检测是否有 chr 前缀。返回 'chr' 或 'nochr'。"""
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
    input_vcf : str   输入 VCF.gz（预期有 chr 前缀）
    dbsnp_vcf : str   dbSNP common VCF.gz（也必须有 chr 前缀）
    output_vcf : str  输出 VCF.gz

    Returns dict with n_before, n_after, n_removed, input_format, dbsnp_format
    """
    input_vcf, dbsnp_vcf, output_vcf = str(input_vcf), str(dbsnp_vcf), str(output_vcf)

    input_fmt = detect_chrom_format(input_vcf)
    dbsnp_fmt = detect_chrom_format(dbsnp_vcf)

    # ── 格式检查 ─────────────────────────────────────────────────────────
    if input_fmt != dbsnp_fmt:
        msg = (
            f"染色体格式不一致，过滤失败！\n"
            f"  Input VCF format:  '{'chr' if input_fmt == 'chr' else 'no chr'} prefix'  (示例: {input_fmt})\n"
            f"  dbSNP VCF format:  '{'chr' if dbsnp_fmt == 'chr' else 'no chr'} prefix'\n"
            f"\n"
            f"请先用 bcftools annotate --rename-chrs 统一 dbSNP 的染色体格式：\n"
            f"\n"
            f"  for i in {{1..22}} X Y; do echo \"$i chr$i\"; done > rename.txt\n"
            f"  bcftools annotate --rename-chrs rename.txt /path/to/your_dbSNP.vcf.gz -Oz -o /path/to/dbSNP_chr.vcf.gz\n"
            f"\n"
            f"然后用处理后的 dbSNP_vcf 重新运行。\n"
        )
        raise ValueError(msg)

    # ── 直接过滤 ─────────────────────────────────────────────────────────
    n_before = count_vcf_variants(input_vcf)

    r = subprocess.run(
        ["bcftools", "view",
         "-T", f"^{dbsnp_vcf}",
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
        "input_format": input_fmt,
        "dbsnp_format": dbsnp_fmt,
    }
