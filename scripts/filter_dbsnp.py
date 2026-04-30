#!/usr/bin/env python3
"""
filter_dbsnp.py — 安全剔除 dbSNP common 变异

自动检测 VCF 和 dbSNP 的染色体命名格式（chr1 vs 1），
自动统一后再执行 bcftools view -T 排除。

用法:
    python filter_dbsnp.py \
        --input /path/to/input.vcf.gz \
        --dbsnp /path/to/dbsnp_common.vcf.gz \
        --output /path/to/output.nocommon.vcf.gz

逻辑:
    1. 读取两个 VCF header，检测染色体命名格式
    2. 如果不一致，自动生成 chr 转换表
    3. 统一 dbSNP 的 chr 前缀（与你 VCF 格式对齐）
    4. bcftools view -T ^dbSNP_common.vcf.gz 执行排除
    5. 报告过滤前后行数
"""

import argparse
import gzip
import re
import subprocess
import sys
from pathlib import Path


def detect_chrom_format(vcf_path: str) -> str:
    """
    检测 VCF 的染色体命名格式。

    扫描前 1000 行非头部的 CHROM 值，判断是 'chr1' 还是 '1'。
    返回 "chr" (带前缀) 或 "nochr" (不带前缀)。
    """
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    chroms = set()

    with opener(vcf_path, "rt") as f:
        for i, line in enumerate(f):
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chroms.add(cols[0])
            if len(chroms) >= 50 or i > 1000:
                break

    chr_count = sum(1 for c in chroms if str(c).startswith("chr"))
    total = len(chroms)

    if total == 0:
        return "nochr"  # safe default

    # 多数原则
    if chr_count >= total * 0.5:
        return "chr"
    return "nochr"


def build_rename_txt(chrom_list: list[str], target_format: str, source_format: str) -> list[tuple[str, str]]:
    """
    生成染色体重命名列表 (from → to)。

    target_format / source_format: "chr" 或 "nochr"
    """
    rename = []
    for chrom in chrom_list:
        if target_format == "chr" and not chrom.startswith("chr"):
            rename.append((chrom, f"chr{chrom}"))
        elif target_format == "nochr" and chrom.startswith("chr"):
            rename.append((chrom, chrom[3:]))  # 去掉 chr
    return rename


def get_chromosomes_from_vcf(vcf_path: str) -> set[str]:
    """从 VCF 中提取所有唯一的染色体名称。"""
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    chroms = set()
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chroms.add(cols[0])
            if len(chroms) >= 100:
                break
    return chroms


def run_bcftools_annotate_rename(input_vcf: str, output_vcf: str, rename_map: list[tuple[str, str]]) -> bool:
    """用 bcftools annotate --rename-chrs 重命名染色体。"""
    if not rename_map:
        return True  # 不需要重命名

    # 写入 rename.txt
    rename_txt = Path(output_vcf).parent / "chr_rename.txt"
    with open(rename_txt, "w") as f:
        for old, new in rename_map:
            f.write(f"{old}\t{new}\n")

    cmd = [
        "bcftools", "annotate",
        "--rename-chrs", str(rename_txt),
        "-Oz", "-o", output_vcf,
        input_vcf,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    bcftools annotate rename-chrs failed: {result.stderr[:500]}", file=sys.stderr)
        return False
    # 索引
    subprocess.run(["bcftools", "index", "-f", output_vcf], capture_output=True)
    return True


def run_filter(input_vcf: str, dbsnp_vcf: str, output_vcf: str) -> dict:
    """
    执行 dbSNP common 变异过滤。

    Returns
    -------
    dict: {"n_before": int, "n_after": int, "n_removed": int, "dbsnp_format": str, "input_format": str}
    """
    input_vcf = str(input_vcf)
    dbsnp_vcf = str(dbsnp_vcf)
    output_vcf = str(output_vcf)
    work_dir = Path(output_vcf).parent

    print(f"Input VCF:       {input_vcf}")
    print(f"dbSNP VCF:       {dbsnp_vcf}")
    print(f"Output VCF:      {output_vcf}")

    # ── 1. 检测染色体命名格式 ──────────────────────────────────────────────
    input_format = detect_chrom_format(input_vcf)
    dbsnp_format = detect_chrom_format(dbsnp_vcf)
    print(f"\nChromosome naming detected:")
    print(f"  Input VCF:  '{'chr' if input_format == 'chr' else 'no chr'} prefix' ({input_format})")
    print(f"  dbSNP VCF:  '{'chr' if dbsnp_format == 'chr' else 'no chr'} prefix' ({dbsnp_format})")

    # ── 2. 统一 dbSNP 格式 ─────────────────────────────────────────────────
    dbsnp_renamed = work_dir / "dbsnp_renamed.vcf.gz"

    if input_format != dbsnp_format:
        print(f"\nRenaming dbSNP chromosomes to match input VCF format...")
        chroms = get_chromosomes_from_vcf(dbsnp_vcf)

        if input_format == "chr":
            # dbSNP 没有 chr，加上
            target = "chr1"; source = "1"
        else:
            # dbSNP 有 chr，去掉
            target = "1"; source = "chr1"

        rename_list = [(c, c[3:] if c.startswith("chr") else f"chr{c}")
                        for c in chroms if (c.startswith("chr") and input_format == "nochr")
                        or (not c.startswith("chr") and input_format == "chr")]

        # 写 rename.txt: from \t to
        rename_txt = work_dir / "dbsnp_chr_rename.txt"
        with open(rename_txt, "w") as f:
            for old_name, new_name in rename_list:
                f.write(f"{old_name}\t{new_name}\n")

        cmd = [
            "bcftools", "annotate",
            "--rename-chrs", str(rename_txt),
            "-Oz", "-o", str(dbsnp_renamed),
            dbsnp_vcf,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  bcftools rename failed: {result.stderr[:300]}", file=sys.stderr)
            print("  Trying without renaming (using bcftools -T anyway)...", file=sys.stderr)
            dbsnp_for_filter = dbsnp_vcf
        else:
            subprocess.run(["bcftools", "index", "-f", str(dbsnp_renamed)], capture_output=True)
            dbsnp_for_filter = str(dbsnp_renamed)
            print(f"  dbSNP renamed: {dbsnp_renamed}")
    else:
        print(f"\nChromosome formats match — no renaming needed")
        dbsnp_for_filter = dbsnp_vcf

    # ── 3. 统计过滤前行数 ──────────────────────────────────────────────────
    n_before = count_vcf_lines(input_vcf)

    # ── 4. bcftools view -T 排除 ──────────────────────────────────────────
    print(f"\nRunning bcftools view -T ^{dbsnp_vcf}...")
    cmd = [
        "bcftools", "view",
        "-T", f"^{dbsnp_for_filter}",
        "-Oz", "-o", output_vcf,
        input_vcf,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  FAILED: {result.stderr[:500]}", file=sys.stderr)
        raise RuntimeError(f"bcftools view -T failed: {result.stderr}")

    # ── 5. 统计过滤后行数 ──────────────────────────────────────────────────
    n_after = count_vcf_lines(output_vcf)
    n_removed = n_before - n_after

    print(f"\nResults:")
    print(f"  Before:  {n_before:,} variants")
    print(f"  After:   {n_after:,} variants")
    print(f"  Removed: {n_removed:,} common variants ({n_removed/n_before*100:.1f}%)")
    print(f"  Output:  {output_vcf}")

    return {
        "n_before": n_before,
        "n_after": n_after,
        "n_removed": n_removed,
        "dbsnp_format": dbsnp_format,
        "input_format": input_format,
    }


def count_vcf_lines(vcf_path: str) -> int:
    """统计 VCF 非头部行数。"""
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    count = 0
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description="安全剔除 dbSNP common 变异")
    parser.add_argument("--input", required=True, help="输入 VCF.gz")
    parser.add_argument("--dbsnp", required=True, help="dbSNP common VCF.gz")
    parser.add_argument("--output", required=True, help="输出 VCF.gz（剔除 common 后）")
    args = parser.parse_args()

    result = run_filter(args.input, args.dbsnp, args.output)
    print(f"\n✅ dbSNP filtering complete")
    return result


if __name__ == "__main__":
    main()
