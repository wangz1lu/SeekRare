#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re
import os
import sys
import time
from collections import defaultdict

if len(sys.argv) != 4:
    print("""
用法:
  python annotate_vcf_csv_by_ncbi_gtf.py <input_csv> <gtf_file> <output_csv>

示例:
  python annotate_vcf_csv_by_ncbi_gtf.py \\
    /mnt/work-dir/changan/vcf2csv/output_csv_test.csv \\
    /mnt/work-dir/changan/ncbi_dataset/data/GCF_000001405.26/genomic.gtf \\
    /mnt/work-dir/changan/vcf2csv/output_csv_test_annotated.csv
""")
    sys.exit(1)

INPUT_CSV = sys.argv[1]
GTF_FILE = sys.argv[2]
OUTPUT_CSV = sys.argv[3]

out_dir = os.path.dirname(OUTPUT_CSV)
if out_dir:
    os.makedirs(out_dir, exist_ok=True)

# ============================================================
# 1. NCBI accession -> chr 映射
# ============================================================

NCBI_TO_CHR = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9": "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9": "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
    "NC_012920.1": "chrM",
}

FEATURE_PRIORITY = {
    "stop_codon": 0,
    "start_codon": 1,
    "CDS": 2,
    "exon": 3,
    "transcript": 4,
    "gene": 5,
}

FEATURE_ORDER = ["stop_codon", "start_codon", "CDS", "exon", "transcript", "gene"]

def normalize_csv_chrom(chrom):
    """
    CSV里面：
      chr11_KI270721v1_random -> chr11
      chr1 -> chr1
      chrX -> chrX
    """
    chrom = str(chrom).strip()

    m = re.match(r"^(chr[0-9XYM]+)", chrom)
    if m:
        return m.group(1)

    return chrom

def normalize_gtf_chrom(chrom):
    """
    GTF里面：
      NC_000001.11 -> chr1
      NC_000023.11 -> chrX
    """
    chrom = str(chrom).strip()
    return NCBI_TO_CHR.get(chrom, chrom)

def parse_gene_name(attr):
    patterns = [
        r'gene_name\s+"([^"]+)"',
        r'gene\s+"([^"]+)"',
        r'Name\s+"([^"]+)"',
        r'gene_id\s+"([^"]+)"',
    ]

    for p in patterns:
        m = re.search(p, attr)
        if m:
            return m.group(1)

    return ""

def feature_rank(feature):
    return FEATURE_PRIORITY.get(feature, 99)

# ============================================================
# 2. 读取 CSV
# ============================================================

print("读取输入CSV...")
snp = pd.read_csv(INPUT_CSV)
snp.columns = [c.lstrip("#") for c in snp.columns]

if "CHROM" not in snp.columns or "POS" not in snp.columns:
    raise ValueError("输入CSV必须包含 CHROM 和 POS 两列")

snp["CHROM"] = snp["CHROM"].astype(str).str.strip()
snp["POS"] = pd.to_numeric(snp["POS"], errors="coerce")

# 用于匹配GTF的染色体
snp["_CHROM_KEY"] = snp["CHROM"].apply(normalize_csv_chrom)

print(f"输入位点数: {len(snp):,}")
print("CSV染色体示例:", snp["_CHROM_KEY"].unique()[:10])

# ============================================================
# 3. 读取 GTF
# ============================================================

print("读取GTF文件...")
t0 = time.time()

gtf_by_chrom = defaultdict(list)

with open(GTF_FILE, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue

        line = line.rstrip("\n")
        if not line:
            continue

        # 原始 GTF 正常情况下是 tab 分隔；如果不是 tab，就用空白分隔兜底
        if "\t" in line:
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            chrom_raw = parts[0]
            feature = parts[2]
            start = parts[3]
            end = parts[4]
            attr = parts[8]
        else:
            parts = line.split()
            if len(parts) < 9:
                continue
            chrom_raw = parts[0]
            feature = parts[2]
            start = parts[3]
            end = parts[4]
            attr = " ".join(parts[8:])

        if feature not in FEATURE_PRIORITY:
            continue

        chrom = normalize_gtf_chrom(chrom_raw)

        # 只保留主染色体，避免 NT/scaffold 影响速度
        if chrom not in {
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
            "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
            "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"
        }:
            continue

        try:
            start = int(start)
            end = int(end)
        except ValueError:
            continue

        gene_name = parse_gene_name(attr)
        rank = feature_rank(feature)

        gtf_by_chrom[chrom].append({
            "start": start,
            "end": end,
            "gene_name": gene_name,
            "feature_type": feature,
            "rank": rank,
        })

# 每条染色体按 start 排序
for chrom in gtf_by_chrom:
    gtf_by_chrom[chrom].sort(key=lambda x: x["start"])

print(f"GTF染色体数量: {len(gtf_by_chrom)}")
print("GTF染色体示例:", list(gtf_by_chrom.keys())[:10])
print(f"GTF读取耗时: {time.time() - t0:.1f}s")

# ============================================================
# 4. 扫描线注释
# ============================================================

print("开始注释位点...")

snp["in_gene"] = False
snp["gene_name"] = ""
snp["feature_type"] = ""

for chrom in sorted(snp["_CHROM_KEY"].dropna().unique()):
    snp_idx = snp["_CHROM_KEY"] == chrom
    snp_sub = snp.loc[snp_idx].copy()

    intervals = gtf_by_chrom.get(chrom, [])

    if len(intervals) == 0:
        print(f"{chrom}: 无GTF记录，跳过")
        continue

    # 按 POS 排序后扫描
    snp_sub = snp_sub.sort_values("POS")
    active = []
    j = 0
    n_intervals = len(intervals)

    hit_in_gene = {}
    hit_gene_name = {}
    hit_feature = {}

    for idx, row in snp_sub.iterrows():
        pos = row["POS"]

        if pd.isna(pos):
            hit_in_gene[idx] = False
            hit_gene_name[idx] = ""
            hit_feature[idx] = ""
            continue

        pos = int(pos)

        # 加入 start <= pos 的区间
        while j < n_intervals and intervals[j]["start"] <= pos:
            active.append(intervals[j])
            j += 1

        # 移除 end < pos 的区间
        active = [x for x in active if x["end"] >= pos]

        if not active:
            hit_in_gene[idx] = False
            hit_gene_name[idx] = ""
            hit_feature[idx] = ""
            continue

        # 当前 active 里面就是覆盖该 pos 的区间
        best = sorted(active, key=lambda x: x["rank"])[0]

        gene_names = []
        seen = set()
        for x in sorted(active, key=lambda x: x["rank"]):
            g = x["gene_name"]
            if g and g not in seen:
                gene_names.append(g)
                seen.add(g)

        hit_in_gene[idx] = True
        hit_gene_name[idx] = ";".join(gene_names)
        hit_feature[idx] = best["feature_type"]

    snp.loc[list(hit_in_gene.keys()), "in_gene"] = list(hit_in_gene.values())
    snp.loc[list(hit_gene_name.keys()), "gene_name"] = list(hit_gene_name.values())
    snp.loc[list(hit_feature.keys()), "feature_type"] = list(hit_feature.values())

    print(f"{chrom}: {snp_idx.sum():,} 个位点，注释到 {sum(hit_in_gene.values()):,} 个")

# ============================================================
# 5. 输出
# ============================================================

snp = snp.drop(columns=["_CHROM_KEY"])

front_cols = ["in_gene", "gene_name", "feature_type"]
other_cols = [c for c in snp.columns if c not in front_cols]
snp = snp[front_cols + other_cols]

# TRUE排前面，然后按 CHROM 和 POS 排序
snp = snp.sort_values(
    by=["in_gene", "CHROM", "POS"],
    ascending=[False, True, True]
).reset_index(drop=True)

snp.to_csv(OUTPUT_CSV, index=False)

n_total = len(snp)
n_in = int(snp["in_gene"].sum())

print("\n========== 注释完成 ==========")
print(f"总位点数: {n_total:,}")
print(f"注释到基因区域: {n_in:,}")
print(f"未注释到: {n_total - n_in:,}")
print(f"输出文件: {OUTPUT_CSV}")

print("\nfeature_type 分布:")
print(snp.loc[snp["in_gene"], "feature_type"].value_counts())

print(f"\n总耗时: {time.time() - t0:.1f}s")
