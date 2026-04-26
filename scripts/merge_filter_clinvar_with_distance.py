#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re
import sys
import os
import bisect
from collections import defaultdict

if len(sys.argv) != 4:
    print("""
用法:
  python merge_filter_clinvar_with_distance.py <annotated_csv> <clinvar_csv> <output_csv>
""")
    sys.exit(1)

ANNOT_FILE = sys.argv[1]
CLINVAR_FILE = sys.argv[2]
OUTPUT_FILE = sys.argv[3]

out_dir = os.path.dirname(OUTPUT_FILE)
if out_dir:
    os.makedirs(out_dir, exist_ok=True)

print("读取注释结果...")
df = pd.read_csv(ANNOT_FILE, low_memory=False)

print("读取ClinVar...")
clin = pd.read_csv(CLINVAR_FILE, low_memory=False)

# ============================================================
# 1. 标准化染色体
# ============================================================

def norm_chrom(x):
    if pd.isna(x):
        return ""
    x = str(x).strip()
    if x == "":
        return ""
    if x.startswith("chr"):
        return x
    if x in ["X", "Y", "M", "MT"]:
        return "chrM" if x in ["M", "MT"] else "chr" + x
    return "chr" + x

df["CHROM_norm"] = df["CHROM"].apply(norm_chrom)
df["POS"] = pd.to_numeric(df["POS"], errors="coerce")

clin["CHROM_norm"] = clin["CHROM"].apply(norm_chrom)
clin["POS"] = pd.to_numeric(clin["POS"], errors="coerce")

# ============================================================
# 2. 提取 ClinVar 基因名
# ============================================================

def extract_gene(geneinfo):
    """
    OR4F5:79501 -> OR4F5
    A:1|B:2 -> A;B
    """
    if pd.isna(geneinfo) or str(geneinfo).strip() == "":
        return ""

    genes = []
    for item in str(geneinfo).split("|"):
        g = item.split(":")[0].strip()
        if g:
            genes.append(g)

    return ";".join(sorted(set(genes)))

clin["gene"] = clin["GENEINFO"].apply(extract_gene)

clin = clin.assign(gene=clin["gene"].str.split(";")).explode("gene")
clin["gene"] = clin["gene"].astype(str).str.strip()

# ============================================================
# 3. 保留必要字段
# ============================================================

clin_sub = clin[[
    "gene",
    "CHROM_norm",
    "POS",
    "CLNSIG",
    "MC",
    "CLNDISDB"
]].copy()

clin_sub = clin_sub[
    (clin_sub["gene"].notna()) &
    (clin_sub["gene"] != "") &
    (clin_sub["gene"] != "nan") &
    (clin_sub["POS"].notna())
].copy()

clin_sub["POS"] = clin_sub["POS"].astype(int)

# ============================================================
# 4. 构建 gene -> ClinVar 信息映射
# ============================================================

def join_unique(x):
    vals = []
    for v in x.astype(str):
        if v and v != "nan" and v != ".":
            vals.append(v)
    return ";".join(sorted(set(vals)))

clin_group = (
    clin_sub.groupby("gene")
    .agg({
        "CLNSIG": join_unique,
        "MC": join_unique,
        "CLNDISDB": join_unique,
    })
    .reset_index()
)

clin_dict = clin_group.set_index("gene").to_dict(orient="index")

# ============================================================
# 5. 构建 gene + chrom -> ClinVar POS列表
# ============================================================

clin_pos_dict = defaultdict(list)

for _, row in clin_sub.iterrows():
    gene = row["gene"]
    chrom = row["CHROM_norm"]
    pos = int(row["POS"])
    clin_pos_dict[(gene, chrom)].append(pos)

for key in clin_pos_dict:
    clin_pos_dict[key] = sorted(set(clin_pos_dict[key]))

# ============================================================
# 6. gene_name 匹配 ClinVar信息
# ============================================================

def match_clinvar(gene_name):
    if pd.isna(gene_name) or str(gene_name).strip() == "":
        return "", "", ""

    genes = [g.strip() for g in str(gene_name).split(";") if g.strip()]

    sigs = []
    mcs = []
    hps = []

    for g in genes:
        if g in clin_dict:
            sigs.append(clin_dict[g]["CLNSIG"])
            mcs.append(clin_dict[g]["MC"])
            hps.append(clin_dict[g]["CLNDISDB"])

    return (
        ";".join(sorted(set([x for x in sigs if x]))),
        ";".join(sorted(set([x for x in mcs if x]))),
        ";".join(sorted(set([x for x in hps if x]))),
    )

# ============================================================
# 7. 计算最近 ClinVar 位点距离
# ============================================================

def nearest_distance(gene_name, chrom, pos):
    if pd.isna(gene_name) or str(gene_name).strip() == "":
        return -1
    if pd.isna(pos):
        return -1

    pos = int(pos)
    genes = [g.strip() for g in str(gene_name).split(";") if g.strip()]

    best = None

    for g in genes:
        key = (g, chrom)

        if key not in clin_pos_dict:
            continue

        pos_list = clin_pos_dict[key]
        i = bisect.bisect_left(pos_list, pos)

        candidates = []

        if i < len(pos_list):
            candidates.append(abs(pos_list[i] - pos))

        if i > 0:
            candidates.append(abs(pos_list[i - 1] - pos))

        if candidates:
            d = min(candidates)
            if best is None or d < best:
                best = d

    return best if best is not None else -1

print("开始匹配ClinVar...")

df[["clinvar_sig", "clinvar_mc", "clinvar_hp"]] = df["gene_name"].apply(
    lambda x: pd.Series(match_clinvar(x))
)

print("计算当前位点到同基因最近ClinVar位点的距离...")

df["clinvar_min_distance"] = df.apply(
    lambda row: nearest_distance(
        row["gene_name"],
        row["CHROM_norm"],
        row["POS"]
    ),
    axis=1
)

# ============================================================
# 8. ClinVar 临床意义排序
# ============================================================

CLNSIG_PRIORITY = {
    "Pathogenic": 1,
    "Likely_pathogenic": 2,
    "Pathogenic/Likely_pathogenic": 3,
    "Conflicting_classifications_of_pathogenicity": 4,
    "Uncertain_significance": 5,
    "Likely_benign": 6,
    "Benign/Likely_benign": 7,
    "Benign": 8,
}

def get_clinvar_rank(sig):
    if pd.isna(sig) or str(sig).strip() == "":
        return 99

    items = re.split(r"[;,]", str(sig))
    ranks = []

    for item in items:
        item = item.strip()
        if item in CLNSIG_PRIORITY:
            ranks.append(CLNSIG_PRIORITY[item])

    if len(ranks) == 0:
        return 98

    return min(ranks)

df["clinvar_rank"] = df["clinvar_sig"].apply(get_clinvar_rank)

# ============================================================
# 9. 过滤
# ============================================================

print("开始过滤...")

in_gene_bool = df["in_gene"].astype(str).str.lower().isin(["true", "1"])

df_filtered = df[
    in_gene_bool &
    df["clinvar_sig"].notna() &
    (df["clinvar_sig"].astype(str).str.strip() != "")
].copy()

# ============================================================
# 10. 排序
#    distance = -1 放到最后
# ============================================================

df_filtered["distance_missing_rank"] = (
    df_filtered["clinvar_min_distance"] == -1
).astype(int)

sort_cols = ["clinvar_rank", "distance_missing_rank"]

if "clinvar_min_distance" in df_filtered.columns:
    sort_cols.append("clinvar_min_distance")
if "gene_name" in df_filtered.columns:
    sort_cols.append("gene_name")
if "CHROM" in df_filtered.columns:
    sort_cols.append("CHROM")
if "POS" in df_filtered.columns:
    sort_cols.append("POS")

df_filtered = df_filtered.sort_values(
    by=sort_cols,
    ascending=[True] * len(sort_cols)
).reset_index(drop=True)

# ============================================================
# 11. 删除辅助列
# ============================================================

drop_cols = ["clinvar_rank", "CHROM_norm", "distance_missing_rank"]
df_filtered = df_filtered.drop(columns=[c for c in drop_cols if c in df_filtered.columns])

# ============================================================
# 12. 把 clinvar_min_distance 放到第一列
# ============================================================

cols = list(df_filtered.columns)

if "clinvar_min_distance" in cols:
    cols.remove("clinvar_min_distance")
    cols.insert(0, "clinvar_min_distance")

df_filtered = df_filtered[cols]

# ============================================================
# 13. 输出
# ============================================================

df_filtered.to_csv(OUTPUT_FILE, index=False)

# ============================================================
# 14. 统计输出
# ============================================================

print("\n========== 完成 ==========")
print(f"原始数据量: {len(df):,}")
print(f"过滤后数据量: {len(df_filtered):,}")
print(f"输出文件: {OUTPUT_FILE}")

print("\nClinVar分类统计 Top10:")
print(df_filtered["clinvar_sig"].value_counts().head(10))

print("\nPathogenic相关位点数:")
print(df_filtered["clinvar_sig"].str.contains("Pathogenic", na=False).sum())

print("\n========== 距离分布统计 ==========")

dist = df_filtered["clinvar_min_distance"]
dist_valid = dist[dist >= 0]

n0 = (dist_valid == 0).sum()
n1 = ((dist_valid > 0) & (dist_valid <= 10)).sum()
n2 = ((dist_valid > 10) & (dist_valid <= 50)).sum()
n3 = ((dist_valid > 50) & (dist_valid <= 100)).sum()
n4 = (dist_valid > 100).sum()
n_missing = (dist == -1).sum()

print(f"距离 = 0 的位点数: {n0:,}")
print(f"0 < 距离 <= 10 的位点数: {n1:,}")
print(f"10 < 距离 <= 50 的位点数: {n2:,}")
print(f"50 < 距离 <= 100 的位点数: {n3:,}")
print(f"距离 > 100 的位点数: {n4:,}")
print(f"无匹配 distance = -1 的位点数: {n_missing:,}")
