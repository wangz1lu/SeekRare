#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
VCF 转 CSV，只保留：
CHROM, POS, REF, ALT, 每个样本的GT

用法：
  python vcf_to_gt_csv.py input.vcf output.csv

支持：
  .vcf
  .vcf.gz
"""

import sys
import gzip
import os
import pandas as pd


if len(sys.argv) != 3:
    print(__doc__)
    sys.exit(1)

VCF_FILE = sys.argv[1]
OUT_CSV = sys.argv[2]

out_dir = os.path.dirname(OUT_CSV)
if out_dir:
    os.makedirs(out_dir, exist_ok=True)


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    else:
        return open(path, "r")


def extract_gt(sample_value, format_fields):
    """
    从样本字段中提取 GT
    例如：
      FORMAT = GT:AD:DP:GQ:PL
      sample = 0|1:10,8:18:99:...
    返回：
      0|1
    """
    if sample_value in [".", "./.", ".|."]:
        return sample_value

    values = sample_value.split(":")

    if "GT" not in format_fields:
        return ""

    gt_idx = format_fields.index("GT")

    if gt_idx >= len(values):
        return ""

    return values[gt_idx]


rows = []
sample_names = []

with open_vcf(VCF_FILE) as f:
    for line in f:
        line = line.rstrip("\n")

        # 跳过注释行
        if line.startswith("##"):
            continue

        # 读取表头
        if line.startswith("#CHROM"):
            header = line.lstrip("#").split("\t")
            sample_names = header[9:]
            continue

        if not line:
            continue

        parts = line.split("\t")

        if len(parts) < 9:
            continue

        chrom = parts[0]
        pos = parts[1]
        ref = parts[3]
        alt = parts[4]
        fmt = parts[8]

        format_fields = fmt.split(":") if fmt else []

        row = {
            "CHROM": chrom,
            "POS": pos,
            "REF": ref,
            "ALT": alt,
        }

        sample_values = parts[9:]

        for sample_name, sample_value in zip(sample_names, sample_values):
            row[sample_name] = extract_gt(sample_value, format_fields)

        rows.append(row)


df = pd.DataFrame(rows)

df["CHROM"] = df["CHROM"].astype(str)
df["POS"] = pd.to_numeric(df["POS"], errors="coerce")

df = df.sort_values(["CHROM", "POS"]).reset_index(drop=True)

df.to_csv(OUT_CSV, index=False)

print(f"VCF输入文件: {VCF_FILE}")
print(f"输出CSV文件: {OUT_CSV}")
print(f"位点数量: {len(df):,}")
print(f"样本数量: {len(sample_names):,}")
print("输出列:")
print(",".join(df.columns))
