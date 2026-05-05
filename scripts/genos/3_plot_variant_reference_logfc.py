#!/usr/bin/env python3
"""
绘制 variant vs reference log2FC 图
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def find_variant_dirs(input_dir):
    variant_dirs = []

    for name in sorted(os.listdir(input_dir)):
        d = os.path.join(input_dir, name)
        if not os.path.isdir(d):
            continue

        need_files = [
            "metadata.csv",
            "hap1_attention_collapsed.csv",
            "hap2_attention_collapsed.csv",
        ]

        if all(os.path.exists(os.path.join(d, f)) for f in need_files):
            variant_dirs.append(d)

    return variant_dirs


def normalize_sample_name(x):
    x = str(x)
    if x == "variant" or x.endswith("__variant"):
        return "variant"
    if x == "reference" or x.endswith("__reference"):
        return "reference"
    return x


def read_attention_csv(path):
    df = pd.read_csv(path)
    df = df.set_index(df.columns[0])
    df.index = [normalize_sample_name(x) for x in df.index]
    df.columns = df.columns.astype(str)

    if "variant" not in df.index or "reference" not in df.index:
        raise ValueError(f"{path} 必须包含 variant 和 reference 两行")

    return df


def get_pos_from_col(col):
    return int(str(col).replace("pos_", ""))


def calc_log2fc(df, eps=1e-6):
    variant = pd.to_numeric(df.loc["variant"], errors="coerce")
    reference = pd.to_numeric(df.loc["reference"], errors="coerce")

    out = pd.DataFrame({
        "position": [get_pos_from_col(c) for c in df.columns],
        "variant": variant.values,
        "reference": reference.values,
    })

    out["log2FC"] = np.log2((out["variant"] + eps) / (out["reference"] + eps))
    out = out.dropna().sort_values("position")
    return out


def plot_log2fc(logfc_df, truth_positions, out_png, title):
    fig, ax = plt.subplots(figsize=(16, 4))

    ax.scatter(
        logfc_df["position"],
        logfc_df["log2FC"],
        s=12,
        alpha=0.7,
        label="log2FC: variant / reference"
    )

    ax.axhline(0, color="gray", linestyle="--", linewidth=1)

    for i, p in enumerate(truth_positions):
        ax.axvline(
            p,
            color="red",
            linestyle="--",
            linewidth=1.2,
            label="truth variant position" if i == 0 else None
        )

    ax.set_xlabel("Genomic position")
    ax.set_ylabel("log2FC")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def process_one_variant_dir(variant_dir, eps):
    variant_name = os.path.basename(variant_dir)

    hap1_path = os.path.join(variant_dir, "hap1_attention_collapsed.csv")
    hap2_path = os.path.join(variant_dir, "hap2_attention_collapsed.csv")
    metadata_path = os.path.join(variant_dir, "metadata.csv")

    output_dir = os.path.join(variant_dir, "logfc_plot_result")
    os.makedirs(output_dir, exist_ok=True)

    metadata = pd.read_csv(metadata_path)

    if "pos" not in metadata.columns:
        raise ValueError(f"{metadata_path} 里必须有 pos 这一列")

    truth_positions = sorted(metadata["pos"].dropna().astype(int).unique().tolist())

    hap1_df = read_attention_csv(hap1_path)
    hap2_df = read_attention_csv(hap2_path)

    common_cols = hap1_df.columns.intersection(hap2_df.columns)
    common_cols = sorted(common_cols, key=get_pos_from_col)

    hap1_df = hap1_df[common_cols].apply(pd.to_numeric, errors="coerce")
    hap2_df = hap2_df[common_cols].apply(pd.to_numeric, errors="coerce")

    merged_df = (hap1_df + hap2_df) / 2

    merge_path = os.path.join(output_dir, "merge_attention_collapsed.csv")
    merged_df.to_csv(merge_path, index=True)

    logfc_df = calc_log2fc(merged_df, eps=eps)

    out_csv = os.path.join(output_dir, "merge_variant_reference_log2FC.csv")
    logfc_df.to_csv(out_csv, index=False)

    out_png = os.path.join(output_dir, "merge_variant_reference_log2FC.png")
    plot_log2fc(
        logfc_df=logfc_df,
        truth_positions=truth_positions,
        out_png=out_png,
        title=f"{variant_name}: variant vs reference log2FC"
    )

    print(f"完成: {variant_name}")
    print(f"  {out_csv}")
    print(f"  {out_png}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True, help="attention_result3 目录")
    parser.add_argument("--eps", type=float, default=1e-6)
    parser.add_argument("--max_variants", type=int, default=None)
    args = parser.parse_args()

    variant_dirs = find_variant_dirs(args.input_dir)

    if args.max_variants is not None:
        variant_dirs = variant_dirs[: args.max_variants]

    if not variant_dirs:
        raise ValueError(f"没有在 {args.input_dir} 下找到有效位点结果文件夹")

    print(f"找到位点文件夹数量: {len(variant_dirs)}")

    success = 0
    failed = 0

    for variant_dir in variant_dirs:
        try:
            process_one_variant_dir(variant_dir, eps=args.eps)
            success += 1
        except Exception as e:
            failed += 1
            print(f"[ERROR] 处理失败: {variant_dir}")
            print(f"        {e}")

    print("=" * 80)
    print(f"全部完成，成功: {success}，失败: {failed}")


if __name__ == "__main__":
    main()
