#!/usr/bin/env python3
"""
检验 variant 正确位置是否有 peak 富集
"""

import argparse
import os
import numpy as np
import pandas as pd


def find_variant_dirs(input_dir):
    variant_dirs = []

    for name in sorted(os.listdir(input_dir)):
        d = os.path.join(input_dir, name)
        if not os.path.isdir(d):
            continue

        metadata_path = os.path.join(d, "metadata.csv")
        logfc_path = os.path.join(
            d, "logfc_plot_result", "merge_variant_reference_log2FC.csv"
        )

        if os.path.exists(metadata_path) and os.path.exists(logfc_path):
            variant_dirs.append(d)

    return variant_dirs


def load_truth_position(metadata_path):
    metadata = pd.read_csv(metadata_path)

    if "pos" not in metadata.columns:
        raise ValueError(f"{metadata_path} 里必须有 pos 这一列")

    if metadata.empty:
        raise ValueError(f"{metadata_path} 是空文件")

    return int(metadata.iloc[0]["pos"])


def process_one_variant(
    variant_dir, window, background_exclude, top_quantile, fold_change,
):
    variant_name = os.path.basename(variant_dir)

    metadata_path = os.path.join(variant_dir, "metadata.csv")
    logfc_csv = os.path.join(
        variant_dir, "logfc_plot_result", "merge_variant_reference_log2FC.csv"
    )

    output_dir = os.path.join(variant_dir, "logfc_plot_result")
    os.makedirs(output_dir, exist_ok=True)

    truth_pos = load_truth_position(metadata_path)

    df = pd.read_csv(logfc_csv)
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["log2FC"] = pd.to_numeric(df["log2FC"], errors="coerce")
    df = df.dropna(subset=["position", "log2FC"])

    df["position"] = df["position"].astype(int)
    df["abs_log2FC"] = df["log2FC"].abs()

    if df.empty:
        raise ValueError(f"{logfc_csv} 没有有效 position/log2FC 数据")

    global_peak_thresh = df["abs_log2FC"].quantile(top_quantile)

    upstream_start = truth_pos - window
    upstream_end = truth_pos - 1

    downstream_start = truth_pos
    downstream_end = truth_pos + window

    near_start = truth_pos - background_exclude
    near_end = truth_pos + background_exclude

    bg_df = df[
        (df["position"] < near_start) | (df["position"] > near_end)
    ].copy()

    bg_mean = bg_df["abs_log2FC"].mean()
    bg_std = bg_df["abs_log2FC"].std()

    up_df = df[
        (df["position"] >= upstream_start) & (df["position"] <= upstream_end)
    ].copy()

    down_df = df[
        (df["position"] >= downstream_start) & (df["position"] <= downstream_end)
    ].copy()

    up_df["abs_log2FC"] = up_df["log2FC"].abs()
    down_df["abs_log2FC"] = down_df["log2FC"].abs()

    up_mean = up_df["abs_log2FC"].mean() if not up_df.empty else np.nan
    down_mean = down_df["abs_log2FC"].mean() if not down_df.empty else np.nan

    up_max = up_df["abs_log2FC"].max() if not up_df.empty else np.nan
    down_max = down_df["abs_log2FC"].max() if not down_df.empty else np.nan

    up_peak_n = int((up_df["abs_log2FC"] >= global_peak_thresh).sum())
    down_peak_n = int((down_df["abs_log2FC"] >= global_peak_thresh).sum())

    if not down_df.empty:
        peak_idx = down_df["abs_log2FC"].idxmax()
        downstream_peak_position = int(down_df.loc[peak_idx, "position"])
        downstream_peak_log2FC = float(down_df.loc[peak_idx, "log2FC"])
    else:
        downstream_peak_position = np.nan
        downstream_peak_log2FC = np.nan

    downstream_enrichment = (
        down_mean / bg_mean if bg_mean and bg_mean > 0 else np.nan
    )

    is_reasonable = (
        down_peak_n >= 3
        and downstream_enrichment >= fold_change
        and down_max >= global_peak_thresh
    )

    result = {
        "variant_dir": variant_name,
        "truth_pos": truth_pos,
        "window": window,
        "background_exclude": background_exclude,
        "global_abs_log2FC_peak_threshold": global_peak_thresh,

        "background_mean_abs_log2FC": bg_mean,
        "background_std_abs_log2FC": bg_std,

        "upstream_mean_abs_log2FC": up_mean,
        "upstream_max_abs_log2FC": up_max,
        "upstream_peak_points_n": up_peak_n,

        "downstream_mean_abs_log2FC": down_mean,
        "downstream_max_abs_log2FC": down_max,
        "downstream_peak_points_n": down_peak_n,
        "downstream_peak_position": downstream_peak_position,
        "downstream_peak_log2FC": downstream_peak_log2FC,
        "downstream_enrichment_vs_background": downstream_enrichment,

        "is_reasonable_truth": is_reasonable,
    }

    out_csv = os.path.join(output_dir, "truth_peak_validation_summary.csv")
    pd.DataFrame([result]).to_csv(out_csv, index=False)

    return result


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input_dir", required=True, help="attention_result3 总目录")
    parser.add_argument(
        "--output_csv",
        default=None,
        help="总汇总结果CSV；默认保存到 input_dir/truth_peak_validation_all_summary.csv"
    )
    parser.add_argument("--window", type=int, default=200)
    parser.add_argument("--background_exclude", type=int, default=500)
    parser.add_argument("--top_quantile", type=float, default=0.95)
    parser.add_argument("--fold_change", type=float, default=2.0)
    parser.add_argument("--max_variants", type=int, default=None)

    args = parser.parse_args()

    variant_dirs = find_variant_dirs(args.input_dir)

    if args.max_variants is not None:
        variant_dirs = variant_dirs[: args.max_variants]

    if not variant_dirs:
        raise ValueError(f"没有找到有效位点目录: {args.input_dir}")

    print(f"找到位点目录数量: {len(variant_dirs)}")

    results = []
    failed = 0

    for variant_dir in variant_dirs:
        try:
            res = process_one_variant(
                variant_dir=variant_dir,
                window=args.window,
                background_exclude=args.background_exclude,
                top_quantile=args.top_quantile,
                fold_change=args.fold_change,
            )
            results.append(res)
            print(f"完成: {os.path.basename(variant_dir)}")
        except Exception as e:
            failed += 1
            print(f"[ERROR] 处理失败: {variant_dir}")
            print(f"        {e}")

    result_df = pd.DataFrame(results)

    output_csv = args.output_csv or os.path.join(
        args.input_dir, "truth_peak_validation_all_summary.csv"
    )

    result_df.to_csv(output_csv, index=False)

    print("=" * 80)
    print(f"成功: {len(results)}")
    print(f"失败: {failed}")
    print(f"总汇总结果: {output_csv}")


if __name__ == "__main__":
    main()
