#!/usr/bin/env python3
"""
按 variant_id 分文件夹导出 attention matrix，不依赖 VCF
"""

import argparse
import json
import math
import os
import re
import time
from collections import defaultdict
from typing import Dict, List, Set

import pandas as pd
import torch
from tqdm import tqdm
from transformers import AutoModel, AutoTokenizer

from calc_flash_attention_v2 import calc_attentions


def safe_name(x: str) -> str:
    x = str(x)
    x = re.sub(r"[^\w\-.]+", "_", x)
    return x.strip("_")


def get_device(gpu_id: int = -1):
    if gpu_id >= 0 and torch.cuda.is_available():
        return torch.device(f"cuda:{gpu_id}")
    return torch.device("cpu")


def determine_position_columns(df: pd.DataFrame) -> Dict[int, int]:
    positions: Set[int] = set()

    for _, row in df.iterrows():
        for col in ["hap1_pos", "hap2_pos"]:
            pos_str = row.get(col, "")
            if isinstance(pos_str, str) and pos_str:
                for p in pos_str.split(";"):
                    if p:
                        positions.add(int(p))

    return {p: 1 for p in sorted(positions)}


def calculate_seq_chunks(seq_len: int, chunk_size: int, overlap: int) -> int:
    if seq_len <= chunk_size:
        return 1
    return math.ceil((seq_len - overlap) / (chunk_size - overlap))


def process_sequence_chunk(
    sample, seq, pos_list, model, tokenizer, device,
    chunk_size, overlap, chunk_idx, total_chunks,
    block_rows=8192, causal=True,
):
    start_pos = 0 if chunk_idx == 0 else chunk_idx * (chunk_size - overlap)
    end_pos = min(start_pos + chunk_size, len(seq))

    if chunk_idx == total_chunks - 1:
        end_pos = len(seq)
        start_pos = max(0, end_pos - chunk_size)

    start_idx = min(start_pos, len(pos_list) - 1)
    end_idx = min(end_pos, len(pos_list))

    if start_idx >= end_idx:
        return {}, 0.0

    seq_chunk = seq[start_pos:end_pos]

    if torch.cuda.is_available():
        torch.cuda.synchronize()

    t0 = time.time()
    attn_scores = calc_attentions(
        seq_chunk, model, tokenizer, device,
        causal=causal, block_rows=block_rows,
    )

    if torch.cuda.is_available():
        torch.cuda.synchronize()

    duration = time.time() - t0

    current_pos_list = pos_list[start_idx:end_idx]
    L_chunk = len(current_pos_list)

    if len(attn_scores) > L_chunk:
        attn_scores = attn_scores[:L_chunk]

    edge = overlap // 2

    valid_start = 0 if chunk_idx == 0 else min(edge, L_chunk)
    valid_end = L_chunk if chunk_idx == total_chunks - 1 else max(0, L_chunk - edge)

    if valid_start >= valid_end:
        return {}, duration

    result = {"sample": sample}

    sum_dict = defaultdict(float)
    cnt_dict = defaultdict(int)

    for local_i in range(valid_start, valid_end):
        pos = current_pos_list[local_i]
        sum_dict[pos] += float(attn_scores[local_i])
        cnt_dict[pos] += 1

    for pos, s in sum_dict.items():
        result[f"pos_{pos}"] = s / cnt_dict[pos]

    return result, duration


def process_sample(row, model, tokenizer, device, chunk_size, overlap, hap_name):
    sample = row["sample"]
    seq = row.get(f"{hap_name}_seq")
    pos_str = row.get(f"{hap_name}_pos")

    if not isinstance(seq, str) or not seq:
        return {"sample": sample}, set(), 0.0, 0, 0

    if not isinstance(pos_str, str) or not pos_str:
        return {"sample": sample}, set(), 0.0, 0, 0

    pos_list = [int(p) for p in pos_str.split(";") if p]

    if len(seq) != len(pos_list):
        raise ValueError(
            f"{sample} {hap_name}: seq长度({len(seq)}) != pos长度({len(pos_list)})"
        )

    total_chunks = calculate_seq_chunks(len(seq), chunk_size, overlap)

    chunk_results = []
    total_time = 0.0
    steps = 0

    for chunk_idx in range(total_chunks):
        chunk_result, duration = process_sequence_chunk(
            sample=sample, seq=seq, pos_list=pos_list,
            model=model, tokenizer=tokenizer, device=device,
            chunk_size=chunk_size, overlap=overlap,
            chunk_idx=chunk_idx, total_chunks=total_chunks,
        )

        if duration > 0:
            total_time += duration
            steps += 1

        if chunk_result:
            chunk_results.append(chunk_result)

    merged = {"sample": sample}
    all_cols = set()
    col_values = defaultdict(list)

    for res in chunk_results:
        for col, val in res.items():
            if col == "sample":
                continue
            all_cols.add(col)
            col_values[col].append(val)

    for col, values in col_values.items():
        merged[col] = sum(values) / len(values)

    return merged, all_cols, total_time, steps, 1


def collapse_haplotype_df(hap_df: pd.DataFrame, pos_occurrences: Dict[int, int]) -> pd.DataFrame:
    cols = ["sample"] + [f"pos_{p}" for p in sorted(pos_occurrences.keys())]
    return hap_df.reindex(columns=cols)


def process_one_variant(
    variant_id, sub_df, model, tokenizer, device,
    output_dir, chunk_size, overlap,
):
    variant_dir = os.path.join(output_dir, safe_name(variant_id))
    os.makedirs(variant_dir, exist_ok=True)

    pos_occurrences = determine_position_columns(sub_df)

    meta_cols = [
        "variant_id", "chrom", "pos", "CSV_REF", "REF", "ALT",
        "indel_delta", "window_start", "window_end", "sample", "sample_type",
    ]
    meta_cols = [c for c in meta_cols if c in sub_df.columns]

    sub_df[meta_cols].drop_duplicates().to_csv(
        os.path.join(variant_dir, "metadata.csv"),
        index=False,
    )

    hap1_results = []
    hap2_results = []
    all_cols = set()

    total_attn_time = 0.0
    total_steps = 0
    total_seqs = 0

    for _, row in sub_df.iterrows():
        h1_res, h1_cols, h1_time, h1_steps, h1_valid = process_sample(
            row, model, tokenizer, device, chunk_size, overlap, "hap1"
        )
        h2_res, h2_cols, h2_time, h2_steps, h2_valid = process_sample(
            row, model, tokenizer, device, chunk_size, overlap, "hap2"
        )

        hap1_results.append(h1_res)
        hap2_results.append(h2_res)

        all_cols.update(h1_cols)
        all_cols.update(h2_cols)

        total_attn_time += h1_time + h2_time
        total_steps += h1_steps + h2_steps
        total_seqs += h1_valid + h2_valid

    position_cols = ["sample"] + sorted(
        list(all_cols), key=lambda x: int(x.replace("pos_", ""))
    )

    hap1_df = pd.DataFrame(hap1_results).reindex(columns=position_cols)
    hap2_df = pd.DataFrame(hap2_results).reindex(columns=position_cols)

    hap1_out = collapse_haplotype_df(hap1_df, pos_occurrences)
    hap2_out = collapse_haplotype_df(hap2_df, pos_occurrences)

    hap1_out.to_csv(os.path.join(variant_dir, "hap1_attention_collapsed.csv"), index=False)
    hap2_out.to_csv(os.path.join(variant_dir, "hap2_attention_collapsed.csv"), index=False)

    stats = {
        "variant_id": variant_id,
        "rows": len(sub_df),
        "total_attention_extraction_time_seconds": total_attn_time,
        "total_attention_steps": total_steps,
        "total_sequences_processed": total_seqs,
        "average_chunks_per_sequence": total_steps / total_seqs if total_seqs else 0,
        "average_time_per_attention_step_seconds": total_attn_time / total_steps if total_steps else 0,
    }

    with open(os.path.join(variant_dir, "timing_stats.json"), "w") as f:
        json.dump(stats, f, indent=4)

    return variant_dir


def main():
    parser = argparse.ArgumentParser(
        description="按 variant_id 分文件夹导出 attention matrix，不依赖 VCF"
    )

    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--model_path", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--seq_chunk_size", type=int, default=4096)
    parser.add_argument("--seq_overlap", type=int, default=1000)
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--max_variants", type=int, default=None)

    args = parser.parse_args()

    if args.seq_overlap >= args.seq_chunk_size:
        raise ValueError("--seq_overlap 必须小于 --seq_chunk_size")

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"读取输入CSV: {args.input_csv}")
    df = pd.read_csv(args.input_csv)

    required_cols = [
        "variant_id", "sample", "sample_type",
        "hap1_seq", "hap1_pos", "hap2_seq", "hap2_pos",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"输入CSV缺少必要列: {missing}")

    variant_ids = list(df["variant_id"].drop_duplicates())
    if args.max_variants is not None:
        variant_ids = variant_ids[: args.max_variants]

    print(f"共检测到 variant_id 数量: {len(variant_ids)}")
    print(f"输出目录: {args.output_dir}")

    device = get_device(args.gpu)

    print(f"加载模型: {args.model_path}")
    print(f"使用设备: {device}")

    tokenizer = AutoTokenizer.from_pretrained(
        args.model_path, trust_remote_code=True,
    )
    model = AutoModel.from_pretrained(
        args.model_path,
        torch_dtype=torch.bfloat16,
        trust_remote_code=True,
    ).to(device)

    model.eval()

    workflow_start = time.time()
    finished_dirs = []

    for variant_id in tqdm(variant_ids, desc="Processing variants"):
        sub_df = df[df["variant_id"] == variant_id].copy()

        variant_dir = process_one_variant(
            variant_id=variant_id,
            sub_df=sub_df,
            model=model,
            tokenizer=tokenizer,
            device=device,
            output_dir=args.output_dir,
            chunk_size=args.seq_chunk_size,
            overlap=args.seq_overlap,
        )
        finished_dirs.append(variant_dir)

    total_time = time.time() - workflow_start

    summary = {
        "input_csv": args.input_csv,
        "output_dir": args.output_dir,
        "total_variants": len(variant_ids),
        "total_time_seconds": total_time,
        "total_time_human": f"{int(total_time // 3600)}h {int((total_time % 3600) // 60)}m {int(total_time % 60)}s",
        "finished_dirs": finished_dirs,
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=4)

    print("=" * 60)
    print(f"全部完成: {args.output_dir}")
    print(f"位点数量: {len(variant_ids)}")
    print(f"总耗时: {summary['total_time_human']}")
    print(f"总汇总文件: {os.path.join(args.output_dir, 'summary.json')}")


if __name__ == "__main__":
    main()
