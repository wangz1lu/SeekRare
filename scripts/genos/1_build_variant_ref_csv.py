#!/usr/bin/env python3
"""
从位点list生成两行CSV：variant 与 reference
"""

import argparse
import csv
import os
from typing import Dict, List, Tuple


def read_fasta(fa: str) -> Dict[str, str]:
    genome = {}
    name = None
    chunks = []

    with open(fa, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seq = "".join(chunks).upper()
                    genome[name] = seq
                    if name.startswith("chr"):
                        genome[name.replace("chr", "", 1)] = seq
                    else:
                        genome["chr" + name] = seq

                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

        if name is not None:
            seq = "".join(chunks).upper()
            genome[name] = seq
            if name.startswith("chr"):
                genome[name.replace("chr", "", 1)] = seq
            else:
                genome["chr" + name] = seq

    return genome


def parse_sites(site_csv: str) -> List[Tuple[str, int, str, str]]:
    sites = []

    with open(site_csv, newline="") as f:
        reader = csv.reader(f)

        for row in reader:
            if not row or len(row) < 4:
                continue

            chrom, pos, ref, alt = [x.strip() for x in row[:4]]

            if chrom.lower() in {"chrom", "chr", "chromosome"}:
                continue

            sites.append((chrom, int(pos), ref.upper(), alt.upper()))

    return sites


def build_one_site(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    genome: Dict[str, str],
    flank: int,
):
    if chrom not in genome:
        raise ValueError(f"参考基因组里找不到染色体: {chrom}")

    chrom_seq = genome[chrom]
    chrom_len = len(chrom_seq)

    ref_from_genome = chrom_seq[pos - 1 : pos - 1 + len(ref)].upper()

    if ref_from_genome != ref.upper():
        print(
            f"[WARN] REF不匹配，使用参考基因组REF替代: "
            f"{chrom}:{pos} CSV_REF={ref}, genome_REF={ref_from_genome}"
        )
        ref = ref_from_genome

    start_1based = max(1, pos - flank)
    end_1based = min(chrom_len, pos + len(ref) - 1 + flank)

    ref_window = chrom_seq[start_1based - 1 : end_1based]

    left_len = pos - start_1based
    right_start = left_len + len(ref)

    var_seq = ref_window[:left_len] + alt + ref_window[right_start:]

    ref_pos = list(range(start_1based, end_1based + 1))

    delta_len = len(alt) - len(ref)

    if delta_len == 0:
        var_pos = ref_pos[:]

    elif delta_len < 0:
        keep_ref_len = len(alt)
        var_pos = (
            list(range(start_1based, pos + keep_ref_len))
            + list(range(pos + len(ref), end_1based + 1))
        )

    else:
        var_pos = (
            list(range(start_1based, pos))
            + [pos] * len(alt)
            + list(range(pos + len(ref), end_1based + 1))
        )

    if len(var_seq) != len(var_pos):
        raise RuntimeError(
            f"序列长度和坐标长度不一致: {chrom}:{pos}, "
            f"seq_len={len(var_seq)}, pos_len={len(var_pos)}"
        )

    if len(ref_window) != len(ref_pos):
        raise RuntimeError(
            f"参考序列长度和坐标长度不一致: {chrom}:{pos}, "
            f"seq_len={len(ref_window)}, pos_len={len(ref_pos)}"
        )

    return ref_window, ref_pos, var_seq, var_pos, delta_len, start_1based, end_1based, ref


def main():
    parser = argparse.ArgumentParser(
        description="从位点list生成两行CSV：variant 与 reference"
    )

    parser.add_argument(
        "--sites_csv",
        required=True,
        help="输入位点CSV: chrom,pos,REF,ALT，无表头或有表头均可",
    )
    parser.add_argument(
        "--fa",
        required=True,
        help="参考基因组FASTA",
    )
    parser.add_argument(
        "--output_csv",
        required=True,
        help="输出CSV",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=2000,
        help="上下游扩展长度，默认2000",
    )
    parser.add_argument(
        "--variant_sample",
        default="variant",
        help="变异序列sample名称",
    )
    parser.add_argument(
        "--ref_sample",
        default="reference",
        help="参考序列sample名称",
    )

    args = parser.parse_args()

    print(f"读取参考基因组: {args.fa}")
    genome = read_fasta(args.fa)

    print(f"读取位点文件: {args.sites_csv}")
    sites = parse_sites(args.sites_csv)

    if not sites:
        raise ValueError("sites_csv 里没有有效位点")

    rows = []
    mismatch_count = 0

    for chrom, pos, ref, alt in sites:
        original_ref = ref

        (
            ref_seq,
            ref_pos,
            var_seq,
            var_pos,
            delta_len,
            win_start,
            win_end,
            genome_ref,
        ) = build_one_site(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            genome=genome,
            flank=args.flank,
        )

        if genome_ref != original_ref:
            mismatch_count += 1

        variant_id = f"{chrom}_{pos}_{genome_ref}_{alt}"

        rows.append(
            {
                "variant_id": variant_id,
                "chrom": chrom,
                "pos": pos,
                "CSV_REF": original_ref,
                "REF": genome_ref,
                "ALT": alt,
                "indel_delta": delta_len,
                "window_start": win_start,
                "window_end": win_end,
                "sample": args.variant_sample,
                "sample_type": 1,
                "hap1_seq": var_seq,
                "hap1_pos": ";".join(map(str, var_pos)),
                "hap2_seq": var_seq,
                "hap2_pos": ";".join(map(str, var_pos)),
            }
        )

        rows.append(
            {
                "variant_id": variant_id,
                "chrom": chrom,
                "pos": pos,
                "CSV_REF": original_ref,
                "REF": genome_ref,
                "ALT": alt,
                "indel_delta": delta_len,
                "window_start": win_start,
                "window_end": win_end,
                "sample": args.ref_sample,
                "sample_type": 0,
                "hap1_seq": ref_seq,
                "hap1_pos": ";".join(map(str, ref_pos)),
                "hap2_seq": ref_seq,
                "hap2_pos": ";".join(map(str, ref_pos)),
            }
        )

    os.makedirs(os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True)

    fieldnames = [
        "variant_id", "chrom", "pos", "CSV_REF", "REF", "ALT",
        "indel_delta", "window_start", "window_end",
        "sample", "sample_type",
        "hap1_seq", "hap1_pos", "hap2_seq", "hap2_pos",
    ]

    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("=" * 60)
    print(f"完成: {args.output_csv}")
    print(f"输入位点数: {len(sites)}")
    print(f"输出行数: {len(rows)}")
    print(f"REF不匹配并替换数量: {mismatch_count}")


if __name__ == "__main__":
    main()
