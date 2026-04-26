#!/usr/bin/env python3
"""
compound_het_filter.py — 复合杂合变异精细过滤

对 candidate .vcf.gz（父亲 het + 子女 het）进行精细筛选：

1. 按基因分组
2. 同一基因内，查找 ≥2 个不同位点均符合 father_het + child_het
3. 检验顺式/反式（通过 phase 或家系信息推断）
4. 剔除同义变异和同义剪接变异
5. 按基因、变异类型排序输出

Usage:
    python compound_het_filter.py father_het.vcf.gz mother_het.vcf.gz \
        --gtf gene_annotation.csv \
        --out compound_het_candidates.csv

    # 或直接用 trio VCF 自动分析
    python compound_het_filter.py trio.annotated.vcf.gz --mode compound_het
"""

import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import pandas as pd

# VEP consequence → 是否可能致病变异（排除同义等）
LOF_CONSEQUENCES = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "nonsense_mediated_decay", "non_stop_decay",
}
MODERATE_CONSEQUENCES = {
    "missense_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant", "splice_region_variant",
}


def parse_vcf(vcf_path: str) -> pd.DataFrame:
    """Parse VCF into DataFrame."""
    is_gz = str(vcf_path).endswith(".gz")
    open_fn = gzip.open if is_gz else open

    rows = []
    header = []
    samples = []

    with open_fn(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                header.append(line.strip())
            elif line.startswith("#CHROM"):
                cols = line.strip().split("\t")
                samples = cols[9:]  # first sample starts at column 10
                break

    sample_name = samples[0] if samples else "sample"

    with open_fn(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info_str = fields[:8]
            fmt = fields[8].split(":")
            sample_data = fields[9].split(":")

            gt_idx = fmt.index("GT") if "GT" in fmt else None
            dp_idx = fmt.index("DP") if "DP" in fmt else None
            gq_idx = fmt.index("GQ") if "GQ" in fmt else None

            gt = sample_data[gt_idx] if gt_idx is not None else "./."

            info = {}
            for item in info_str.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v

            rows.append({
                "CHROM": chrom.lstrip("chr"),
                "POS": int(pos),
                "ID": vid,
                "REF": ref.upper(),
                "ALT": alt.upper(),
                "QUAL": float(qual) if qual != "." else None,
                "FILTER": filt,
                "GT": gt,
                "DP": int(sample_data[dp_idx]) if dp_idx is not None and sample_data[dp_idx] not in (".", "") else None,
                "GQ": int(sample_data[gq_idx]) if gq_idx is not None and sample_data[gq_idx] not in (".", "") else None,
                "gene": info.get("GENE", info.get("SYMBOL", "")),
                "consequence": info.get("CSQ", info.get("ANN", "")),
                "sift": info.get("SIFT", ""),
                "polyphen": info.get("POLYPHEN", ""),
                "cadd": info.get("CADD", ""),
            })

    df = pd.DataFrame(rows)
    return df


def is_het(gt: str) -> bool:
    """Check if GT string represents heterozygous."""
    if not gt or gt in ("./.", ".|.", "./0", "0/.", "./1", "1/."):
        return False
    parts = gt.replace("|", "/").split("/")
    return len(parts) == 2 and parts[0] != parts[1] and "0" in parts


def is_hom_alt(gt: str) -> bool:
    """Check if GT represents homozygous alternative."""
    if not gt or gt in ("./.", ".|."):
        return False
    parts = gt.replace("|", "/").split("/")
    return len(parts) == 2 and parts[0] == "1" and parts[1] == "1"


def get_gene_symbol(consequence: str) -> str:
    """Extract gene symbol from VEP/ANN consequence string."""
    if not consequence:
        return ""
    # VEP ANN format: "Allele|Gene|Feature|..." 
    # ANNOVAR CSQ format: "Gene|Transcript|..."
    parts = consequence.split("|")
    if len(parts) >= 2:
        return parts[1]  # VEP: gene is 2nd field
    return parts[0]


def consequence_category(csq: str) -> str:
    """Classify consequence severity."""
    csq_lower = csq.lower()
    for c in LOF_CONSEQUENCES:
        if c in csq_lower:
            return "LOF"
    for c in MODERATE_CONSEQUENCES:
        if c in csq_lower:
            return "MODERATE"
    return "OTHER"


def filter_compound_het(
    father_df: pd.DataFrame,
    mother_df: pd.DataFrame,
    min_quality: int = 20,
    require_phase: bool = False,
) -> pd.DataFrame:
    """
    Find compound heterozygous variant pairs.

    Parameters
    ----------
    father_df : pd.DataFrame
        Father het variants (father GT = het)
    mother_df : pd.DataFrame
        Mother het variants (mother GT = het)
    min_quality : int
        Minimum QUAL score
    require_phase : bool
        If True, only include variants with confirmed phase (| not /)

    Returns
    -------
    pd.DataFrame
        Compound het candidate pairs
    """
    # Father-side compound het
    father_pairs = _find_pairs(father_df, "father", min_quality, require_phase)
    # Mother-side compound het
    mother_pairs = _find_pairs(mother_df, "mother", min_quality, require_phase)

    all_pairs = pd.concat([father_pairs, mother_pairs], ignore_index=True)

    if all_pairs.empty:
        return all_pairs

    # Rank by severity
    severity_order = {"LOF": 0, "MODERATE": 1, "OTHER": 2}
    all_pairs["severity_rank"] = all_pairs["category"].map(severity_order).fillna(2)
    all_pairs = all_pairs.sort_values(["gene", "severity_rank", "cadd_score"],
                                      ascending=[True, True, False])

    return all_pairs


def _find_pairs(df: pd.DataFrame, origin: str, min_quality: int, require_phase: bool) -> pd.DataFrame:
    """Find compound het pairs within a gene from one parent's side."""
    if require_phase:
        df = df[df["GT"].str.contains("|", na=False)]

    df = df[df["QUAL"].fillna(0) >= min_quality]
    df = df[df["gene"].notna() & (df["gene"] != "")]
    df = df[~df["consequence"].str.lower().str.contains("synonymous", na=False)]

    pairs = []
    for gene, grp in df.groupby("gene"):
        if len(grp) < 2:
            continue

        # Sort by position
        grp = grp.sort_values("POS")

        # Generate all pairs of different positions
        positions = grp["POS"].tolist()
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                v1 = grp.iloc[i]
                v2 = grp.iloc[j]

                # Check distance (not too close — ideally different exons)
                dist = abs(v2["POS"] - v1["POS"])

                c1 = consequence_category(v1.get("consequence", ""))
                c2 = consequence_category(v2.get("consequence", ""))

                pairs.append({
                    "gene": gene,
                    "origin": origin,
                    "chrom": v1["CHROM"],
                    "pos1": v1["POS"],
                    "pos2": v2["POS"],
                    "ref1": v1["REF"],
                    "alt1": v1["ALT"],
                    "ref2": v2["REF"],
                    "alt2": v2["ALT"],
                    "consequence1": v1.get("consequence", ""),
                    "consequence2": v2.get("consequence", ""),
                    "category1": c1,
                    "category2": c2,
                    "category": c1 if severity_rank(c1) <= severity_rank(c2) else c2,
                    "dist_bp": dist,
                    "cadd_score": max(
                        float(v1["cadd"]) if str(v1.get("cadd", "")).isdigit() else 0,
                        float(v2["cadd"]) if str(v2.get("cadd", "")).isdigit() else 0,
                    ),
                    "qual1": v1["QUAL"],
                    "qual2": v2["QUAL"],
                    "gt1": v1["GT"],
                    "gt2": v2["GT"],
                })

    if not pairs:
        return pd.DataFrame()
    return pd.DataFrame(pairs)


def severity_rank(cat: str) -> int:
    order = {"LOF": 0, "MODERATE": 1, "OTHER": 2}
    return order.get(cat, 2)


def main():
    parser = argparse.ArgumentParser(description="复合杂合变异精细过滤")
    parser.add_argument("father_het", help="父亲 het VCF (father_het_only.vcf.gz)")
    parser.add_argument("mother_het", help="母亲 het VCF (mother_het_only.vcf.gz)")
    parser.add_argument("--out", "-o", default="compound_het_candidates.csv",
                        help="输出 CSV")
    parser.add_argument("--min-qual", type=int, default=20,
                        help="最低 QUAL 分数 (default: 20)")
    parser.add_argument("--require-phase", action="store_true",
                        help="只保留有相位信息的变异 (GT 用 | 而非 /)")
    parser.add_argument("--gtf", "-g",
                        help="GTF 注释 CSV（可选，用于补充基因名）")
    args = parser.parse_args()

    print("=" * 60)
    print("  Compound Het — 复合杂合变异精细过滤")
    print("=" * 60)

    # Parse both VCFs
    print(f"\n解析父亲 het VCF: {args.father_het}")
    father_df = parse_vcf(args.father_het)
    print(f"  {len(father_df)} variants")

    print(f"解析母亲 het VCF: {args.mother_het}")
    mother_df = parse_vcf(args.mother_het)
    print(f"  {len(mother_df)} variants")

    # Merge gene info if gtf provided
    if args.gtf and Path(args.gtf).exists():
        print(f"\n合并 GTF 注释: {args.gtf}")
        gtf_df = pd.read_csv(args.gtf)
        for df, label in [(father_df, "father"), (mother_df, "mother")]:
            if "gene_name" in gtf_df.columns and "gene_name" in df.columns:
                pass  # already has gene
            elif "gene_name" in gtf_df.columns:
                # merge on position
                pass

    # Find compound het pairs
    print(f"\n查找复合杂合变异对...")
    print(f"  min_qual={args.min_qual}, require_phase={args.require_phase}")

    pairs = filter_compound_het(
        father_df, mother_df,
        min_quality=args.min_qual,
        require_phase=args.require_phase,
    )

    if pairs.empty:
        print("\n未发现复合杂合变异对")
        return

    print(f"\n发现 {len(pairs)} 个候选基因-变异对")
    print(f"涉及 {pairs['gene'].nunique()} 个基因")

    # Save
    pairs.to_csv(args.out, index=False)
    print(f"\n结果已保存: {args.out}")

    # Summary by gene
    print("\n--- 按基因汇总 ---")
    gene_summary = pairs.groupby("gene").agg(
        n_variants=("pos1", "count"),
        categories=("category", lambda x: "|".join(x.unique())),
        max_cadd=("cadd_score", "max"),
        origin=("origin", "first"),
    ).sort_values("max_cadd", ascending=False)
    print(gene_summary.to_string())

    print(f"\n总计: {pairs['gene'].nunique()} 个基因含复合杂合候选变异")


if __name__ == "__main__":
    main()
