#!/usr/bin/env bash
#
# bcftools_preprocess.sh - VCF Preprocessing Pipeline
#
# Steps:
#   1. Normalize VCFs (bcftools norm -m -both)
#   2. Merge trio (bcftools merge)
#   3. Quality filter (QUAL, DP, GQ)
#   4. Split by inheritance pattern (de novo / recessive / father_tree)
#   5. Exclude common dbSNP variants
#
# Usage:
#   bash bcftools_preprocess.sh /path/to/outdir father.vcf.gz mother.vcf.gz child.vcf.gz
#
# Requirements: bcftools, GRCh38 reference fa
#

set -euo pipefail

OUTDIR="${1:?Usage: bcftools_preprocess.sh <outdir> <father.vcf.gz> <mother.vcf.gz> <child.vcf.gz>}"
FATHER_VCF="${2:?}"
MOTHER_VCF="${3:?}"
CHILD_VCF="${4:?}"

REF="${REF:-/path/to/GRCh38_no_alt_analysis_set.fa}"   # Set REF env var
DBSNP="${DBSNP:-/mnt/dbSNP/00-common_all_chr.vcf.gz}"  # Set DBSNP env var

mkdir -p "$OUTDIR"

NORM_DIR="$OUTDIR/norm"
MERGE_DIR="$OUTDIR/merge"
FILTER_DIR="$OUTDIR/filtered
mkdir -p "$NORM_DIR" "$MERGE_DIR" "$FILTER_DIR"

echo "=== Step 1: Normalize VCFs ==="
for sample_vcf in "$FATHER_VCF" "$MOTHER_VCF" "$CHILD_VCF"; do
    name=$(basename "$sample_vcf" .vcf.gz)
    bcftools norm -m -both -f "$REF" "$sample_vcf" -Oz -o "$NORM_DIR/${name}.norm.vcf.gz"
    bcftools index "$NORM_DIR/${name}.norm.vcf.gz"
    echo "  Normalized: $name"
done

FATHER_NORM="$NORM_DIR/$(basename "$FATHER_VCF" .vcf.gz).norm.vcf.gz"
MOTHER_NORM="$NORM_DIR/$(basename "$MOTHER_VCF" .vcf.gz).norm.vcf.gz"
CHILD_NORM="$NORM_DIR/$(basename "$CHILD_VCF" .vcf.gz).norm.vcf.gz"

echo "=== Step 2: Merge trio ==="
bcftools merge "$FATHER_NORM" "$MOTHER_NORM" "$CHILD_NORM" \
    -Oz -o "$MERGE_DIR/trio.vcf.gz"
bcftools index "$MERGE_DIR/trio.vcf.gz"
echo "  Merged trio: $MERGE_DIR/trio.vcf.gz"

echo "=== Step 3: Quality filter ==="
# WGS: DP > 10-15, WES: DP > 8-10
bcftools filter -i 'QUAL>30 && FMT/DP>10 && FMT/GQ>20' \
    "$MERGE_DIR/trio.vcf.gz" -Oz -o "$FILTER_DIR/filtered.vcf.gz"
bcftools index "$FILTER_DIR/filtered.vcf.gz"
echo "  Quality filtered: $FILTER_DIR/filtered.vcf.gz"

echo "=== Step 4: Split by inheritance ==="

# De novo: parents ref/ref, child het or hom_alt
bcftools view -i ' \
    (GT[0]="0/0" || GT[0]="0|0") && \
    (GT[1]="0/0" || GT[1]="0|0") && \
    (GT[2]~"0[\/|]1" || GT[2]~"1[\/|]0" || GT[2]~"1[\/|]1") \
' "$FILTER_DIR/filtered.vcf.gz" -Oz -o "$FILTER_DIR/denovo.vcf.gz"

# AR: both parents het, child hom_alt
bcftools view -i ' \
    (GT[0]~"0[\/|]1" || GT[0]~"1[\/|]0") && \
    (GT[1]~"0[\/|]1" || GT[1]~"1[\/|]0") && \
    (GT[2]="1/1" || GT[2]="1|1") \
' "$FILTER_DIR/filtered.vcf.gz" -Oz -o "$FILTER_DIR/recessive.vcf.gz"

# Parent-only het (for compound het analysis)
bcftools view -i ' \
    (GT[0]~"0[\/|]1" || GT[0]~"1[\/|]0") && \
    (GT[1]="0/0" || GT[1]="0|0") && \
    (GT[2]="0/0" || GT[2]="0|0") \
' "$FILTER_DIR/filtered.vcf.gz" -Oz -o "$FILTER_DIR/father_tree.vcf.gz"

echo "  De novo: $FILTER_DIR/denovo.vcf.gz"
echo "  Recessive: $FILTER_DIR/recessive.vcf.gz"
echo "  Father tree: $FILTER_DIR/father_tree.vcf.gz"

if [[ -n "$DBSNP" && -f "$DBSNP" ]]; then
    echo "=== Step 5: Exclude common dbSNP ==="
    for vcf in "$FILTER_DIR"/denovo.vcf.gz "$FILTER_DIR"/recessive.vcf.gz "$FILTER_DIR"/father_tree.vcf.gz; do
        name=$(basename "$vcf" .vcf.gz)
        bcftools view -T "^$DBSNP" "$vcf" -Oz -o "$FILTER_DIR/${name}.nocommon.vcf.gz"
        echo "  Excluded common: $FILTER_DIR/${name}.nocommon.vcf.gz"
    done
else
    echo "=== Step 5: dbSNP not provided, skipping ==="
fi

echo ""
echo "=== Preprocessing Complete ==="
echo "Output directory: $OUTDIR"
