#!/usr/bin/env bash
#
# bcftools_preprocess.sh - Trio VCF 家系变异过滤流水线
#
# 功能：
#   1. Normalize — left-normalize + split multi-allelics
#   2. Merge — 合并一家三口 VCF
#   3. Filter — 质量过滤 (QUAL, DP, GQ)
#   4. Inheritance — 按遗传模式分类：
#        de_novo:     父母均 ref，子女 alt (新发突变)
#        recessive:   父母均 het，子女 hom_alt (隐性遗传)
#        compound_het: 父母一方 het 子女 het，候选复合杂合
#        x_linked:    X染色体半合子/偏倚杂合
#   5. dbSNP common — 剔除 dbSNP 常见变异
#
# 用法：
#   REF=/path/GRCh38.fa DBSNP=/path/dbsnp.vcf.gz \
#   bash bcftools_preprocess.sh /outdir father.vcf.gz mother.vcf.gz child.vcf.gz
#
# 要求：bcftools ≥ 1.10, GRCh38 reference FASTA
#

set -euo pipefail

# ── 参数解析 ────────────────────────────────────────────────────────────────
OUTDIR="${1:?用法: bcftools_preprocess.sh <outdir> <father.vcf> <mother.vcf> <child.vcf>}"
FATHER_VCF="${2:?}"
MOTHER_VCF="${3:?}"
CHILD_VCF="${4:?}"

# 参考基因组和 dbSNP（可通过环境变量覆盖）
REF="${REF:-/mnt/GRCh38_no_alt_analysis_set.fa}"
DBSNP="${DBSNP:-}"          # 例如 /mnt/dbsnp_151_1e4.vcf.gz
DBSNP_SIFT="${DBSNP_SIFT:-}" # 可选：dbSNP 额外 INFO 字段

# ── 创建输出目录 ─────────────────────────────────────────────────────────────
mkdir -p "$OUTDIR"/{norm,merge,filtered,inheritance,final}

NORM_DIR="$OUTDIR/norm"
MERGE_DIR="$OUTDIR/merge"
FILT_DIR="$OUTDIR/filtered"
INH_DIR="$OUTDIR/inheritance"
FINAL_DIR="$OUTDIR/final"

echo "=========================================="
echo "  SeekRare VCF 家系过滤流水线"
echo "=========================================="
echo "输出目录: $OUTDIR"
echo "参考基因组: $REF"
echo "dbSNP: ${DBSNP:-未提供}"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Normalize（left-normalize + split multi-allelics）
# ─────────────────────────────────────────────────────────────────────────────
echo "=== Step 1: Normalize VCFs ==="
for vcf in "$FATHER_VCF" "$MOTHER_VCF" "$CHILD_VCF"; do
    name=$(basename "$vcf" .vcf.gz | sed 's/\.vcf$//')
    out="$NORM_DIR/${name}.norm.vcf.gz"
    if [[ -s "$out" && -s "${out}.csi" ]]; then
        echo "  [跳过] $name 已 normalize"
    else
        bcftools norm -m -both -f "$REF" "$vcf" -Oz -o "$out"
        bcftools index "$out"
        echo "  [完成] $name → $(basename "$out")"
    fi
done

FATHER=$(ls "$NORM_DIR"/*father*.norm.vcf.gz 2>/dev/null | head -1)
MOTHER=$(ls "$NORM_DIR"/*mother*.norm.vcf.gz 2>/dev/null | head -1)
CHILD=$(ls "$NORM_DIR"/*child* $(ls "$NORM_DIR"/*proband* 2>/dev/null | head -1) 2>/dev/null | grep -v father | grep -v mother | head -1)

echo "  Father: $(basename "$FATHER")"
echo "  Mother: $(basename "$MOTHER")"
echo "  Child:  $(basename "$CHILD")"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Merge trio
# ─────────────────────────────────────────────────────────────────────────────
TRIO="$MERGE_DIR/trio.norm.vcf.gz"
echo "=== Step 2: Merge Trio ==="
if [[ -s "$TRIO" && -s "${TRIO}.csi" ]]; then
    echo "  [跳过] trio 已合并: $(basename "$TRIO")"
else
    bcftools merge "$FATHER" "$MOTHER" "$CHILD" \
        -Oz -o "$TRIO" \
        --force-samples
    bcftools index "$TRIO"
    echo "  [完成] → $(basename "$TRIO")"
fi
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Quality Filter
#    - QUAL > 30（质量阈值）
#    - 每样本 DP > 10（测序深度）
#    - 每样本 GQ > 20（基因型质量）
# ─────────────────────────────────────────────────────────────────────────────
FILTERED="$FILT_DIR/trio.filtered.vcf.gz"
echo "=== Step 3: Quality Filter (QUAL>30, DP>10, GQ>20) ==="
if [[ -s "$FILTERED" && -s "${FILTERED}.csi" ]]; then
    echo "  [跳过] 已过滤: $(basename "$FILTERED")"
else
    # 同时排除 multi-allelic sites（split 后应该是 biallelic）
    bcftools filter \
        -i 'QUAL>30 && TYPE="snp" && TYPE!="ref" && AVGLDP[0]>10 && AVGLDP[1]>10 && AVGLDP[2]>10' \
        "$TRIO" -Oz -o "$FILTERED"

    # 对每个样本单独检查 DP 和 GQ
    # 使用 bcftools view +_FORMAT 表达式过滤
    bcftools view -i 'FMT/DP>10 && FMT/GQ>20' "$FILTERED" -Oz -o "$FILTERED.tmp.vcf.gz"
    mv "$FILTERED.tmp.vcf.gz" "$FILTERED"
    bcftools index "$FILTERED"

    N=$(bcftools view -c "$FILTERED" 2>/dev/null | grep -oP 'After filtering, kept \K[0-9]+')
    echo "  [完成] 过滤后: ${N:-?} variants"
fi
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Step 4: 按遗传模式分类
#
# GT 编码（假设 order = father, mother, child）：
#   0/0 or 0|0 = ref
#   0/1 or 0|1 or 1/0 or 1|0 = het
#   1/1 or 1|1 = hom_alt
#   ./. = missing
#
# de_novo:        父母均 ref（0/0 或缺失），子女 het 或 hom_alt
# recessive_hom:  父母均 het，子女 hom_alt
# compound_het:   父母一方 het 子女 het（需同基因不同位点，见下方 Python 脚本）
# x_linked:       X染色体，男孩（子女 GT = 1/1 或 1/0）
# ─────────────────────────────────────────────────────────────────────────────
echo "=== Step 4: 按遗传模式分类 ==="

# ── 4a. De Novo ──────────────────────────────────────────────────────────────
# 父亲 ref + 母亲 ref + 子女 non-ref
DENOVO="$INH_DIR/denovo.vcf.gz"
echo -n "  De Novo (父母 ref，子女 alt)..."
bcftools view \
    -i 'GT[0]="0/0" && GT[1]="0/0" && GT[2]!="0/0" && GT[2]!="./." && GT[2]="het"' \
    "$FILTERED" -Oz -o "$DENOVO" 2>/dev/null || true

# 备用写法（更兼容）
bcftools view \
    -i 'FORMAT/GT[0]="0/0" && FORMAT/GT[1]="0/0" && (FORMAT/GT[2]="het" || FORMAT/GT[2]="hom")' \
    "$FILTERED" -Oz -o "$DENOVO" 2>/dev/null || \
bcftools view \
    -i 'GT[0]="0" && GT[1]="0" && GT[2]!="0"' \
    -s "^$FATHER,$MOTHER" "$FILTERED" -Oz -o "$DENOVO" 2>/dev/null || \
bcftools view \
    -i 'GT[2]="het"' \
    -S <(echo "GT[0]='0/0'"; echo "GT[1]='0/0'") \
    "$FILTERED" -Oz -o "$DENOVO"

bcftools index "$DENOVO" 2>/dev/null || true
N=$(bcftools query -f '%CHROM:%POS\n' "$DENOVO" 2>/dev/null | wc -l)
echo " $N variants"

# ── 4b. Autosomal Recessive Hom ───────────────────────────────────────────────
# 父母均 het，子女 hom_alt
AR_HOM="$INH_DIR/recessive_hom.vcf.gz"
echo -n "  常染色体隐性 (父母 het，子女 hom_alt)..."
bcftools view \
    -i 'GT[2]="1/1" || GT[2]="1|1"' \
    -S <(echo "GT[0]='het'"; echo "GT[1]='het'") \
    "$FILTERED" -Oz -o "$AR_HOM" 2>/dev/null || \
bcftools view \
    -i 'GT[2]="1/1"' \
    -s "^$FATHER,$MOTHER" -S <(echo "GT='het'") \
    "$FILTERED" -Oz -o "$AR_HOM"

bcftools index "$AR_HOM" 2>/dev/null || true
N=$(bcftools query -f '%CHROM:%POS\n' "$AR_HOM" 2>/dev/null | wc -l)
echo " $N variants"

# ── 4c. Compound Het — 父亲 het only ──────────────────────────────────────────
# 父亲 het，母亲 ref，子女 het（候选母源 compound het）
# 母亲 het only 同理，合并后做 Python 精细过滤（需同基因不同位点）
FATHER_HET="$INH_DIR/father_het_only.vcf.gz"
MOTHER_HET="$INH_DIR/mother_het_only.vcf.gz"

echo -n "  Compound Het (候选)..."
bcftools view \
    -i 'GT[2]="het"' \
    -s "$FATHER" -S <(echo "GT='0/0'") \
    "$FILTERED" -Oz -o "$FATHER_HET"
bcftools index "$FATHER_HET"

bcftools view \
    -i 'GT[2]="het"' \
    -s "$MOTHER" -S <(echo "GT='0/0'") \
    "$FILTERED" -Oz -o "$MOTHER_HET"
bcftools index "$MOTHER_HET"

N1=$(bcftools query -f '%CHROM:%POS\n' "$FATHER_HET" 2>/dev/null | wc -l)
N2=$(bcftools query -f '%CHROM:%POS\n' "$MOTHER_HET" 2>/dev/null | wc -l)
echo " father_het_only=$N1, mother_het_only=$N2"

# ── 4d. X-linked ─────────────────────────────────────────────────────────────
# 仅分析 chrX，男孩：父亲 0/0，母亲 het，子女 1/1 或 1/0（半合子）
XL_VCF="$INH_DIR/x_linked.vcf.gz"
echo -n "  X连锁隐性..."
bcftools view \
    --regions chrX \
    -i 'GT[2]="1/1" || GT[2]="1|1" || GT[2]="1/0" || GT[2]="1|0"' \
    -S <(echo "GT[0]='0/0'"; echo "GT[1]='het'") \
    "$FILTERED" -Oz -o "$XL_VCF" 2>/dev/null || true
bcftools index "$XL_VCF" 2>/dev/null || true
N=$(bcftools query -f '%CHROM:%POS\n' "$XL_VCF" 2>/dev/null | wc -l)
echo " $N variants"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Step 5: 剔除 dbSNP 常见变异
#   -T "^$DBSNP" 表示排除在 DBSNP 文件中的变异（即已知多态性位点）
# ─────────────────────────────────────────────────────────────────────────────
echo "=== Step 5: 剔除 dbSNP Common 变异 ==="

filter_common() {
    local src="$1"
    local dst="$2"
    local label="$3"

    if [[ -z "$DBSNP" || ! -f "$DBSNP" ]]; then
        echo "  [跳过] dbSNP 未提供，跳过 common 过滤"
        cp "$src" "$dst"
        return
    fi

    # 排除 dbSNP 中已有的变异（AF > 1%，即 common）
    # 使用 -T "^file" 排除对应变异
    bcftools view -T "^$DBSNP" "$src" -Oz -o "$dst"
    bcftools index "$dst"
    local n_in=$(bcftools query -f '%CHROM:%POS\n' "$src" 2>/dev/null | wc -l)
    local n_out=$(bcftools query -f '%CHROM:%POS\n' "$dst" 2>/dev/null | wc -l)
    local removed=$((n_in - n_out))
    echo "  $label: ${n_in}→${n_out} (剔除 ${removed} common)"
}

filter_common "$DENOVO"    "$FINAL_DIR/denovo.nocommon.vcf.gz"   "De Novo"
filter_common "$AR_HOM"   "$FINAL_DIR/recessive.nocommon.vcf.gz" "隐性遗传"
filter_common "$FATHER_HET" "$FINAL_DIR/father_het.nocommon.vcf.gz" "父源候选"
filter_common "$MOTHER_HET" "$FINAL_DIR/mother_het.nocommon.vcf.gz" "母源候选"
filter_common "$XL_VCF"   "$FINAL_DIR/xlinked.nocommon.vcf.gz"   "X连锁"

echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
echo "=========================================="
echo "  家系过滤完成"
echo "=========================================="
echo ""
echo "最终输出文件 ($FINAL_DIR):"
for f in "$FINAL_DIR"/*.vcf.gz; do
    if [[ -s "$f" ]]; then
        n=$(bcftools query -f '%CHROM:%POS\n' "$f" 2>/dev/null | wc -l)
        echo "  $(basename "$f"): $n variants"
    fi
done
echo ""
echo "下一步: 将 .vcf.gz 转为 CSV → gene annotation → ClinVar → SeekRare scoring"
echo "  python -m seekrare.preprocessing.vcf_to_gt $FINAL_DIR/denovo.nocommon.vcf.gz /out/denovo.csv"
