# SeekRare
**三阶段罕见病诊断系统 — Three-Stage Rare Disease Diagnosis System**

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
[![develop](https://img.shields.io/badge/branch-develop-blue.svg)](https://github.com/wangz1lu/SeekRare/tree/develop)

---

## 🎯 核心创新

SeekRare 采用**三阶段架构**，从家系 VCF 输入到排序候选变异输出：

```
患者症状（自由文本）
        │
        ▼
┌─────────────────────────────────────────────────────────────┐
│  Stage 3: LLM Symptom Interpreter                           │
│  1. LLM 解析患者症状 → 相关 HPO terms + relevance score     │
│  2. LLM 输出动态 weight vector（每个患者不同）               │
│  3. HPOMatcher 语义匹配，扩展 HPO 覆盖                       │
└─────────────────────────────────────────────────────────────┘
        │
        ▼
Annotated CSV ──→ Dual-Dynamic Scorer ──→ Personalized Ranking
                   │
                   ├── Column weights: LLM-adjusted per symptom
                   └── HPO relevance: semantic similarity to patient
```

---

## 🏗️ 三阶段架构

### Stage 1 — VCF 家系预处理 + 基本注释（必须）

```
Trio VCFs (father + mother + child)
        │
        ▼
scripts/bcftools_preprocess.sh
  ├── bcftools norm (left-normalize, split multi-allelics)
  ├── bcftools merge (trio)
  ├── Quality filter (QUAL>30, DP>10, GQ>20)
  ├── Inheritance 分类: denovo / recessive / compound_het / xlinked
  └── Exclude common dbSNP variants
        │
        ▼
scripts/compound_het_filter.py   (Python 精细复合杂合过滤)
        │
        ▼
preprocess/vcf_to_gt.py           → CHROM, POS, REF, ALT, GT per sample
preprocess/gene_annotation.py     → gene_name, feature_type (GTF sweep-line)
preprocess/clinvar_annotation.py   → clinvar_sig, clinvar_mc, min_distance
                                     + omim_diseases, omim_inheritance, omim_mim_number
                                     + hpo_terms (基因→疾病映射)
        │
        ▼
Stage 1 输出: annotated_variants.csv
（基础列 + 必须注释列）
```

### Stage 2 — 高级注释（可选）

```
Stage 1 CSV + 资源配置
        │
        ▼
annotation/gtex_eqtl.py         → eqtl_gene, eqtl_pval, eqtl_tissue (GTEx eQTL)
annotation/alphafold3.py        → alphafold_predicted, alphafold_pdb_url (STUB)
annotation/alphafold3.GenosAnnotationStub → genos_pathogenicity, genos_acmg_criteria (STUB)
        │
        ▼
Stage 2 输出: 在 Stage 1 基础上追加高级注释列
```

### Stage 3 — LLM 分析

```
Annotated CSV + Patient Symptoms
        │
        ▼
LLMSymptomParser → {relevant_hpos: [...], weight_vector: {...}}
        │
        ▼
HPOMatcher (optional) → additional HPO semantic matching
        │
        ▼
DualDynamicScorer → seekrare_score per variant
        │
        ▼
Ranked Candidate Variants (personalized top-K)
```

---

## 📁 项目结构

```
seekrare/
├── README.md
├── pyproject.toml
├── requirements.txt
├── .gitignore
│
├── scripts/                               # 原始脚本（Stage 1）
│   ├── bcftools_preprocess.sh             # 家系 VCF 过滤 (bash)
│   ├── compound_het_filter.py              # 复合杂合精细过滤 (Python)
│   ├── vcf_to_gt_csv.py                   # VCF → GT CSV
│   ├── annotate_vcf_csv_by_ncbi_gtf.py    # GTF 基因注释
│   ├── merge_filter_clinvar_with_distance.py  # ClinVar + 距离
│   └── annotate_wgs_gtex.py               # GTEx eQTL 注释 (原始脚本)
│
└── src/seekrare/
    ├── __init__.py
    ├── pipeline.py                        # 三阶段流水线编排器
    │
    ├── preprocess/                         # Stage 1 Python 封装
    │   ├── vcf_to_gt.py                   # vcf_to_gt_csv()
    │   ├── gene_annotation.py              # annotate_by_gtf()
    │   ├── clinvar_annotation.py          # merge_filter_clinvar()
    │   └── bcftools_wrapper.py            # run_bcftools_preprocess()
    │
    ├── annotation/                         # 注释加载器
    │   ├── clinvar_loader.py              # ClinVar VCF/CSV loader (Stage 1)
    │   ├── hpo_matcher.py                 # HPO ontology 语义匹配 (Stage 3)
    │   ├── gtex_eqtl.py                   # GTEx eQTL 注释 (Stage 2)
    │   ├── alphafold3.py                  # AlphaFold3 + GenosAnnotationStub (Stage 2, STUB)
    │   └── combiner.py                    # 多源注释合并工具
    │
    ├── llm/                               # LLM 模块 (Stage 3)
    │   └── symptom_parser.py              # LLMSymptomParser
    │
    └── scoring/                           # 打分与排序 (Stage 3)
        ├── engine.py                      # DualDynamicScorer
        └── ranker.py                     # rank_variants()
```

---

## 🔧 安装（当前版本未发布 PyPI，请从源码安装）

```bash
# 克隆 GitHub 仓库
git clone https://github.com/wangz1lu/SeekRare.git
cd SeekRare

# 创建 Python 3.10+ 虚拟环境（推荐）
python3 -m venv .venv
source .venv/bin/activate   # Linux/macOS
# .venv\Scripts\activate     # Windows

# 安装（开发模式，editable）
pip install -e .

# 或只安装核心依赖
pip install -r requirements.txt
```

**requirements.txt 内容：**
```
pysam>=0.21
pandas>=2.0
numpy>=1.24
pyarrow>=14.0
pyyaml>=6.0
openai>=1.0
anthropic>=0.18
tqdm>=4.66
loguru>=0.7
biopython>=1.83
requests>=2.31
```

**可选依赖：**
```bash
pip install pyhpo>=3.0    # HPO ontology 语义匹配（需要科学上网）
```

---

## 🚀 快速开始

### 1. 准备文件

Stage 1 需要以下四个文件（请从对应地址下载）：

| 文件 | 推荐下载地址 |
|------|------------|
| GRCh38 参考基因组 | https://ftp.ensembl.org/pub/fasta/homo_sapiens/dna/ |
| NCBI genomic.gtf | https://ftp.ensembl.org/pub/gtf/homo_sapiens/ |
| ClinVar VCF | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/ |
| dbSNP VCF | https://ftp.ncbi.nlm.nih.gov/snp/ |

### 2. 基本使用

```python
from seekrare import SeekRarePipeline

# 初始化（填入你的文件路径）
pipeline = SeekRarePipeline(
    # ── Stage 1 ──────────────────────────────────────────────
    vcf_proband="/path/to/child.vcf.gz",
    vcf_father="/path/to/father.vcf.gz",   # 可选，无则跳过家系过滤
    vcf_mother="/path/to/mother.vcf.gz",
    ref_fasta="/path/to/GRCh38.fa",
    gtf_file="/path/to/genomic.gtf",
    clinvar_vcf="/path/to/clinvar.vcf.gz",
    dbSNP_vcf="/path/to/dbsnp.vcf.gz",     # 可选

    # ── Stage 3 ──────────────────────────────────────────────
    llm_provider="openai",
    llm_model="gpt-4o",
    api_key=os.getenv("OPENAI_API_KEY"),
)

# Stage 1 only（最快验证流程是否通）
df = pipeline.stage1_preprocess()
df.to_csv("variants_annotated.csv", index=False)

# 完整三阶段
result = pipeline.run(
    symptoms="智力障碍，癫痫，全身肌张力低，脑电图显示全面性棘慢波发放",
)
result.to_csv("candidates.csv", index=False)
```

### 3. Stage 2 高级注释（可选）

```python
# Stage 1 + Stage 2
pipeline.stage1_preprocess()
df = pipeline.stage2_advanced_annotation()  # 需要配置 gtex_tissue_dir 等

# Stage 1 + Stage 2 + Stage 3
result = pipeline.run(
    symptoms="...",
    skip_stage2=False,  # 启用 Stage 2（需配置 gtex_tissue_dir 等）
)
```

### 4. Stage 1 脚本独立使用

```bash
# bcftools 预处理（需要 bcftools 和参考基因组）
REF=/path/GRCh38.fa DBSNP=/path/dbsnp.vcf.gz \
bash scripts/bcftools_preprocess.sh /outdir father.vcf.gz mother.vcf.gz child.vcf.gz

# VCF → GT CSV
python -m seekrare.preprocessing.vcf_to_gt input.vcf output.csv

# GTF 基因注释
python -m seekrare.preprocessing.gene_annotation input.csv genomic.gtf output.csv

# ClinVar + OMIM + HPO 注释
python -m seekrare.preprocessing.clinvar_annotation annotated.csv clinvar.vcf.gz output.csv
```

---

## ⚙️ 配置参数

### SeekRareConfig 主要参数

| 参数 | 说明 | 必填 | 默认值 |
|------|------|------|--------|
| `vcf_proband` | 先证者 VCF | ✅ | — |
| `vcf_father` | 父亲 VCF | ❌ | None |
| `vcf_mother` | 母亲 VCF | ❌ | None |
| `ref_fasta` | GRCh38 参考基因组 FASTA | ❌* | None |
| `gtf_file` | NCBI genomic.gtf | ❌* | None |
| `clinvar_vcf` | ClinVar VCF | ❌* | None |
| `dbSNP_vcf` | dbSNP VCF | ❌ | None |
| `gtex_tissue_dir` | GTEx eQTL parquet 目录 | ❌ | None |
| `alphafold3_mode` | "server" 或 "colabfold" | ❌ | None |
| `llm_provider` | "openai" / "anthropic" / "local" | ❌ | "openai" |
| `llm_model` | 模型名 | ❌ | "gpt-4o" |
| `top_k` | 返回 top-K 候选数 | ❌ | 50 |

*无家系 VCF 时，ref_fasta/gtf_file/clinvar_vcf 可选（跳过对应步骤）

---

## 📊 输出列说明

### Stage 1 输出（基本注释）

| 列名 | 来源 | 说明 |
|------|------|------|
| CHROM, POS, REF, ALT | VCF | 变异坐标 |
| GT_proband, GT_father, GT_mother | VCF | 基因型 |
| inheritance_type | bcftools | denovo/recessive/compound_het/xlinked |
| gene_name | GTF | 基因符号 |
| feature_type | GTF | exon/CDS/UTR... |
| clinvar_sig | ClinVar | 致病性评级 |
| clinvar_mc | ClinVar | 分子 consequence |
| min_distance | ClinVar | 距最近 ClinVar 变异距离 |
| omim_diseases | ClinVar/OMIM | OMIM 疾病 |
| omim_inheritance | OMIM | 遗传模式 |
| hpo_terms | ClinVar/OMIM | 基因关联的 HPO terms |

### Stage 2 输出（追加列）

| 列名 | 来源 | 说明 |
|------|------|------|
| eqtl_gene | GTEx | eQTL 关联基因 |
| eqtl_pval | GTEx | eQTL p-value |
| eqtl_tissue | GTEx | eQTL 组织来源 |
| alphafold_predicted | AlphaFold3 | 是否已预测（STUB） |
| genos_pathogenicity | Genos | Genos 模型评分（STUB） |

### Stage 3 输出（追加列）

| 列名 | 说明 |
|------|------|
| seekrare_score | 双动态打分最终分数 |
| hpo_similarity | HPO 语义相似度 |
| rank | 最终排名 |

---

## 🔬 科学原理

### 双动态打分（Dual-Dynamic Scoring）

| 维度 | 传统方法 | SeekRare |
|------|---------|---------|
| Column weights | 固定权重 | **LLM 按患者症状动态调整** |
| HPO relevance | 二元匹配 | **语义相似度到患者症状** |
| 输出 | 通用排序 | **个性化排名** |

---

## 📄 License

MIT
