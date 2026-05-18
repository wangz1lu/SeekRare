# SeekRare

**罕见病候选变异诊断系统 — 三阶段 VCF-to-Ranking**

```
VCF (proband ± father + mother)
         │
  ┌──────┴─────────────────────────────────┐
  │         Stage 1（必须）                 │
  │  家系 trio: bcftools merge +           │
  │    de_novo/recessive/xlinked 分类      │
  │  单样本: VCF → GT → GTF → ClinVar      │
  │  → simplify_clinvar_csv（精简列）      │
  └──────┬─────────────────────────────────┘
         │  3_clinvar_annotated.csv
  ┌──────┴─────────────────────────────────┐
  │         Stage 2（可选，多步）            │
  │  - stage2_eqtl_annotation: GTEx eQTL   │
  │  - stage2_splicevardb_annotation        │
  │  - stage2_omim_hpo_annotation          │
  └──────┬─────────────────────────────────┘
         │  6_omim_hpo_annotated.csv（示例）
  ┌──────┴─────────────────────────────────┐
  │         Stage 3（必须）                 │
  │  LLM 动态评分 → 排序 Top-K              │
  └──────┬─────────────────────────────────┘
         │  stage3_ranked.csv
  ┌──────┴─────────────────────────────────┐
  │         Stage 4（可选）                 │
  │  Genos 模型分析 / AlphaFold3 预测        │
  └─────────────────────────────────────────┘
```

---

## 目录

- [安装](#安装)
- [架构概览](#架构概览)
- [Stage 1 — VCF 预处理 & 基本注释](#stage-1--vcf-预处理--基本注释)
- [Stage 2 — 高级注释（可选，多步）](#stage-2--高级注释可选多步)
- [Stage 3 — LLM 动态评分排序](#stage-3--llm-动态评分排序)
- [Stage 4 — Genos / AlphaFold3（可选）](#stage-4--genos--alphafold3可选)
- [完整 Pipeline 用法](#完整-pipeline-用法)
- [输出文件说明](#输出文件说明)

---

## 安装

```bash
pip install seekrare

# 可选：AlphaFold3 支持
pip install seekrare[alphafold]

# 可选：Genos 模型支持
pip install seekrare[genomodel]

# 源码开发模式
git clone https://github.com/wangz1lu/SeekRare.git
cd SeekRare
pip install -e .
```

**Stage 1 运行所需的参考文件（用户自行准备）：**

| 文件 | 说明 |
|------|------|
| `ref_fasta` | 参考基因组 FASTA（GRCh38） |
| `gtf_file` | 基因注释 GTF 文件 |
| `clinvar_vcf` | ClinVar VCF.gz（ClinVar 2024+） |
| `dbsnp_vcf` | dbSNP VCF.gz（可选，common 过滤用） |

---

## 架构概览

### 核心原则

- **Stage 1 必须**，负责 VCF 标准化 + 基本注释
- **Stage 2 可选**，在 Stage 1 基础上追加高级注释（eQTL / SpliceVARDB / OMIM-HPO）
- **Stage 3 必须**，LLM 动态打分 + 排序
- **Stage 4 可选**，Genos 模型 / AlphaFold3 结构验证

### 模式检测（自动）

| 参数配置 | 模式 |
|----------|------|
| 只有 `vcf_proband` | 单样本模式 |
| `vcf_proband` + `vcf_father` + `vcf_mother` | 家系 trio 模式 |

---

## Stage 1 — VCF 预处理 & 基本注释

### 家系模式流程

```
三口 VCF → bcftools norm（左规范化 + 拆分 multi-allelic）
       ↓
  bcftools merge（三口合并）
       ↓
  bcftools filter（QUAL>30, DP>10, GQ>20）
       ↓
  bcftools norm --check-ref w（严格规范化）
       ↓
  去除任意样本基因型缺失的位点
       ↓
  dbSNP common 过滤（AF>0.01 的位点移除）← 若提供了 dbsnp_vcf
       ↓
  ┌──────────────────────────────────┐
  │ de_novo:     父母均 0/0，子女 alt  │
  │ recessive:   父母均 het，子女 hom  │
  │ xlinked:     父亲 het，母亲 0/0    │
  └──────────────────────────────────┘
       ↓
  各模式分别：GT CSV → GTF 注释 → ClinVar 注释
       ↓
  simplify_clinvar_csv（精简列）
       ↓
  合并所有模式 → 3_clinvar_annotated.csv
```

### 单样本模式流程

```
VCF → bcftools norm（左规范化）
       ↓
  dbSNP common 过滤（若提供 dbsnp_vcf）
       ↓
  VCF → GT CSV（提取 GT 信息）
       ↓
  GTF 注释（基因名、功能区）
       ↓
  ClinVar 注释（致病性、疾病、星级）
       ↓
  simplify_clinvar_csv（精简列）
       ↓
  3_clinvar_annotated.csv
```

### Stage 1 精简列规则

`simplify_clinvar_csv` 执行以下 9 步：

| 步骤 | 操作 | 说明 |
|------|------|------|
| 1 | 删除 `ORIGIN` 列 | 不需要 |
| 2 | `CLNREVSTAT` → `clinvarstar` | 0~4 星映射 |
| 3 | 删除 `CLNHGVS` 列 | 不需要 |
| 4 | `CLNDN` → `diseasename` | 疾病名称 |
| 5 | `CLNSIG` → `significance` | 过滤 benign/likely_benign |
| 6 | 删除 `CLNVC` 列 | 不需要 |
| 7 | `MC` 列保留 `|` 右侧 | 如 `SO:0001583` |
| 8 | `CLNDISDB` → `HPO` / `OMIM` / `Orphanet` | 拆分 ontology 标签 |
| 9 | 删除 `CLNDISDB` 和 `CLNREVSTAT` 列 | 已拆分，原始列移除 |

### Stage 1 输出列

```
CPRA, CHROM, POS, REF, ALT, in_gene, gene_name, feature_type,
inheritance_mode, clinvarstar, significance, diseasename,
HPO, OMIM, Orphanet, MC, rank
```

### Stage 1 调用方式

**方式 1：通过 SeekRarePipeline（自动检测家系/单样本）**

```python
from seekrare import SeekRarePipeline

pipeline = SeekRarePipeline(
    vcf_proband="child.vcf.gz",
    vcf_father="father.vcf.gz",     # 三者都提供 → 家系模式
    vcf_mother="mother.vcf.gz",
    ref_fasta="/path/to/hg38.fa",
    gtf_file="/path/to/genomic.gtf",
    clinvar_vcf="/path/to/clinvar.vcf.gz",
    dbsnp_vcf="/path/to/dbsnp.vcf.gz",   # 可选
    work_dir="seekrare_output",
)
pipeline.stage1_preprocess()
# 输出: seekrare_output/3_clinvar_annotated.csv
```

**方式 2：直接调用 `run_family_preprocess`（独立函数）**

```python
from seekrare.preprocess import run_family_preprocess

result = run_family_preprocess(
    work_dir="seekrare_output",
    child_vcf="child.vcf.gz",
    father_vcf="father.vcf.gz",
    mother_vcf="mother.vcf.gz",
    ref_fasta="/path/to/hg38.fa",
    gtf_file="/path/to/genomic.gtf",
    clinvar_vcf="/path/to/clinvar.vcf.gz",
    dbsnp_vcf="/path/to/dbsnp.vcf.gz",
)
# result = {"modes": {...}, "combined": DataFrame, "output_csv": str}
```

**方式 3：直接调用 Stage 1 各步骤函数（独立使用）**

```python
from seekrare.preprocess import (
    stage1_vcf_to_gt_csv,      # VCF → GT CSV
    stage1_annotate_by_gtf,     # GT CSV → GTF 注释
    stage1_merge_filter_clinvar, # ClinVar 注释
    stage1_dbsnp_filter,        # dbSNP common 过滤
    simplify_clinvar_csv,      # 精简 CSV
)
```

---

## Stage 2 — 高级注释（可选，多步）

Stage 2 由多步可选模块组成，顺序执行（均为独立函数，可单独使用）：

```
Stage 1 输出 (3_clinvar_annotated.csv)
       │
  ┌────┴────────────────────────────┐
  │  step 1: stage2_eqtl_annotation │
  │  GTEx eQTL 组织注释              │
  └────┬────────────────────────────┘
       │ (若有 tissue_dir)
  ┌────┴────────────────────────────┐
  │  step 2: stage2_splicevardb_... │
  │  SpliceVARDB 剪接注释           │
  └────┬────────────────────────────┘
       │ (若提供 SpliceVARDB TSV)
  ┌────┴────────────────────────────┐
  │  step 3: stage2_omim_hpo_...   │
  │  OMIM + HPO 二次注释           │
  └────┴────────────────────────────┘
         ↓
  6_omim_hpo_annotated.csv（示例）
```

### Step 1: GTEx eQTL 注释（stage2_eqtl_annotation）

**功能**：根据患者症状，LLM 选择相关组织，从 GTEx parquets 中匹配 eQTL 数据。

```python
from seekrare import stage2_eqtl_annotation

stage2_eqtl_annotation(
    stage1_csv="seekrare_output/3_clinvar_annotated.csv",
    tissue_dir="/path/to/GTEx_v11_eQTL_parquets",
    symptoms="地中海贫血",
    output_csv="seekrare_output/4_eqtl_annotated.csv",
    llm_model="deepseek-v4-flash",
    api_key="sk-xxxxxxxx",
    base_url="https://api.deepseek.com",
)
```

**新增列**：`eqtl_gene`, `eqtl_slope`, `eqtl_pval`, `eqtl_tissue`, `n_eqtl_tissues`

---

### Step 2: SpliceVARDB 注释（stage2_splicevardb_annotation）

**功能**：用 SpliceVARDB TSV 的 hg38 列匹配变异，注释 classification。

```python
from seekrare import stage2_splicevardb_annotation

stage2_splicevardb_annotation(
    input_csv="seekrare_output/4_eqtl_annotated.csv",
    splicevardb_tsv="/path/to/splicevardb.download.tsv",
    output_csv="seekrare_output/5_splicevardb_annotated.csv",
)
```

**新增列**：`splicevardb`（值如 `Splice-altering` / `Low-frequency` 等）

---

### Step 3: OMIM + HPO 二次注释（stage2_omim_hpo_annotation）

**功能**：
- 用 `genemap2.txt` 根据 gene_name 补充 OMIM 号（已有值不覆盖）
- 用 `phenotype.hpoa` 根据 OMIM 追加 HPO 标签（不覆盖已有）
- 用 `phenotype.hpoa` 的 disease_name 补充 diseasename（已有值不覆盖）

```python
from seekrare import stage2_omim_hpo_annotation

stage2_omim_hpo_annotation(
    input_csv="seekrare_output/5_splicevardb_annotated.csv",
    genemap2_path="/path/to/2022_05_05-genemap2.txt",
    mimtitles_path="/path/to/mimTitles.txt",
    phenotype_hpoa_path="/path/to/phenotype.hpoa",
    output_csv="seekrare_output/6_omim_hpo_annotated.csv",
)
```

**更新列**：`OMIM`（补充空值）、`HPO`（追加新标签）、`diseasename`（补充空值）

---

### 通过 Pipeline 统一调用 Stage 2

```python
pipeline = SeekRarePipeline(
    vcf_proband="child.vcf.gz",
    work_dir="seekrare_output",
    gtex_tissue_dir="/path/to/GTEx_v11_eQTL_parquets",    # eQTL 必填
    splicevardb_tsv="/path/to/splicevardb.download.tsv", # SpliceVARDB 必填
    genemap2_path="/path/to/genemap2.txt",                # OMIM-HPO 必填
    mimtitles_path="/path/to/mimTitles.txt",
    phenotype_hpoa_path="/path/to/phenotype.hpoa",
)
pipeline.stage1_preprocess()

# 逐个调用（更清晰）
pipeline.stage2_eqtl_annotation(symptoms="地中海贫血")
pipeline.stage2_splicevardb_annotation()
pipeline.stage2_omim_hpo_annotation()
```

---

## Stage 3 — LLM 动态评分排序

### 评分体系

**静态列（内置映射表）**：

| 列 | 映射规则 |
|----|---------|
| `feature_type` | CDS=1.0, exon=0.9, gene=0.7, start_codon=0.8, stop_codon=0.8, transcript=0.5 |
| `significance` | `Benign/Likely_benign` 取最坏（→ -1.0），`Pathogenic` → 1.0，`/` 分隔按最坏端；空 → -0.5 |
| `clinvarstar` | 0~5 直接映射 |
| `eqtl_tissue` | 有内容=0.5，空=0 |
| `splicevardb` | `Splice-altering`=1.0, `Low-frequency`=0.6, `Conflicting`=0.5, `Normal`=0.2 |

**动态列（LLM 根据症状给每个取值打分）**：

| 列 | 说明 |
|----|------|
| `gene_name` | LLM 根据症状决定每个基因的相关性分数 |
| `diseasename` | LLM 根据症状决定每个疾病名的相关性分数 |
| `HPO` | HP:xxxx 标签与症状的语义相关性 |
| `OMIM` | OMIM:xxxxx 与症状的相关性 |
| `Orphanet` | Orphanet:xxxxx 与症状的相关性 |
| `inheritance_mode` | LLM 根据疾病特征决定 de_novo/recessive/xlinked 谁更高分 |
| `MC` | SO:xxxx 标签的功能影响打分 |

**一致性奖励**：`gene_name` 和 `diseasename` 同时高分（>0.6）时，额外 +0.1 奖励。

### 调用方式

**方式 1：独立函数**

```python
from seekrare import stage3_score_and_rank

top = stage3_score_and_rank(
    csv_path="seekrare_output/6_omim_hpo_annotated.csv",
    symptoms="地中海贫血",
    top_k=50,
    output_csv="seekrare_output/stage3_ranked.csv",
    llm_model="deepseek-v4-flash",
    api_key="sk-xxxxxxxx",
    base_url="https://api.deepseek.com",
)
print(top[["CPRA", "gene_name", "seekrare_score", "rank"]].head(10))
```

**方式 2：通过 Pipeline**

```python
pipeline.stage3_score_and_rank(symptoms="地中海贫血")
```

**方式 3：直接使用 Stage3Scorer**

```python
from seekrare.scoring import Stage3Scorer

scorer = Stage3Scorer(
    csv_path="seekrare_output/6_omim_hpo_annotated.csv",
    symptoms="地中海贫血",
    top_k=50,
    llm_model="deepseek-v4-flash",
    api_key="sk-xxxxxxxx",
)
df = scorer.run()
scorer.save(df, "seekrare_output/stage3_ranked.csv")
```

### Stage 3 输出列

```
原有所有列 + seekrare_score（加权总分）+ rank（排名）
```

---

## Stage 4 — Genos / AlphaFold3（可选）

### Stage 4A: Genos 模型分析

**功能**：对 Stage 3 排序结果 top-K 位点做 Genos 模型预测。

```python
from seekrare import stage4_genos_analysis

stage4_genos_analysis(
    sites="top:10",                        # 取 top 10
    stage3_csv="seekrare_output/stage3_ranked.csv",
    genome_fa="/path/to/hg38.fa",
    model_path="/path/to/Genos-1.2B",
    output_dir="seekrare_output/genos_result",
)
```

**依赖**：Genos 模型权重目录（需提前部署）。

---

### Stage 4B: AlphaFold3 结构预测

```python
from seekrare import stage4_alphafold_prediction

stage4_alphafold_prediction(
    csv_path="seekrare_output/stage3_ranked.csv",
    ref_fasta="/path/to/hg38.fa",
    gtf_file="/path/to/genomic.gtf",
    top_n=5,
    output_dir="seekrare_output/alphafold_results",
)
```

---

## 完整 Pipeline 用法

```python
from seekrare import SeekRarePipeline

# ── 配置 ──────────────────────────────────────────────────
pipeline = SeekRarePipeline(
    # Stage 1
    vcf_proband="/path/to/proband.vcf.gz",
    vcf_father="/path/to/father.vcf.gz",
    vcf_mother="/path/to/mother.vcf.gz",
    ref_fasta="/path/to/hg38.fa",
    gtf_file="/path/to/genomic.gtf",
    clinvar_vcf="/path/to/clinvar.vcf.gz",
    dbsnp_vcf="/path/to/dbsnp.vcf.gz",
    work_dir="seekrare_output",
    # Stage 2（全部可选）
    gtex_tissue_dir="/path/to/GTEx_v11_eQTL_parquets",
    splicevardb_tsv="/path/to/splicevardb.download.tsv",
    genemap2_path="/path/to/genemap2.txt",
    mimtitles_path="/path/to/mimTitles.txt",
    phenotype_hpoa_path="/path/to/phenotype.hpoa",
    # Stage 3 LLM
    llm_model="deepseek-v4-flash",
    api_key="sk-xxxxxxxx",
    base_url="https://api.deepseek.com",
)

# ── Stage 1 ──────────────────────────────────────────────
pipeline.stage1_preprocess()

# ── Stage 2（多步可选）────────────────────────────────────
pipeline.stage2_eqtl_annotation(symptoms="地中海贫血")
pipeline.stage2_splicevardb_annotation()
pipeline.stage2_omim_hpo_annotation()

# ── Stage 3 ──────────────────────────────────────────────
result = pipeline.stage3_score_and_rank(symptoms="地中海贫血")
print(result[["CPRA", "gene_name", "seekrare_score", "rank"]].head(20))

# ── Stage 4（可选）────────────────────────────────────────
pipeline.stage4_genos_analysis(
    sites="top:10",
    genome_fa="/path/to/hg38.fa",
    model_path="/path/to/Genos-1.2B",
    output_dir="seekrare_output/genos_result",
)
```

---

## 输出文件说明

| 文件 | 对应 Stage | 说明 |
|------|-----------|------|
| `3_clinvar_annotated.csv` | Stage 1 结束 | 家系/单样本合并结果，含 inheritance_mode |
| `4_eqtl_annotated.csv` | Stage 2 Step 1 | GTEx eQTL 注释后 |
| `5_splicevardb_annotated.csv` | Stage 2 Step 2 | SpliceVARDB 注释后 |
| `6_omim_hpo_annotated.csv` | Stage 2 Step 3 | OMIM + HPO 二次注释后 |
| `stage3_ranked.csv` | Stage 3 | 排序完成，含 `seekrare_score` 和 `rank` |
| `genos_result/` | Stage 4A | Genos 模型预测结果目录 |
| `alphafold_results/` | Stage 4B | AlphaFold3 结构输出目录 |

---

## 依赖

**核心依赖**（`pip install seekrare` 自动安装）：
- pandas, numpy, pyarrow, pyyaml
- openai, anthropic
- tqdm, loguru, requests

**可选依赖**：
- `[alphafold]` — alphafold3
- `[genomodel]` — torch, transformers, pygromaix
- `[dev]` — pytest, black, ruff