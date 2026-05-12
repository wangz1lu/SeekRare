# SeekRare

**罕见病候选变异诊断系统**
<img width="1692" height="930" alt="seekrare" src="https://github.com/user-attachments/assets/99a4fe38-d278-46a9-bd2b-eceac98ae4b2" />

```
VCF (proband + father + mother)
         │
  ┌──────┴─────────────────────────────────┐
  │         Stage 1 (必须)                  │
  │  家系 trio: bcftools merge +           │
  │    de_novo / recessive / xlinked 分类  │
  │  单样本: VCF → GT → GTF → ClinVar     │
  └──────┬─────────────────────────────────┘
         │  3_clinvar_annotated.csv
  ┌──────┴─────────────────────────────────┐
  │         Stage 2 (可选)                  │
  │  GTEx eQTL 高级注释（LLM组织选择）        │
  └──────┬─────────────────────────────────┘
         │  4_eqtl_annotated.csv（若未配置 tissue_dir 则跳过）
  ┌──────┴─────────────────────────────────┐
  │         Stage 3 (必须)                  │
  │  LLM 双重动态评分 → 排序 Top-K           │
  └──────┬─────────────────────────────────┘
         │  stage3_ranked.csv
  ┌──────┴─────────────────────────────────┐
  │         Stage 4 (可选)                  │
  │  Genos 模型分析 / AlphaFold3 结构预测    │
  └───────────────────────────────────────┘
```

---

## 目录

- [快速开始](#快速开始)
- [安装](#安装)
- [四阶段架构](#四阶段架构)
  - [Stage 1: 家系预处理 & 基本注释](#stage-1-vcf-家系预处理--基本注释)
  - [Stage 2: GTEx eQTL 高级注释](#stage-2-eqtl-高级注释可选)
  - [Stage 3: LLM 双重动态评分排序](#stage-3-llm-双重动态评分排序)
  - [Stage 4: Genos & AlphaFold3](#stage-4-可选-高级验证)
- [文件结构](#文件结构)
- [输出文件说明](#输出文件说明)
- [配置文件](#配置文件)
- [依赖](#依赖)

---

## 快速开始

```bash
# pip 安装
pip install seekrare

# Stage 1/2/3/4 均可独立调用
```

```python
from seekrare import SeekRarePipeline, stage3_score_and_rank
from seekrare.preprocess import (
    stage1_vcf_to_gt_csv,
    stage1_annotate_by_gtf,
    stage1_merge_filter_clinvar,
    stage2_eqtl_annotation,
    stage4_genos_analysis,
    stage4_alphafold_prediction,
)

# ── Stage 1 ──────────────────────────────────────────────
pipeline = SeekRarePipeline(
    vcf_proband="child.vcf.gz",
    vcf_father="father.vcf.gz",
    vcf_mother="mother.vcf.gz",
    ref_fasta="/path/to/GRCh38.fa",
    gtf_file="/path/to/genomic.gtf",
    clinvar_vcf="/path/to/clinvar.vcf.gz",
    dbSNP_vcf="/path/to/dbsnp.vcf.gz",
    work_dir="seekrare_output",
)
pipeline.stage1_preprocess()

# ── Stage 3 ──────────────────────────────────────────────
top = stage3_score_and_rank(
    csv_path="seekrare_output/3_clinvar_annotated.csv",
    symptoms="眼部病变，视网膜色素变性",
    llm_model="deepseek-v4-flash",
    api_key="sk-xxxxxxxx",
    base_url="https://api.deepseek.com",
)
print(top[["CHROM","POS","REF","ALT","gene_name","final_score"]].head(10))

# ── Stage 4A Genos ───────────────────────────────────────
stage4_genos_analysis(
    sites="top:10",
    stage3_csv="seekrare_output/stage3_ranked.csv",
    genome_fa="/path/to/GRCh38.fa",
    model_path="/path/to/Genos-1.2B",
    output_dir="seekrare_output/genos_result",
)

# ── Stage 4B AlphaFold ──────────────────────────────────
stage4_alphafold_prediction(
    csv_path="seekrare_output/stage3_ranked.csv",
    ref_fasta="/path/to/GRCh38.fa",
    gtf_file="/path/to/genomic.gtf",
    top_n=5,
    output_dir="seekrare_output/alphafold_results",
)
```

---

## 安装

### 方式一: pip 直接安装（推荐）

```bash
pip install seekrare

# 可选：Stage 4 Genos 模型支持
pip install seekrare[genomodel]

# 可选：AlphaFold3 结构预测
pip install seekrare[alphafold]
```

### 方式二: 从源码安装（开发模式）

```bash
# 克隆 GitHub 仓库
git clone https://github.com/wangz1lu/SeekRare.git
cd SeekRare

# 创建 Python 3.10+ 虚拟环境（推荐）
python3 -m venv .venv
source .venv/bin/activate        # Linux / macOS
# .venv\Scripts\activate         # Windows

# 安装（开发模式，editable）
pip install -e .

# 可选分组
pip install -e ".[genomodel]"    # Genos 模型支持
pip install -e ".[alphafold]"    # AlphaFold3 支持
```

### 文件下载

Stage 1 运行需要以下参考文件（由用户自行准备，路径在 config 中指定）：

| 文件 | 来源 | 说明 |
|------|------|------|
| 参考基因组 FASTA | [GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) | `.fa` 或 `.fa.gz` |
| GTF 基因注释 | [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/) | 例如 `genomic.gtf` |
| ClinVar VCF | [ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) | `.vcf.gz` + `.tbi` |
| dbSNP VCF | [dbSNP](https://ftp.ncbi.nlm.nih.gov/snp/) | `.vcf.gz` + `.tbi`（可选） |

Stage 2 GTEx eQTL 数据：
```bash
# 下载 GTEx Analysis v11 eQTL parquets
# https://gtexportal.org/home/downloads/Adult
# 文件: GTEx_Analysis_v11_eQTL.tar（选 parquet 格式）
```

Stage 4 Genos 模型：
- 模型路径由用户配置（部署机器上已有的 Genos 模型目录）
- 模型路径示例: `/mnt/zzbnew/peixunban/group_01_model/chenjh356/hf_ckpts/Genos-1.2B`

---

## 四阶段架构

### Stage 1: VCF 家系预处理 & 基本注释

**自动检测家系模式：**

| 情况 | 模式 | 说明 |
|------|------|------|
| 只有 `vcf_proband` | 单样本模式 | VCF → GT → GTF → ClinVar |
| `vcf_proband` + `father_vcf` + `mother_vcf` | 家系 trio 模式 | bcftools merge + 遗传模式分类 |

**家系模式流程（纯 Python + bcftools subprocess）:**
```
每个样本: bcftools norm -m -both（拆分 multi-allelic + left-normalize）
         ↓
  bcftools merge（三口合并）
         ↓
  bcftools filter（QUAL>30, DP>10, GQ>20）
         ↓
  bcftools norm --check-ref w（严格规范化）
         ↓
  去除缺失基因型
         ↓
  ┌──────────────────────────────────────────┐
  │ de_novo:     父母均 0/0，子女为 alt       │
  │ recessive:   父母均 het，子女为 hom_alt    │
  │ xlinked:    父亲 het，母亲 0/0，子女 alt  │
  └──────────────────────────────────────────┘
         ↓
  各遗传模式 VCF → GT CSV → GTF → ClinVar → 合并 CSV
```

**单样本模式流程：**
```
VCF → GT CSV → GTF gene annotation → ClinVar 注释
```

**输出**: `seekrare_output/3_clinvar_annotated.csv`（含 `inheritance_mode` 列）

**用法**:
```python
from seekrare import SeekRarePipeline

# ── 家系 trio 模式 ──────────────────────────────────────────
p = SeekRarePipeline(
    vcf_proband="child.vcf.gz",
    vcf_father="father.vcf.gz",
    vcf_mother="mother.vcf.gz",
    ref_fasta="/path/to/GRCh38.fa",   # 必须（bcftools 需要）
    gtf_file="/path/to/genomic.gtf",
    clinvar_vcf="/path/to/clinvar.vcf.gz",
)
p.stage1_preprocess()

# ── 单样本模式 ─────────────────────────────────────────────
p = SeekRarePipeline(
    vcf_proband="child.vcf.gz",
    gtf_file="/path/to/genomic.gtf",
    clinvar_vcf="/path/to/clinvar.vcf.gz",
)
p.stage1_preprocess()
```

---

### Stage 2: eQTL 高级注释（可选）

**触发条件**: `gtex_tissue_dir` 路径已配置

**流程**:
```
Stage 1 CSV
    ↓
  LLM 分析症状 → 选择相关 3-8 个 GTEx 组织
    ↓
  读取各组织 .parquet 文件（variant_id 匹配）
    ↓
  追加列: eqtl_gene, eqtl_slope, eqtl_pval, eqtl_tissue, n_eqtl_tissues
```

**用法**:
```python
from seekrare.preprocess import stage2_eqtl_annotation

# 方式A — 手动指定组织（直接传入）
stage2_eqtl_annotation(
    stage1_csv="seekrare_output/3_clinvar_annotated.csv",
    tissue_dir="/path/to/GTEx_Analysis_v11_eQTL_parquet/",
    selected_tissues=["Brain_Cortex", "Eye_Retina", "Whole_Blood"],
)

# 方式B — LLM 自动选组织（根据症状自动筛选）
stage2_eqtl_annotation(
    stage1_csv="seekrare_output/3_clinvar_annotated.csv",
    tissue_dir="/path/to/GTEx_Analysis_v11_eQTL_parquet/",
    symptoms="眼部病变，视网膜色素变性",
    llm_provider="openai",
    llm_model="deepseek-v4-flash",
    api_key="sk-xxxxxxxx",           # 替换为你的 key
    base_url="https://api.deepseek.com",
)
```

**LLM 筛选逻辑**: LLM 从 56 个 GTEx 组织中，根据患者症状选择 3-8 个最相关的组织进行 eQTL 注释，避免全量扫描。

---

### Stage 3: LLM 双重动态评分排序

**输入**: `symptoms`（患者症状文字描述）

**流程**:
```
症状文本
    ↓
  HPO term 匹配（症状 → HPO 标准术语）
    ↓
  Stage 1/2 CSV 中提取: feature_type, CLNSIG, CLNDISDB, CLNDN, CLNDN 等
    ↓
  LLM 评判: 固定权重（CLNSIG/CLNREVSTAT/feature_type 值 → 分数）
             + 动态权重（CLNDISDB/CLNDN tags → 症状相关度评分）
    ↓
  加权求和 → max(tag_scores)  → final_score
    ↓
  排序 Top-K → stage3_ranked.csv
```

**固定权重（LLM 一次评判，多次复用）**:
| 列 | 值 | 权重 |
|----|-----|------|
| `CLNSIG` | Pathogenic | +4.0 |
| `CLNSIG` | Likely_pathogenic | +3.0 |
| `CLNSIG` | Uncertain_significance | +1.0 |
| `CLNSIG` | ... | ... |
| `CLNREVSTAT` | criteria_provided | +1.5 |
| `feature_type` | CDS | +2.0 |
| `feature_type` | intron | -1.0 |

**动态权重**: LLM 分析每个 `CLNDISDB` / `CLNDN` tag 与患者症状的相关性，输出 tag-specific 分数。

**用法**:
```python
from seekrare import stage3_score_and_rank

top = stage3_score_and_rank(
    csv_path="seekrare_output/3_clinvar_annotated.csv",  # 或 Stage 2 输出
    symptoms="眼部病变，视网膜色素变性",
    llm_model="deepseek-v4-flash",
    api_key="sk-xxxxxxxx",
    base_url="https://api.deepseek.com",
)
print(top.head(10))
```

---

### Stage 4: （可选）高级验证

#### 4A. Genos 模型分析

**功能**: 对候选位点进行 attention-based 序列验证，检测 variant 是否在正确基因组位置产生富集 peak。

**位点选择方式**:
```python
from seekrare import stage4_genos_analysis

# 方式A: 取 Stage 3 排序前 N 个
stage4_genos_analysis(
    sites="top:10",
    stage3_csv="seekrare_output/stage3_ranked.csv",
    genome_fa="/path/to/GRCh38.fa",
    model_path="/path/to/Genos-1.2B",
    output_dir="seekrare_output/genos_result",
)

# 方式B: 直接指定位点列表
stage4_genos_analysis(
    sites=[("chr1", 123456, "A", "G"), ("chr2", 789012, "CG", "T")],
    genome_fa="/path/to/GRCh38.fa",
    model_path="/path/to/Genos-1.2B",
    output_dir="seekrare_output/genos_result",
)

# 方式C: 指定 Stage 3 CSV 行号（1-based，支持范围）
stage4_genos_analysis(
    sites="rows:1,3,5-8",
    stage3_csv="seekrare_output/stage3_ranked.csv",
    genome_fa="/path/to/GRCh38.fa",
    model_path="/path/to/Genos-1.2B",
    output_dir="seekrare_output/genos_result",
)
```

**Genos 4 步 pipeline**:
```
Step 1: 1_build_variant_ref_csv.py
          提取 variant + reference 序列（±flank bp）

Step 2: 2_export_attention_by_variant_no_vcf.py
          Genos 模型计算 attention matrix（每 variant 一个子目录）

Step 3: 3_plot_variant_reference_logfc.py
          绘制 variant vs reference log2FC 散点图

Step 4: 4_check_truth_peak.py
          检验 variant 正确位置是否有富集 peak
          输出: is_reasonable_truth (True/False)
```

**输出文件**:
- `seekrare_output/genos_result/variant_ref_input.csv`
- `seekrare_output/genos_result/attention_result/{variant_id}/`
  - `metadata.csv`, `hap1/2_attention_collapsed.csv`, `timing_stats.json`
  - `logfc_plot_result/` — log2FC 图 + CSV
  - `truth_peak_validation_summary.csv`
- `seekrare_output/genos_result/truth_peak_validation_all_summary.csv`

---

#### 4B. AlphaFold3 蛋白结构预测

**功能**: 对候选位点进行蛋白结构预测，评估变异对蛋白结构的影响。

**用法**:
```python
from seekrare import stage4_alphafold_prediction

results = stage4_alphafold_prediction(
    csv_path="seekrare_output/stage3_ranked.csv",
    ref_fasta="/path/to/GRCh38.fa",
    gtf_file="/path/to/genomic.gtf",
    top_n=5,
    alphafold_mode="server",   # "server" (EBI) 或 "colabfold"
    # alphafold_url="https://your-colabfold-server/predict",
    # alphafold_api_key="your-key",
    output_dir="seekrare_output/alphafold_results",
)
```

**AlphaFold 两种模式**:
| 模式 | 说明 |
|------|------|
| `server` | EBI AlphaFold Server（免费，提交后轮询结果） |
| `colabfold` | 自部署 ColabFold（需 api_key 和 server URL） |

**流程**:
```
Stage 3 CSV (chrom, pos, ref, alt, gene_name)
    ↓
  GTF 解析 → 找基因 CDS 转录本
    ↓
  参考基因组提取 CDS 序列（验证 REF 匹配）
    ↓
  应用变异 ALT → 突变 CDS → 翻译为突变蛋白 AA 序列
    ↓
  AlphaFold3 提交（单链模式）
    ↓
  输出: mutant_proteins.fa + alphafold_summary.csv
```

**⚠️ 蛋白复合体限制**: 当前为单链预测模式。若变异涉及复合体蛋白（如 ABCA4 + CRB1 同属视网膜复合体），建议使用 AlphaFold3 multimer 模式单独处理。复合体定义表为待开发功能。

**输出文件**:
- `seekrare_output/alphafold_results/mutant_proteins.fa` — 突变蛋白序列
- `seekrare_output/alphafold_results/alphafold_summary.csv` — 汇总表

---

## 文件结构

```
seekrare/
├── src/seekrare/
│   ├── __init__.py              → SeekRarePipeline, SeekRareConfig
│   ├── pipeline.py              → SeekRarePipeline 主流程编排
│   ├── scoring/
│   │   ├── __init__.py
│   │   ├── stage3_prep.py       → CLNDISDB/CLNDN tag 解析，LLM prompt 构建
│   │   └── stage3_annotate.py   → Stage3Scorer（LLM → JSON → Python评分 → 排序）
│   ├── preprocess/
│   │   ├── __init__.py
│   │   ├── vcf_to_gt.py         → VCF → GT CSV
│   │   ├── bcftools_wrapper.py  → bcftools norm/merge/filter
│   │   ├── compound_het_filter.py → 复合杂合过滤
│   │   ├── gene_annotation.py   → GTF gene annotation
│   │   ├── clinvar_annotation.py → ClinVar 注释
│   │   ├── dbsnp_filter.py      → dbSNP common 过滤
│   │   ├── eqtl_annotation.py   → Stage 2: GTEx eQTL + LLM 组织选择
│   │   ├── stage4_genos.py      → Stage 4A: Genos 4步 pipeline
│   │   └── stage4_alphafold.py  → Stage 4B: AlphaFold3 蛋白结构预测
│   ├── annotation/
│   │   ├── __init__.py          → HPOMatcher, ClinVarLoader, GTExEQTLAnnotator,
│   │   │                           AlphaFold3Annotator, GenosAnnotationStub
│   │   ├── hpo_matcher.py       → HPO term 语义匹配
│   │   ├── clinvar_loader.py    → ClinVar VCF 解析
│   │   ├── gtex_eqtl.py         → GTEx eQTL parquet 读取
│   │   └── alphafold3.py        → AlphaFold3 API 封装
│   └── llm/
│       └── __init__.py          → LLM 调用封装
├── scripts/
│   ├── genos/
│   │   ├── 1_build_variant_ref_csv.py
│   │   ├── 2_export_attention_by_variant_no_vcf.py
│   │   ├── 3_plot_variant_reference_logfc.py
│   │   └── 4_check_truth_peak.py
│   ├── vcf_to_gt_csv.py
│   ├── annotate_vcf_csv_by_ncbi_gtf.py
│   ├── annotate_wgs_gtex.py
│   └── ...
├── tests/
│   ├── test_preprocess.py
│   ├── test_scoring.py
│   └── ...
├── README.md
└── pyproject.toml
```

---

## 输出文件说明

### Stage 1 输出（3_clinvar_annotated.csv）

| 列名 | 说明 |
|------|------|
| `CHROM` | 染色体 |
| `POS` | 基因组位置（1-based） |
| `REF` | 参考碱基 |
| `ALT` | 变异碱基 |
| `in_gene` | 是否在基因区间（0/1） |
| `gene_name` | 基因名 |
| `feature_type` | 变异区域类型（CDS / intron / ...） |
| `CLNDISDB` | ClinVar 数据库来源（tag 格式） |
| `CLNDN` | 疾病名（tag 格式） |
| `CLNREVSTAT` | ClinVar 审核状态 |
| `CLNSIG` | ClinVar 致病性评级 |
| `CLNVC` | 变异类型（SNV / insertion / deletion / ...） |
| `ORIGIN` | 变异来源（germline/somatic） |

### Stage 2 新增列

| 列名 | 说明 |
|------|------|
| `eqtl_gene` | eQTL 关联基因 |
| `eqtl_slope` | eQTL 效应方向（正/负） |
| `eqtl_pval` | eQTL 显著性 p 值 |
| `eqtl_tissue` | eQTL 来源组织 |
| `n_eqtl_tissues` | 该 variant 关联的 eQTL 组织数 |

### Stage 3 输出（stage3_ranked.csv）

| 列名 | 说明 |
|------|------|
| `CHROM`, `POS`, `REF`, `ALT` | 变异位点 |
| `gene_name` | 基因名 |
| `feature_type` | 区域类型 |
| `CLNSIG` | ClinVar 评级 |
| `clnsig_score` | 固定权重评分 |
| `dyn_score` | 动态症状相关评分 |
| `final_score` | 加权总分（= clnsig_score + dyn_score） |
| `rank` | 排名（1 = 最高分） |
| `tag_scores` | 各 CLNDISDB/CLNDN tag 的动态评分详情 |

### Stage 4A Genos 输出

| 文件 | 说明 |
|------|------|
| `variant_ref_input.csv` | Step 1 输出 |
| `attention_result/{variant_id}/metadata.csv` | variant 元信息 |
| `attention_result/{variant_id}/hap1/2_attention_collapsed.csv` | attention matrix |
| `attention_result/{variant_id}/logfc_plot_result/merge_variant_reference_log2FC.png` | log2FC 散点图 |
| `truth_peak_validation_all_summary.csv` | 全局 peak 验证汇总 |
| `is_reasonable_truth` | 该 variant 是否通过 peak 验证（True/False） |

---

## 配置文件

### SeekRareConfig

```python
from dataclasses import dataclass

@dataclass
class SeekRareConfig:
    # Stage 1 必需
    vcf_proband: str                    # 先证者 VCF
    vcf_father: str | None = None       # 父亲 VCF
    vcf_mother: str | None = None       # 母亲 VCF
    ref_fasta: str | None = None        # 参考基因组 FASTA
    gtf_file: str | None = None         # GTF 注释文件
    clinvar_vcf: str | None = None      # ClinVar VCF
    dbSNP_vcf: str | None = None        # dbSNP VCF（可选）
    work_dir: str = "seekrare_output"   # 输出目录

    # Stage 2（可选）
    gtex_tissue_dir: str | None = None  # GTEx eQTL parquet 目录

    # Stage 3
    llm_provider: str = "openai"
    llm_model: str = "deepseek-v4-flash"
    api_key: str | None = None
    base_url: str | None = None
    llm_temperature: float = 0.0
    top_k: int = 30

    # Stage 4A Genos（可选）
    genome_fa: str | None = None        # 参考基因组 FASTA
    genos_model_path: str | None = None # Genos 模型路径
```

---

## 依赖

```
Python >= 3.10
bcftools >= 1.18
pandas >= 2.0
numpy
loguru
requests
```

Stage 1 运行需要 `bcftools` 在 PATH 中。

Stage 2 运行需要 `openai` / `requests` 库。

Stage 4 Genos 运行需要部署机器有 Genos 模型（脚本已在 `scripts/genos/`）。

Stage 4 AlphaFold3 运行需要网络访问 EBI AlphaFold Server（或自部署 ColabFold）。

---

## 科学原理

### 双重动态评分（Stage 3）

传统 CLNSIG 评分是静态的，无法反映患者具体症状。SeekRare Stage 3 通过：

1. **HPO 症状标准化**: 将自由文本症状映射到 HPO 标准术语
2. **固定权重**: LLM 学习 ClinVar 各值的致病强度（一次性评判）
3. **动态权重**: LLM 结合患者具体症状，对每个 `CLNDISDB`/`CLNDN` tag 打出症状相关度分数
4. **最终排序**: `max(dynamic_scores) + sum(fixed_scores)`

### Genos Attention 验证（Stage 4A）

Genos 模型通过 attention matrix 判断 variant 序列与参考序列的差异是否在正确基因组位置产生富集。若 `downstream_peak_n >= 3` 且 `enrichment >= 2x background`，则 `is_reasonable_truth = True`。

### AlphaFold3 结构预测（Stage 4B）

通过在基因 CDS 层面应用变异，直接获取突变后蛋白序列（而非仅依赖预测的结构域），提交 AlphaFold3 进行结构预测，辅助判断变异对蛋白结构的影响。
