# SeekRare
**LLM-powered Rare Disease Diagnosis System with Dual-Dynamic Variant Scoring**

[![PyPI](https://img.shields.io/pypi/v/seekrare.svg)](https://pypi.org/project/seekrare/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 📦 Download & Installation

```bash
# Install from PyPI (recommended)
pip install seekrare

# Or install latest development version from GitHub
pip install git+https://github.com/wangz1lu/SeekRare.git

# With all optional dependencies
pip install seekrare[dev]

# Development install (editable)
git clone https://github.com/wangz1lu/SeekRare.git
cd SeekRare
pip install -e .
```

**Requirements:**
- Python ≥ 3.10
- pandas ≥ 2.0, numpy ≥ 1.24
- pysam ≥ 0.21 (VCF handling)
- openai ≥ 1.0 or anthropic ≥ 0.18 (LLM APIs)

---

## 🎯 Core Innovation

Traditional variant prioritization tools apply **fixed, static weights** to annotation columns. SeekRare uses an **LLM-driven dual-dynamic scoring system**:

```
Patient Symptoms (free text)
        │
        ▼
┌─────────────────────────────────────────────────────────────┐
│  LLM Symptom Interpreter                                     │
│  1. Extracts relevant HPO terms with semantic relevance     │
│  2. Outputs dynamic weight vector for annotation columns     │
│     e.g., {"hpo_weight": 0.35, "clinvar_weight": 0.40, ...}  │
└─────────────────────────────────────────────────────────────┘
        │
        ▼
Variant CSV ──→ Dual-Dynamic Scoring Engine ──→ Personalized Ranking
                   │
                   ├── Column weights: LLM-adjusted per symptom
                   └── HPO relevance: semantic similarity to patient
```

---

## 🏗️ Full Module Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         SeekRare                                  │
│                                                                  │
│  ┌────────────────────── Stage 1 ───────────────────────┐         │
│  │          VCF Preprocessing & Annotation               │         │
│  │                                                        │         │
│  │  preprocess/vcf_to_gt.py     VCF → GT CSV            │         │
│  │  preprocess/gene_annotation.py  GTF gene annotation  │         │
│  │  preprocess/clinvar_annotation.py ClinVar merge       │         │
│  │                                                        │         │
│  └──────────────────────────────────────────────────────┘         │
│                            ↓                                    │
│  ┌────────────────────── Stage 2 ───────────────────────┐         │
│  │          LLM-Powered Analysis                       │         │
│  │                                                        │         │
│  │  llm/symptom_parser.py   LLM symptom → HPO + weights│         │
│  │  llm/genos_client.py    Genos 临床遗传学专用 LLM   │         │
│  │  llm/alphafold_client.py AlphaFold 结构预测       │         │
│  │                                                        │         │
│  └──────────────────────────────────────────────────────┘         │
│                            ↓                                    │
│  ┌────────────────────── Stage 3 ───────────────────────┐         │
│  │          Scoring & Ranking                           │         │
│  │                                                        │         │
│  │  scoring/engine.py   Dual-dynamic scoring engine     │         │
│  │  scoring/ranker.py   Personalized ranking            │         │
│  │                                                        │         │
│  └──────────────────────────────────────────────────────┘         │
│                            ↓                                    │
│  ┌────────────────────── Stage 4 ───────────────────────┐         │
│  │          Annotation Extensions                         │         │
│  │                                                        │         │
│  │  annotation/hpo_matcher.py   HPO 语义匹配              │         │
│  │  annotation/combiner.py    多源注释合并               │         │
│  └──────────────────────────────────────────────────────┘         │
└─────────────────────────────────────────────────────────────────┘
```

---

## 🚀 Quick Start

### Python API

```python
from seekrare import SeekRarePipeline

# Initialize pipeline
pipeline = SeekRarePipeline(
    vcf_proband="child.vcf.gz",
    vcf_father="father.vcf.gz",
    vcf_mother="mother.vcf.gz",
    ref_fasta="/ref/GRCh38.fa",
    gtf_file="/ref/genomic.gtf",
    clinvar_csv="/ref/clinvar.csv",
    llm_provider="openai",
    llm_model="gpt-4o",
    api_key=os.getenv("OPENAI_API_KEY"),
)

# Run full pipeline
result = pipeline.run(
    symptoms="Patient presents with intellectual disability, seizures, "
             "hypotonia, and characteristic facial features. "
             "EEG shows generalized spike-wave discharges."
)
result.to_csv("candidate_variants.csv", index=False)
```

### Stage-by-Stage Usage

```python
from seekrare.preprocess import vcf_to_gt_csv, annotate_by_gtf, merge_filter_clinvar
from seekrare.llm import LLMSymptomParser, GenosClient
from seekrare.annotation import HPOMatcher
from seekrare.scoring import DualDynamicScorer, rank_variants

# Stage 1: VCF → Annotated CSV
vcf_to_gt_csv("trio.vcf.gz", "1_gt.csv")
annotate_by_gtf("1_gt.csv", "genomic.gtf", "2_annotated.csv")
merge_filter_clinvar("2_annotated.csv", "clinvar.csv", "3_clinvar.csv")

# Stage 2: HPO Matching
matcher = HPOMatcher(use_ontology=True)
hpo_results = matcher.match("intellectual disability, seizures, hypotonia")

# Stage 3: LLM Interpretation
llm = LLMSymptomParser(provider="openai", model="gpt-4o")
llm_out = llm.interpret("intellectual disability, seizures...")

# Stage 4: Scoring & Ranking
scorer = DualDynamicScorer(weight_vector=llm_out["weight_vector"])
df = pd.read_csv("3_clinvar.csv")
scored = scorer.score(df, relevant_hpos=llm_out["relevant_hpos"])
ranked = rank_variants(scored, top_k=50)
```

---

## 🔬 New Modules (Extended Capabilities)

### Genos LLM — 临床遗传学专用大模型

专门针对临床遗传学的 LLM，可进行综合分析和临床报告生成。

```python
from seekrare.llm import GenosClient

# Initialize Genos
genos = GenosClient(
    api_key="yourgenosapikey",      # Set GENOS_API_KEY env var
    base_url="https://api.genos.tech/v1",
    model="genos-clinical-v1",
)

# Analyze a single variant in clinical context
variant_info = {
    "gene": "SCN1A",
    "cDNA_change": "c.2545C>T",
    "protein_change": "p.Arg849Cys",
    "clinvar_sig": "Pathogenic",
    "cadd_score": 35.0,
    "gnomad_af": 0.00001,
    "impact": "HIGH",
}
analysis = genos.analyze_variant(
    variant_info=variant_info,
    patient_phenotype="Dravet syndrome: febrile seizures, developmental regression",
    hpo_terms=["HP:0001250", "HP:0012469"],
)
print(analysis["pathogenicity_assessment"])

# Batch analyze top candidates
results = genos.batch_analyze(candidate_df, patient_phenotype, top_n=20)

# Generate clinical report
report = genos.generate_clinical_report(
    candidate_variants=results,
    patient_info={"age": "3y", "sex": "M", "family_history": "None"},
    phenotype="febrile seizures, developmental regression",
)
print(report)
```

### AlphaFold Server — 蛋白质结构预测

对候选基因进行蛋白质结构预测，验证变异对蛋白结构的影响。

```python
from seekrare.llm import AlphaFoldClient

af = AlphaFoldClient(output_dir="alphafold_results")

# Submit prediction (non-blocking)
job = af.predict_sequence(
    sequence="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPH...",
    gene_name="HBB",
    wait=True,         # Poll until complete
    poll_interval=30,
)
print(f"PDB URL: {job['pdb_url']}")

# Download PDB
pdb_path = af.download_pdb(job["pdb_url"], "HBB_structure.pdb")

# Quick submit + download in one call
pdb_path = af.predict_and_download(sequence, gene_name="BRCA1")
```

### HPO Matcher — 症状语义匹配

将自由文本症状描述匹配到 HPO ontology 术语，支持 pyhpo 语义相似度。

```python
from seekrare.annotation import HPOMatcher

matcher = HPOMatcher(use_ontology=True)

# Match symptoms to HPO terms
results = matcher.match(
    "intellectual disability, seizures, generalized spike-wave on EEG, hypotonia",
    top_k=15,
)
# Returns: [{"hpo_id": "HP:0001250", "hpo_name": "Seizure", "score": 0.95}, ...]

# Batch matching
batch = matcher.batch_match(
    ["seizures + hypotonia", "cardiomyopathy + ichthyosis"],
    top_k=10,
)
```

---

## 🏗️ Two-Stage Architecture

### Stage 1: VCF Preprocessing & Annotation

```
Trio VCFs (father + mother + child)
        │
        ▼
bcftools preprocessing ────────────────────────────────
  • bcftools norm (left-normalize, split multi-allelics)
  • bcftools merge (trio)
  • bcftools filter (QUAL>30, DP>10, GQ>20)
  • Split by inheritance: de novo / recessive
  • Exclude common dbSNP variants
        │
        ▼
vcf_to_gt_csv.py       → CHROM, POS, REF, ALT, GT per sample
        │
        ▼
annotate_gtf.py        → gene_name, feature_type (NCBI GTF sweep-line)
        │
        ▼
merge_clinvar.py       → clinvar_sig, clinvar_mc, clinvar_min_distance
        │
        ▼
Annotated Variant CSV (ready for LLM scoring)
```

### Stage 2: LLM-Powered Analysis

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
GenosClient (optional) → deep clinical genetics analysis
        │
        ▼
DualDynamicScorer → seekrare_score per variant
        │
        ▼
AlphaFoldClient (optional) → protein structure validation
        │
        ▼
Ranked Candidate Variants (personalized top-K)
```

---

## 📁 Project Structure

```
seekrare/
├── README.md
├── pyproject.toml
├── requirements.txt
├── scripts/
│   ├── bcftools_preprocess.sh        # VCF bcftools 预处理
│   ├── vcf_to_gt_csv.py              # VCF → GT CSV
│   ├── annotate_vcf_csv_by_ncbi_gtf.py # GTF 基因注释
│   └── merge_filter_clinvar_with_distance.py  # ClinVar 注释
└── src/seekrare/
    ├── pipeline.py                   # 全流程编排器
    ├── preprocess/
    │   ├── vcf_to_gt.py
    │   ├── gene_annotation.py
    │   └── clinvar_annotation.py
    ├── annotation/
    │   ├── combiner.py              # 多源注释合并
    │   └── hpo_matcher.py           # HPO 语义匹配
    ├── llm/
    │   ├── symptom_parser.py         # LLM 症状解析
    │   ├── genos_client.py          # Genos 临床 LLM
    │   └── alphafold_client.py      # AlphaFold 结构预测
    └── scoring/
        ├── engine.py                # 双动态打分引擎
        └── ranker.py                # 排序输出
```

---

## ⚙️ Configuration

| Parameter | Description | Default |
|-----------|-------------|---------|
| `gtf_file` | NCBI genomic.gtf 基因注释文件 | Required |
| `clinvar_csv` | ClinVar CSV 致病性注释 | Optional |
| `ref_fasta` | GRCh38 参考基因组 | Required for bcftools |
| `dbSNP_vcf` | dbSNP VCF 常见变异过滤 | Optional |
| `max_af` | Max gnomAD allele frequency | 0.01 |
| `llm_provider` | "openai", "anthropic", "local", "genos" | "openai" |
| `llm_model` | Model name | "gpt-4o" |
| `top_k` | Return top-K candidates | 50 |

---

## 🔬 Scientific Rationale

### Why dual-dynamic?

| Dimension | Traditional | SeekRare |
|-----------|------------|----------|
| **Column weights** | Fixed at pipeline design | **LLM-adjusted per patient** |
| **HPO relevance** | Binary match | **Semantic similarity to patient symptoms** |
| **Output** | One-size-fits-all | **Personalized ranking** |

---

## 📄 License

MIT
