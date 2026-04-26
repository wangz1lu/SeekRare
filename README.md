# SeekRare
**LLM-powered Rare Disease Diagnosis System with Dual-Dynamic Variant Scoring**

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

## 🏗️ Two-Stage Architecture

### Stage 1: VCF Preprocessing & Annotation (bioinformatics pipeline)

```
Trio VCFs (father + mother + child)
        │
        ▼
┌─ bcftools preprocessing ─────────────────────────────┐
│  • bcftools norm (left-normalize, split multi-allelics)
│  • bcftools merge (trio)
│  • bcftools filter (QUAL>30, DP>10, GQ>20)
│  • Split by inheritance: de novo / recessive / etc.
│  • Exclude common dbSNP variants
└──────────────────────────────────────────────────┘
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
DualDynamicScorer → seekrare_score per variant
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
├── .gitignore
├── scripts/
│   ├── vcf_to_gt_csv.py                 # VCF → GT CSV
│   ├── annotate_vcf_csv_by_ncbi_gtf.py  # Gene annotation (GTF)
│   ├── merge_filter_clinvar_with_distance.py  # ClinVar annotation
│   └── bcftools_preprocess.sh            # bcftools preprocessing wrapper
└── src/seekrare/
    ├── pipeline.py         # Full two-stage orchestrator
    ├── preprocess/
    │   ├── __init__.py
    │   ├── vcf_to_gt.py
    │   ├── gene_annotation.py
    │   └── clinvar_annotation.py
    ├── annotation/
    │   └── combiner.py     # General annotation combiner
    ├── llm/
    │   └── symptom_parser.py  # LLM symptom → HPO + weights
    └── scoring/
        ├── engine.py       # Dual-dynamic scoring engine
        └── ranker.py      # Ranking utilities
```

---

## 🚀 Quick Start

```python
from seekrare import SeekRarePipeline

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

result = pipeline.run(
    symptoms="Patient presents with intellectual disability, seizures, "
             "hypotonia, and characteristic facial features. "
             "EEG shows generalized spike-wave discharges."
)

result.to_csv("candidate_variants.csv", index=False)
```

---

## ⚙️ Configuration

Key settings in `SeekRareConfig`:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `gtf_file` | NCBI genomic.gtf for gene annotation | Required |
| `clinvar_csv` | ClinVar CSV for pathogenicity annotation | Optional |
| `ref_fasta` | GRCh38 reference for bcftools normalization | Required for bcftools |
| `dbSNP_vcf` | dbSNP VCF for common SNP exclusion | Optional |
| `max_af` | Max gnomAD allele frequency filter | 0.01 |
| `llm_provider` | "openai", "anthropic", or "local" | "openai" |
| `llm_model` | Model name | "gpt-4o" |
| `top_k` | Number of top candidates to return | 50 |

---

## 🔬 Scientific Rationale

### Why dual-dynamic?

| Dimension | Traditional | SeekRare |
|-----------|------------|----------|
| **Column weights** | Fixed at pipeline design | **LLM-adjusted per patient** |
| **HPO relevance** | Binary match | **Semantic similarity to patient symptoms** |
| **Output** | One-size-fits-all | **Personalized ranking** |

### LLM as Clinical Geneticist

The LLM acts as a clinical geneticist in the loop:
- Reads the patient's actual symptom description
- Knows which phenotypic terms are semantically related
- Can weight ClinVar evidence differently for a seizure phenotype vs. a cardiac phenotype

---

## 📋 Pipeline Steps (Stage 1 Scripts)

```bash
# Step 1a: bcftools preprocessing (bash)
bash scripts/bcftools_preprocess.sh /outdir father.vcf.gz mother.vcf.gz child.vcf.gz

# Step 1b: VCF → GT CSV
python -m seekrare.preprocess.vcf_to_gt input.vcf output.csv

# Step 1c: Gene annotation (GTF sweep-line)
python -m seekrare.preprocess.gene_annotation input.csv genomic.gtf output.csv

# Step 1d: ClinVar annotation + distance filter
python -m seekrare.preprocess.clinvar_annotation annotated.csv clinvar.csv output.csv
```

---

## 📄 License

MIT
