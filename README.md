# SeekRare

**Rare Disease Diagnosis System powered by LLM**

> Given a patient's symptom description, family VCF files, and a novel dual-dynamic scoring pipeline, SeekRare delivers personalized, ranked candidate variants — tailored to each patient's unique clinical presentation.

---

## 🎯 Core Idea

Traditional variant prioritization tools apply **fixed, static weights** to annotation columns (ClinVar, CADD, HPO, gnomAD, etc.). This is a poor fit for rare disease diagnostics, where the *same* symptom can have wildly different phenotypic interpretations depending on the patient.

SeekRare solves this with **LLM-driven dual-dynamic scoring**:

```
Patient Symptoms (free text)
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  LLM Symptom Interpreter                                 │
│  1. Extracts relevant HPO terms (semantic search)       │
│  2. Outputs dynamic weight vector W = [w_clinvar, w_hpo, │
│     w_cadd, w_gnomad, ...]                              │
│     dynamically calibrated for this patient                │
└─────────────────────────────────────────────────────────┘
        │
        ▼
Variant CSV (rows = variants, cols = annotations)
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  Dynamic Scoring Engine                                  │
│  score(v) = Σ  W[col] × semantic_similarity(HPO[col], │
│                        relevant_HPOs) × normalize(val)   │
│                                                         │
│  → Column weights change per patient                     │
│  → HPO relevance changes per symptom                   │
│  → Both dimensions are dynamic ("dual-dynamic")          │
└─────────────────────────────────────────────────────────┘
        │
        ▼
Ranked Candidate Variants (personalized to patient)
```

---

## 📂 Project Structure

```
seekrare/
├── README.md
├── requirements.txt
├── pyproject.toml
├── .gitignore
├── src/
│   └── seekrare/
│       ├── __init__.py
│       ├── pipeline.py          # Main orchestrator
│       ├── io/
│       │   ├── vcf_parser.py    # VCF loading (trio supported)
│       │   └── csv_writer.py    # Results export
│       ├── annotation/
│       │   ├── clinvar.py       # ClinVar annotation
│       │   ├── vep.py           # VEP-style annotation
│       │   ├── hpo_matcher.py   # HPO ↔ symptom semantic matching
│       │   └── combiner.py      # Merge all annotations into CSV
│       ├── llm/
│       │   ├── symptom_parser.py # LLM symptom → HPO + weight vector
│       │   └── scorer.py        # LLM-instructed scoring planner
│       └── scoring/
│           ├── engine.py        # Dual-dynamic scoring engine
│           └── ranker.py        # Ranking and output
└── tests/
    └── ...
```

---

## 🔑 Key Concepts

### Dual-Dynamic Scoring

| Dimension | Static (traditional) | Dynamic (SeekRare) |
|-----------|-------------------|-------------------|
| **Column weights** | Fixed at pipeline design | **LLM-adjusted per symptom** |
| **HPO relevance** | Binary match | **Semantic similarity to patient symptoms** |
| **Output** | One-size-fits-all ranking | **Personalized to this patient** |

### LLM Symptom Interpreter

The LLM receives:
- Patient's free-text symptom description
- List of available HPO terms (ontology)

The LLM outputs:
1. **Relevant HPO terms** with semantic relevance scores (0–1)
2. **Annotation column weights** — e.g., `"hpo_weight": 0.35, "clinvar_weight": 0.40, "cadd_weight": 0.15, ...`

### Variant CSV Schema (input to scoring engine)

| column | description |
|--------|-------------|
| `chrom` | Chromosome |
| `pos` | Position |
| `ref` | Reference allele |
| `alt` | Alternate allele |
| `gene` | Gene symbol |
| `clinvar_significance` | ClinVar pathogenicity |
| `clinvar_stars` | ClinVar review status |
| `cadd_score` | CADD phred score |
| `gnomad_af` | gnomAD allele frequency |
| `sift_score` | SIFT prediction |
| `polyphen_score` | PolyPhen prediction |
| `hpo_terms` | Associated HPO terms (pipe-separated) |
| `hgvs_c` | HGVS coding change |
| `impact` | Variant impact (HIGH/MODERATE/LOW) |
| ... | (extensible) |

---

## 🚀 Quick Start

```python
from seekrare import SeekRarePipeline

pipeline = SeekRarePipeline(
    vcf_proband="data/proband.vcf",
    vcf_father="data/father.vcf",
    vcf_mother="data/mother.vcf",
    llm_provider="openai",        # or "anthropic", "local"
    llm_model="gpt-4o",
    api_key=os.getenv("OPENAI_API_KEY"),
)

result = pipeline.run(
    symptoms="Patient presents with intellectual disability, seizures, "
             "hypotonia, and characteristic facial features. EEG shows "
             "generalized spike-wave discharges.",
    top_k=50,                    # Return top 50 candidate variants
)

result.to_csv("output/candidate_variants.csv", index=False)
```

---

## ⚙️ Configuration

See `config/default.yaml` for full configuration options.

Key settings:
- `annotation.clinvar_cache`: Path to ClinVar VCF cache
- `annotation.vep_cache`: VEP annotation cache
- `llm.temperature`: LLM sampling temperature
- `scoring.normalization`: Per-column normalization strategy
- `filtering.max_af`: Max gnomAD allele frequency (default: 0.01)

---

## 📋 Pipeline Steps

```
1. load_vcf()          — Parse and merge trio VCF files
2. annotate()           — ClinVar, VEP, HPO annotation
3. filter_variants()   — Quality + frequency filters
4. llm_interpret()     — LLM: symptom → HPO + weights
5. score_variants()    — Dual-dynamic scoring engine
6. rank_and_export()   — Sort and export CSV
```

---

## 🔬 Scientific Rationale

### Why LLM for scoring?

Rare disease diagnostics requires integrating **heterogeneous, high-dimensional annotation data** — many columns, each capturing different aspects of pathogenicity. Static weight schemes (e.g., ACMG criteria) are widely applicable but fail to capture **patient-specific phenotypic context**.

The LLM acts as a "clinical geneticist in the loop":
- It reads the patient's actual symptoms
- It knows which phenotypic terms are semantically related
- It can weight ClinVar evidence differently for a seizure phenotype vs. a cardiac phenotype

### Why dual-dynamic?

A single dynamic adjustment (either weights or HPO matching) still leaves one dimension static. True personalization requires both:
1. **Which annotation columns matter most** (LLM weight allocation)
2. **Which HPO terms are clinically relevant** (semantic symptom matching)

---

## 📄 License

MIT
