"""
SeekRare Main Pipeline — Full rare disease diagnosis orchestration.

Two-stage design:
  Stage 1 (Preprocess): VCF → GT CSV → Gene Annotation → ClinVar Annotation
  Stage 2 (Analysis):   LLM Symptom Interpretation → Dual-Dynamic Scoring → Rank
"""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger

from seekrare.io.vcf_parser import load_trio_vcf
from seekrare.preprocess.vcf_to_gt import vcf_to_gt_csv
from seekrare.preprocess.gene_annotation import annotate_by_gtf
from seekrare.preprocess.clinvar_annotation import merge_filter_clinvar
from seekrare.llm.symptom_parser import LLMSymptomParser
from seekrare.scoring.engine import DualDynamicScorer
from seekrare.scoring.ranker import rank_variants


@dataclass
class SeekRareConfig:
    """Full configuration for SeekRare."""

    # ── Stage 1: Preprocessing ────────────────────────────────
    vcf_proband: str | Path
    vcf_father: Optional[str | Path] = None
    vcf_mother: Optional[str | Path] = None

    # bcftools preprocessing
    use_bcftools: bool = True
    ref_fasta: Optional[str | Path] = None        # GRCh38 reference
    dbSNP_vcf: Optional[str | Path] = None       # Common SNP exclusion
    inheritance_modes: tuple[str, ...] = ("denovo", "recessive")  # what to call

    # Annotation
    gtf_file: Optional[str | Path] = None
    clinvar_csv: Optional[str | Path] = None

    # ── Stage 2: Analysis ───────────────────────────────────
    llm_provider: str = "openai"
    llm_model: str = "gpt-4o"
    api_key: Optional[str] = None
    base_url: Optional[str] = None
    llm_temperature: float = 0.0

    # Filtering
    max_af: float = 0.01
    min_cadd: Optional[float] = None

    # Output
    top_k: int = 50
    work_dir: str | Path = "seekrare_output"


class SeekRarePipeline:
    """
    Full SeekRare pipeline: Preprocess → Annotate → LLM Score → Rank

    Example:
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
            symptoms="Intellectual disability, seizures, hypotonia, "
                     "characteristic facial features, generalized spike-wave on EEG."
        )
        result.to_csv("candidates.csv", index=False)
    """

    def __init__(self, config: SeekRareConfig):
        self.cfg = config
        self.work_dir = Path(config.work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self._init_llm()

    def _init_llm(self) -> None:
        self.llm = LLMSymptomParser(
            provider=self.cfg.llm_provider,
            model=self.cfg.llm_model,
            api_key=self.cfg.api_key,
            base_url=self.cfg.base_url,
            temperature=self.cfg.llm_temperature,
        )
        logger.info(f"LLM ready: {self.cfg.llm_provider}/{self.cfg.llm_model}")

    # ── Stage 1: Preprocessing ──────────────────────────────────────────────

    def stage1_preprocess(self) -> pd.DataFrame:
        """
        Run Stage 1: VCF → GT CSV → Gene Annotation → ClinVar Annotation.

        Returns
        -------
        pd.DataFrame
            Annotated variants DataFrame with columns from all annotation steps.
        """
        logger.info("=" * 60)
        logger.info("Stage 1: VCF Preprocessing & Annotation")
        logger.info("=" * 60)

        wd = self.work_dir

        # 1a. bcftools preprocessing
        if self.cfg.use_bcftools and self.cfg.ref_fasta:
            logger.info("Step 1a: bcftools preprocessing...")
            vcf_trio = self._run_bcftools_preprocess()
        else:
            logger.info("Step 1a: Skipping bcftools (using direct VCF loading)...")
            vcf_trio = self._load_direct_vcf()

        # 1b. VCF → GT CSV
        gt_csv = wd / "1_gt.csv"
        logger.info(f"Step 1b: VCF → GT CSV: {gt_csv}")
        vcf_to_gt_csv(vcf_trio, gt_csv)

        # 1c. Gene annotation (NCBI GTF)
        annot_csv = wd / "2_annotated.csv"
        if self.cfg.gtf_file:
            logger.info(f"Step 1c: Gene annotation (GTF): {annot_csv}")
            annotate_by_gtf(gt_csv, self.cfg.gtf_file, annot_csv)
        else:
            raise FileNotFoundError("gtf_file is required for gene annotation")

        # 1d. ClinVar annotation
        clinvar_csv = wd / "3_clinvar.csv"
        if self.cfg.clinvar_csv:
            logger.info(f"Step 1d: ClinVar annotation: {clinvar_csv}")
            df = merge_filter_clinvar(annot_csv, self.cfg.clinvar_csv, clinvar_csv)
        else:
            logger.warning("clinvar_csv not provided, skipping ClinVar annotation")
            df = pd.read_csv(annot_csv)

        logger.info(f"Stage 1 complete: {len(df):,} variants, columns: {list(df.columns)}")
        return df

    def _run_bcftools_preprocess(self) -> Path:
        """Run bcftools preprocessing pipeline, return path to filtered VCF."""
        wd = self.work_dir / "bcftools"
        wd.mkdir(parents=True, exist_ok=True)
        norm_dir = wd / "norm"
        merge_dir = wd / "merge"
        filt_dir = wd / "filtered"
        for d in [norm_dir, merge_dir, filt_dir]:
            d.mkdir(exist_ok=True)

        ref = str(self.cfg.ref_fasta)
        father = str(self.cfg.vcf_father)
        mother = str(self.cfg.vcf_mother)
        child = str(self.cfg.vcf_proband)

        # Normalize each sample
        for vcf, label in [(father, "father"), (mother, "mother"), (child, "child")]:
            out = norm_dir / f"{label}.norm.vcf.gz"
            if not out.exists():
                subprocess.run(
                    ["bcftools", "norm", "-m", "-both", "-f", ref, vcf, "-Oz", "-o", str(out)],
                    check=True, capture_output=True,
                )
                subprocess.run(["bcftools", "index", str(out)], check=True, capture_output=True)
            logger.info(f"  Normalized {label}: {out}")

        # Merge trio
        trio_vcf = merge_dir / "trio.vcf.gz"
        if not trio_vcf.exists():
            subprocess.run(
                ["bcftools", "merge",
                 str(norm_dir / "father.norm.vcf.gz"),
                 str(norm_dir / "mother.norm.vcf.gz"),
                 str(norm_dir / "child.norm.vcf.gz"),
                 "-Oz", "-o", str(trio_vcf)],
                check=True, capture_output=True,
            )
            subprocess.run(["bcftools", "index", str(trio_vcf)], check=True, capture_output=True)
        logger.info(f"  Merged trio: {trio_vcf}")

        # Quality filter
        filt_vcf = filt_dir / "filtered.vcf.gz"
        if not filt_vcf.exists():
            subprocess.run(
                ["bcftools", "filter",
                 "-i", "QUAL>30 && FMT/DP>10 && FMT/GQ>20",
                 str(trio_vcf), "-Oz", "-o", str(filt_vcf)],
                check=True, capture_output=True,
            )
            subprocess.run(["bcftools", "index", str(filt_vcf)], check=True, capture_output=True)
        logger.info(f"  Quality filtered: {filt_vcf}")

        return filt_vcf

    def _load_direct_vcf(self) -> Path:
        """Load VCF directly without bcftools."""
        return Path(self.cfg.vcf_proband)

    # ── Stage 2: LLM + Scoring ───────────────────────────────────────────

    def stage2_analyze(self, df: pd.DataFrame, symptoms: str) -> pd.DataFrame:
        """Run Stage 2: LLM interpretation + dual-dynamic scoring + ranking."""
        logger.info("=" * 60)
        logger.info("Stage 2: LLM Analysis & Dual-Dynamic Scoring")
        logger.info("=" * 60)

        # Step 2a: LLM symptom interpretation
        logger.info("Step 2a: LLM interpreting symptoms...")
        llm_out = self.llm.interpret(symptoms)
        logger.info(f"  HPO terms: {len(llm_out['relevant_hpos'])}, "
                    f"weights: {llm_out['weight_vector']}")

        # Step 2b: Dual-dynamic scoring
        logger.info("Step 2b: Dual-dynamic scoring...")
        scorer = DualDynamicScorer(weight_vector=llm_out["weight_vector"])
        scored = scorer.score(df, relevant_hpos=llm_out["relevant_hpos"])

        # Step 2c: Filter
        if "gnomad_af" in scored.columns:
            scored = scored[scored["gnomad_af"].fillna(0) <= self.cfg.max_af]

        # Step 2d: Rank
        logger.info("Step 2c: Ranking and exporting...")
        ranked = rank_variants(scored, top_k=self.cfg.top_k)

        out_path = self.work_dir / "candidate_variants.csv"
        ranked.to_csv(out_path, index=False)
        logger.info(f"  Results saved: {out_path}")

        return ranked

    def run(self, symptoms: str) -> pd.DataFrame:
        """
        Run the full SeekRare pipeline.

        Parameters
        ----------
        symptoms : str
            Free-text patient symptom description

        Returns
        -------
        pd.DataFrame
            Ranked candidate variants with all annotation columns + seekrare_score
        """
        df = self.stage1_preprocess()
        return self.stage2_analyze(df, symptoms)
