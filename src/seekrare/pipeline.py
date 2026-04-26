"""
SeekRare — full two-stage rare disease diagnosis pipeline.

Stage 1: VCF preprocessing & annotation (bioinformatics)
  bcftools → VCF→GT → GTF annotation → VEP/CADD/SpliceAI/ClinVar/OMIM annotation

Stage 2: LLM-powered analysis (AI)
  LLM symptom interpretation → Dual-dynamic scoring → Model analysis → Report

Usage:
    from seekrare import SeekRarePipeline

    pipeline = SeekRarePipeline(
        vcf_proband="child.vcf.gz",
        vcf_father="father.vcf.gz",
        vcf_mother="mother.vcf.gz",
        ref_fasta="/ref/GRCh38.fa",
        gtf_file="/ref/genomic.gtf",
        vep_cache="/ref/vep_cache",
        clinvar_csv="/ref/clinvar.vcf.gz",
        cadd_tsv="/ref/cadd.tsv.gz",
        omim_dir="/ref/omim",
        llm_provider="openai",
        llm_model="gpt-4o",
        api_key=os.getenv("OPENAI_API_KEY"),
    )

    result = pipeline.run(symptoms="intellectual disability, seizures...")
    result.to_csv("candidates.csv", index=False)
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger

from seekrare.preprocessing import (
    vcf_to_gt_csv,
    annotate_by_gtf,
    merge_filter_clinvar,
)
from seekrare.annotation import (
    VEPAnnotationLoader,
    ClinVarLoader,
    OMIMLoader,
    CADDLoader,
    SpliceAILoader,
    HPOMatcher,
)
from seekrare.llm import LLMSymptomParser, GenosClient
from seekrare.scoring import (
    DualDynamicScorer,
    rank_variants,
    ModelAnalyzer,
    ModelAnalyzerConfig,
)


@dataclass
class SeekRareConfig:
    """SeekRare pipeline configuration."""

    # Stage 1: VCF preprocessing
    vcf_proband: str
    vcf_father: Optional[str] = None
    vcf_mother: Optional[str] = None
    ref_fasta: Optional[str] = None
    gtf_file: Optional[str] = None

    # Annotation resources
    vep_cache: Optional[str] = None
    clinvar_vcf: Optional[str] = None
    cadd_tsv: Optional[str] = None
    omim_dir: Optional[str] = None
    spliceai_vcf: Optional[str] = None
    dbscsnv_txt: Optional[str] = None

    # Filtering thresholds
    max_af: float = 0.01
    min_quality: float = 30.0

    # Stage 2: LLM
    llm_provider: str = "openai"
    llm_model: str = "gpt-4o"
    api_key: Optional[str] = None
    base_url: Optional[str] = None

    # Scoring
    top_k: int = 50

    # Model analysis (after scoring)
    enable_model_analysis: bool = False
    model_analysis_top_k: int = 10
    model_analysis_language: str = "cn"   # "cn" or "en"

    # Output
    output_dir: str = "seekrare_output"

    def __post_init__(self):
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)


class SeekRarePipeline:
    """
    Full two-stage rare disease diagnosis pipeline.

    Stage 1: VCF → Annotated CSV
    Stage 2: LLM Interpretation → Scoring → Model Analysis → Report
    """

    def __init__(self, config: Union[SeekRareConfig, dict]):
        if isinstance(config, dict):
            config = SeekRareConfig(**config)
        self.config = config

        self._stage1_df: Optional[pd.DataFrame] = None
        self._llm_interpretation: Optional[dict] = None

        # Init components lazily
        self._llm_parser: Optional[LLMSymptomParser] = None
        self._hpo_matcher: Optional[HPOMatcher] = None
        self._model_analyzer: Optional[ModelAnalyzer] = None

    # ── Stage 1 ────────────────────────────────────────────────────────────────

    def stage1_preprocess(self) -> pd.DataFrame:
        """
        Run Stage 1: VCF preprocessing + multi-source annotation.

        Returns
        -------
        pd.DataFrame
            Fully annotated variants
        """
        cfg = self.config
        out = cfg.output_dir
        logger.info("=== Stage 1: VCF Preprocessing & Annotation ===")

        # Step 1: VCF → GT CSV
        gt_csv = out / "1_gt.csv"
        if not gt_csv.exists():
            vcf_to_gt_csv(cfg.vcf_proband, str(gt_csv))
        else:
            logger.info(f"Using existing: {gt_csv}")

        df = pd.read_csv(gt_csv, dtype=str)
        df["POS"] = df["POS"].astype(int)
        logger.info(f"Loaded {len(df)} variants from GT CSV")

        # Step 2: GTF gene annotation
        gtf_annotated = out / "2_gtf_annotated.csv"
        if cfg.gtf_file and not gtf_annotated.exists():
            annotate_by_gtf(str(gt_csv), cfg.gtf_file, str(gtf_annotated))
        else:
            logger.info(f"Using existing: {gtf_annotated}")

        df = pd.read_csv(str(gtf_annotated), dtype=str) if gtf_annotated.exists() else df
        logger.info(f"After GTF: {len(df)} variants")

        # Step 3: ClinVar annotation
        clinvar_annotated = out / "3_clinvar.csv"
        if cfg.clinvar_vcf and not clinvar_annotated.exists():
            merge_filter_clinvar(str(gtf_annotated), cfg.clinvar_vcf, str(clinvar_annotated))
        else:
            logger.info(f"Using existing: {clinvar_annotated}")

        df = pd.read_csv(str(clinvar_annotated), dtype=str) if clinvar_annotated.exists() else df
        logger.info(f"After ClinVar: {len(df)} variants")

        # Step 4: VEP annotation
        if cfg.vep_cache:
            vep_annotated = out / "4_vep.csv"
            if not vep_annotated.exists():
                vep = VEPAnnotationLoader(cfg.vep_cache)
                vep.build_index()
                df = vep.annotate_variants(df)
                df.to_csv(vep_annotated, index=False)
            else:
                df = pd.read_csv(str(vep_annotated), dtype=str)
                logger.info(f"Using existing VEP: {vep_annotated}")
            logger.info(f"After VEP: {len(df)} variants")

        # Step 5: CADD annotation
        if cfg.cadd_tsv:
            cadd_annotated = out / "5_cadd.csv"
            if not cadd_annotated.exists():
                cadd = CADDLoader(cfg.cadd_tsv)
                df = cadd.annotate_variants(df)
                df.to_csv(cadd_annotated, index=False)
            else:
                df = pd.read_csv(str(cadd_annotated), dtype=str)
                logger.info(f"Using existing CADD: {cadd_annotated}")
            logger.info(f"After CADD: {len(df)} variants")

        # Step 6: SpliceAI annotation
        if cfg.spliceai_vcf:
            splice_annotated = out / "6_spliceai.csv"
            if not splice_annotated.exists():
                splice = SpliceAILoader(cfg.spliceai_vcf)
                df = splice.annotate_variants(df)
                df.to_csv(splice_annotated, index=False)
            else:
                df = pd.read_csv(str(splice_annotated), dtype=str)
                logger.info(f"Using existing SpliceAI: {splice_annotated}")
            logger.info(f"After SpliceAI: {len(df)} variants")

        # Step 7: OMIM annotation
        if cfg.omim_dir:
            omim_annotated = out / "7_omim.csv"
            if not omim_annotated.exists():
                omim = OMIMLoader(cfg.omim_dir)
                omim.load()
                df = omim.annotate_genes(df)
                df.to_csv(omim_annotated, index=False)
            else:
                df = pd.read_csv(str(omim_annotated), dtype=str)
                logger.info(f"Using existing OMIM: {omim_annotated}")
            logger.info(f"After OMIM: {len(df)} variants")

        # Step 8: Allele frequency filter
        if "gnomad_af" in df.columns:
            n_before = len(df)
            df["gnomad_af"] = pd.to_numeric(df["gnomad_af"], errors="coerce")
            df = df[df["gnomad_af"].isna() | (df["gnomad_af"] <= cfg.max_af)]
            logger.info(f"AF filter: {n_before} → {len(df)} (max_af={cfg.max_af})")

        self._stage1_df = df
        logger.info(f"Stage 1 complete: {len(df)} annotated variants")
        return df

    # ── Stage 2 ───────────────────────────────────────────────────────────────

    def stage2_analyze(
        self,
        symptoms: str,
        hpo_terms: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        """
        Run Stage 2: LLM interpretation → scoring → ranking.

        Parameters
        ----------
        symptoms : str
            Patient phenotype description
        hpo_terms : list[str], optional
            Pre-identified HPO terms (auto-discovery also attempted)

        Returns
        -------
        pd.DataFrame
            Scored and ranked candidate variants
        """
        logger.info("=== Stage 2: LLM-Powered Analysis ===")

        df = self._stage1_df
        if df is None:
            raise RuntimeError("Stage 1 not run. Call stage1_preprocess() first.")

        # ── LLM Symptom Interpretation ────────────────────────────────────────
        if self._llm_parser is None:
            self._llm_parser = LLMSymptomParser(
                provider=self.config.llm_provider,
                model=self.config.llm_model,
                api_key=self.config.api_key,
                base_url=self.config.base_url,
            )

        llm_out = self._llm_parser.interpret(symptoms)
        self._llm_interpretation = llm_out
        logger.info(f"LLM weights: {llm_out.get('weight_vector', {})}")
        logger.info(f"Relevant HPOs: {[h['hpo_id'] for h in llm_out.get('relevant_hpos', [])]}")

        # ── HPO Matcher (optional augmentation) ───────────────────────────────
        relevant_hpos = llm_out.get("relevant_hpos", [])
        if hpo_terms:
            # Merge user-provided HPOs with LLM output
            hpo_ids = {h["hpo_id"] for h in relevant_hpos}
            for hpo_id in hpo_terms:
                if hpo_id not in hpo_ids:
                    relevant_hpos.append({"hpo_id": hpo_id, "score": 0.9})

        # ── Dual-Dynamic Scoring ──────────────────────────────────────────────
        weight_vector = llm_out.get("weight_vector", {})
        if not weight_vector:
            # Fallback default weights
            weight_vector = {
                "clinvar_sig": 0.3,
                "cadd_score": 0.2,
                "impact": 0.15,
                "gnomad_af": 0.2,
                "hpo_terms": 0.15,
            }

        scorer = DualDynamicScorer(weight_vector)
        df = scorer.score(df, relevant_hpos=relevant_hpos)

        # ── Ranking ────────────────────────────────────────────────────────────
        df = rank_variants(df, top_k=self.config.top_k)

        logger.info(f"Stage 2 complete: top {len(df)} candidates")
        return df

    def run(
        self,
        symptoms: str,
        hpo_terms: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        """
        Run full two-stage pipeline.

        Parameters
        ----------
        symptoms : str
            Patient phenotype description
        hpo_terms : list[str], optional
            Pre-identified HPO terms

        Returns
        -------
        pd.DataFrame
            Final ranked candidate variants
        """
        self.stage1_preprocess()
        df = self.stage2_analyze(symptoms, hpo_terms)

        # ── Model Analysis (optional) ───────────────────────────────────────
        if self.config.enable_model_analysis:
            df = self._run_model_analysis(df, symptoms, hpo_terms)

        # Save output
        out_path = self.config.output_dir / "final_candidates.csv"
        df.to_csv(out_path, index=False)
        logger.info(f"Results saved to {out_path}")

        return df

    def _run_model_analysis(
        self,
        df: pd.DataFrame,
        symptoms: str,
        hpo_terms: Optional[list[str]],
    ) -> pd.DataFrame:
        """Run LLM model analysis on top candidates."""
        logger.info("=== Model Analysis (Stage 2b) ===")

        if self._model_analyzer is None:
            self._model_analyzer = ModelAnalyzer(ModelAnalyzerConfig(
                provider=self.config.llm_provider,
                model=self.config.llm_model,
                api_key=self.config.api_key,
                base_url=self.config.base_url,
                max_tokens=8192,
            ))

        analysis = self._model_analyzer.analyze(
            variants_df=df,
            patient_phenotype=symptoms,
            hpo_terms=hpo_terms,
            top_k=self.config.model_analysis_top_k,
            language=self.config.model_analysis_language,
        )

        # Attach model analysis to DataFrame
        if analysis.get("candidates"):
            acmg_classes = {}
            for cand in analysis["candidates"]:
                rank = cand.get("rank", 0)
                acmg_classes[rank] = cand.get("acmg_class", "VUS")

            df = df.copy()
            df["model_acmg_class"] = df.index.map(
                lambda i: acmg_classes.get(i + 1, "VUS")
            )

        # Save analysis report
        import json
        report_path = self.config.output_dir / "model_analysis_report.json"
        with open(report_path, "w", encoding="utf-8") as f:
            json.dump(analysis, f, ensure_ascii=False, indent=2)
        logger.info(f"Model analysis report saved to {report_path}")

        return df
