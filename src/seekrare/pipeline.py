"""
SeekRare Main Pipeline Orchestrator

Coordinates the full rare disease diagnosis pipeline:
1. Load VCF (trio: proband + parents)
2. Annotate variants (ClinVar, HPO, CADD, gnomAD, etc.)
3. LLM symptom interpretation (HPO terms + weight vector)
4. Dual-dynamic scoring engine
5. Rank and export
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd
from loguru import logger

from seekrare.io.vcf_parser import load_trio_vcf
from seekrare.annotation.combiner import annotate_variants
from seekrare.llm.symptom_parser import LLMSymptomParser
from seekrare.scoring.engine import DualDynamicScorer
from seekrare.scoring.ranker import rank_variants


@dataclass
class SeekRareConfig:
    """Configuration for SeekRare pipeline."""

    # VCF inputs
    vcf_proband: str | Path
    vcf_father: Optional[str | Path] = None
    vcf_mother: Optional[str | Path] = None

    # LLM settings
    llm_provider: str = "openai"          # "openai" | "anthropic" | "local"
    llm_model: str = "gpt-4o"
    api_key: Optional[str] = None
    base_url: Optional[str] = None         # For OpenAI-compatible APIs
    llm_temperature: float = 0.0

    # Annotation
    clinvar_vcf: Optional[str | Path] = None
    vep_cache_dir: Optional[str | Path] = None

    # Filtering
    max_af: float = 0.01                 # Max gnomAD allele frequency
    filter_synonymous: bool = True
    min_cadd: Optional[float] = None

    # Scoring
    top_k: int = 50

    # Output
    output_dir: str | Path = "output"


class SeekRarePipeline:
    """
    Main SeekRare pipeline.

    Example:
        pipeline = SeekRarePipeline(
            vcf_proband="proband.vcf",
            vcf_father="father.vcf",
            vcf_mother="mother.vcf",
            llm_provider="openai",
            llm_model="gpt-4o",
            api_key=os.getenv("OPENAI_API_KEY"),
        )
        result = pipeline.run(
            symptoms="Intellectual disability, seizures, hypotonia, "
                     "characteristic facial features, generalized spike-wave on EEG."
        )
        result.to_csv("output/candidates.csv", index=False)
    """

    def __init__(self, config: SeekRareConfig):
        self.config = config
        self._setup_logger()
        self._init_llm_parser()

    def _setup_logger(self) -> None:
        logger.add(
            os.path.join(str(self.config.output_dir), "seekrare.log"),
            rotation="10 MB",
            level="INFO",
        )

    def _init_llm_parser(self) -> None:
        self.llm_parser = LLMSymptomParser(
            provider=self.config.llm_provider,
            model=self.config.llm_model,
            api_key=self.config.api_key,
            base_url=self.config.base_url,
            temperature=self.config.llm_temperature,
        )
        logger.info(f"LLM parser initialized: {self.config.llm_provider}/{self.config.llm_model}")

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
            Ranked candidate variants, one row per variant,
            with annotation columns and final `seekrare_score` column.
        """
        logger.info("=" * 60)
        logger.info("SeekRare Pipeline Started")
        logger.info("=" * 60)

        # ── Step 1: Load VCF ────────────────────────────────────────
        logger.info("Step 1: Loading VCF files...")
        variants_df = load_trio_vcf(
            proband_vcf=self.config.vcf_proband,
            father_vcf=self.config.vcf_father,
            mother_vcf=self.config.vcf_mother,
        )
        logger.info(f"  Raw variants: {len(variants_df)}")

        # ── Step 2: Annotate ──────────────────────────────────────
        logger.info("Step 2: Annotating variants...")
        annotated_df = annotate_variants(
            variants_df,
            clinvar_vcf=self.config.clinvar_vcf,
            vep_cache_dir=self.config.vep_cache_dir,
        )
        logger.info(f"  Annotated variants: {len(annotated_df)}")

        # ── Step 3: Filter ────────────────────────────────────────
        logger.info("Step 3: Filtering variants...")
        filtered_df = self._filter_variants(annotated_df)
        logger.info(f"  After filter: {len(filtered_df)}")

        if len(filtered_df) == 0:
            logger.warning("No variants pass filters. Returning empty DataFrame.")
            return filtered_df

        # ── Step 4: LLM Symptom Interpretation ─────────────────────
        logger.info("Step 4: LLM interpreting symptoms...")
        llm_output = self.llm_parser.interpret(symptoms)
        logger.info(f"  LLM output: {llm_output}")

        # ── Step 5: Dual-Dynamic Scoring ───────────────────────────
        logger.info("Step 5: Running dual-dynamic scoring...")
        scorer = DualDynamicScorer(weight_vector=llm_output["weight_vector"])
        scored_df = scorer.score(filtered_df, relevant_hpos=llm_output["relevant_hpos"])
        logger.info(f"  Scored variants: {len(scored_df)}")

        # ── Step 6: Rank and Export ────────────────────────────────
        logger.info("Step 6: Ranking and exporting...")
        ranked_df = rank_variants(scored_df, top_k=self.config.top_k)

        os.makedirs(self.config.output_dir, exist_ok=True)
        out_path = os.path.join(self.config.output_dir, "candidate_variants.csv")
        ranked_df.to_csv(out_path, index=False)
        logger.info(f"  Results saved to: {out_path}")

        logger.info("=" * 60)
        logger.info("SeekRare Pipeline Completed")
        logger.info("=" * 60)

        return ranked_df

    def _filter_variants(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply quality and frequency filters."""
        mask = pd.Series(True, index=df.index)

        if self.config.max_af < 1.0 and "gnomad_af" in df.columns:
            mask &= df["gnomad_af"].fillna(0) <= self.config.max_af

        if self.config.filter_synonymous and "impact" in df.columns:
            mask &= df["impact"].isin(["HIGH", "MODERATE", "LOW"])

        if self.config.min_cadd is not None and "cadd_score" in df.columns:
            mask &= df["cadd_score"].fillna(0) >= self.config.min_cadd

        return df[mask].reset_index(drop=True)
