"""
pipeline.py — SeekRare 全流程编排器

三阶段 + Stage 4 架构:
    Stage 1: VCF 家系预处理 → 基本注释（必须）
    Stage 2: GTEx eQTL 高级注释（可选）
    Stage 3: LLM 排序
    Stage 4: Genos 模型分析（可选）

Usage:
    from seekrare import SeekRarePipeline
    pipeline = SeekRarePipeline(vcf_proband="child.vcf.gz", ...)
    pipeline.stage1_preprocess()
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger

from seekrare.preprocess import (
    stage1_vcf_to_gt_csv,
    stage1_annotate_by_gtf,
    stage1_merge_filter_clinvar,
    stage1_dbsnp_filter,
    stage2_eqtl_annotation,
    stage4_genos_analysis,
)


# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class SeekRareConfig:
    """SeekRare 全局配置。所有文件路径由用户配置。"""
    vcf_proband: Union[str, Path]
    vcf_father: Optional[Union[str, Path]] = None
    vcf_mother: Optional[Union[str, Path]] = None
    ref_fasta: Optional[Union[str, Path]] = None
    gtf_file: Optional[Union[str, Path]] = None
    clinvar_vcf: Optional[Union[str, Path]] = None
    dbSNP_vcf: Optional[Union[str, Path]] = None
    work_dir: Union[str, Path] = "seekrare_output"

    # Stage 2
    gtex_tissue_dir: Optional[Union[str, Path]] = None

    # Stage 3
    llm_provider: str = "openai"
    llm_model: str = "deepseek-v4-flash"
    api_key: Optional[str] = None
    base_url: Optional[str] = None

    # Stage 4
    genome_fa: Optional[Union[str, Path]] = None
    genos_model_path: Optional[Union[str, Path]] = None

    def __post_init__(self):
        self.work_dir = Path(self.work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# SeekRarePipeline
# ─────────────────────────────────────────────────────────────────────────────

class SeekRarePipeline:
    """
    SeekRare 三阶段流水线。

    Example:
        from seekrare import SeekRarePipeline

        p = SeekRarePipeline(
            vcf_proband="child.vcf.gz",
            vcf_father="father.vcf.gz",
            vcf_mother="mother.vcf.gz",
            ref_fasta="/ref/GRCh38.fa",
            gtf_file="/ref/genomic.gtf",
            clinvar_vcf="/ref/clinvar.vcf.gz",
            dbSNP_vcf="/ref/dbsnp.vcf.gz",
        )
        p.stage1_preprocess()
        p.stage2_eqtl_annotation(symptoms="眼部病变")
        p.stage3_score_and_rank(symptoms="视网膜色素变性")
        p.stage4_genos_analysis(top_n=5)
    """

    def __init__(self, config: Union[SeekRareConfig, dict]):
        if isinstance(config, dict):
            config = SeekRareConfig(**config)
        self.cfg = config
        self._stage1_df: Optional[pd.DataFrame] = None
        self._llm_interpretation: Optional[dict] = None

    # ── Stage 1 ──────────────────────────────────────────────────────────

    def stage1_preprocess(self) -> pd.DataFrame:
        """
        Stage 1: VCF → 基本注释 CSV。

        流程:
        1. dbSNP common 过滤（如果有 dbSNP_vcf）
        2. VCF → GT CSV（CHROM/POS/REF/ALT/GT）
        3. GTF gene annotation（gene_name, feature_type）
        4. ClinVar 注释（CLNDISDB, CLNDN, CLNREVSTAT, CLNSIG, CLNVC, ORIGIN）

        Returns
        -------
        pd.DataFrame: 基本注释完成的 variants
        """
        logger.info("=" * 60)
        logger.info("Stage 1: VCF Preprocessing & Basic Annotation")
        logger.info("=" * 60)

        cfg = self.cfg
        wd = cfg.work_dir

        # dbSNP 过滤
        proband_vcf = cfg.vcf_proband
        if cfg.dbSNP_vcf and cfg.ref_fasta:
            dbSNP_filtered = wd / "proband.nocommon.vcf.gz"
            if not dbSNP_filtered.exists():
                logger.info(f"Step 1: dbSNP common 过滤 → {dbSNP_filtered}")
                result = stage1_dbsnp_filter(
                    input_vcf=str(proband_vcf),
                    dbsnp_vcf=str(cfg.dbSNP_vcf),
                    output_vcf=str(dbSNP_filtered),
                )
                logger.info(f"  {result['n_removed']}/{result['n_before']} removed")
                proband_vcf = str(dbSNP_filtered)
            else:
                proband_vcf = str(dbSNP_filtered)
                logger.info(f"  [跳过] {dbSNP_filtered} 已存在")

        # VCF → GT CSV
        gt_csv = wd / "1_gt.csv"
        if not gt_csv.exists():
            logger.info(f"Step 2: VCF → GT CSV")
            stage1_vcf_to_gt_csv(str(proband_vcf), str(gt_csv))
        df = pd.read_csv(gt_csv, dtype=str)
        logger.info(f"  {len(df)} variants")

        # GTF annotation
        if cfg.gtf_file:
            gtf_csv = wd / "2_gtf_annotated.csv"
            if not gtf_csv.exists():
                logger.info(f"Step 3: GTF annotation")
                stage1_annotate_by_gtf(str(gt_csv), str(cfg.gtf_file), str(gtf_csv))
            df = pd.read_csv(str(gtf_csv), dtype=str)

        # ClinVar annotation
        if cfg.clinvar_vcf:
            clinvar_csv = wd / "3_clinvar_annotated.csv"
            if not clinvar_csv.exists():
                logger.info(f"Step 4: ClinVar annotation")
                stage1_merge_filter_clinvar(str(gtf_csv), str(cfg.clinvar_vcf), str(clinvar_csv))
            df = pd.read_csv(str(clinvar_csv), dtype=str)

        self._stage1_df = df
        logger.info(f"Stage 1 完成: {len(df)} variants, {len(df.columns)} columns")
        return df

    # ── Stage 2 ──────────────────────────────────────────────────────────

    def stage2_eqtl_annotation(
        self,
        symptoms: str,
        tissue_dir: Optional[str] = None,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Stage 2: GTEx eQTL 高级注释。

        Parameters
        ----------
        symptoms : str
            患者症状描述（用于 LLM 筛选相关组织）
        tissue_dir : str, optional
            GTEx eQTL parquet 目录
        **kwargs : 其他参数传递给 stage2_eqtl_annotation

        Returns
        -------
        pd.DataFrame: Stage 1 基础上追加 eQTL 列
        """
        if tissue_dir is None:
            tissue_dir = self.cfg.gtex_tissue_dir

        df = stage2_eqtl_annotation(
            stage1_csv=str(self.cfg.work_dir / "3_clinvar_annotated.csv"),
            tissue_dir=str(tissue_dir),
            symptoms=symptoms,
            llm_model=self.cfg.llm_model,
            api_key=self.cfg.api_key,
            base_url=self.cfg.base_url,
            **kwargs,
        )
        self._stage2_df = df
        return df

    # ── Stage 3 ──────────────────────────────────────────────────────────

    def stage3_score_and_rank(
        self,
        symptoms: str,
        stage2_csv: Optional[str] = None,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Stage 3: LLM 排序。

        Parameters
        ----------
        symptoms : str
            患者症状描述
        stage2_csv : str, optional
            Stage 2 输出 CSV（不填则用 Stage 1 输出）

        Returns
        -------
        pd.DataFrame: top-K 排序后的 candidates
        """
        from seekrare.scoring import Stage3Scorer

        csv_path = stage2_csv or str(self.cfg.work_dir / "3_clinvar_annotated.csv")

        scorer = Stage3Scorer(
            csv_path=csv_path,
            symptoms=symptoms,
            llm_provider=self.cfg.llm_provider,
            llm_model=self.cfg.llm_model,
            api_key=self.cfg.api_key,
            base_url=self.cfg.base_url,
            **kwargs,
        )
        result = scorer.run()
        self._stage3_df = result
        return result

    # ── Stage 4 ──────────────────────────────────────────────────────────

    def stage4_genos_analysis(
        self,
        stage3_csv: Optional[str] = None,
        top_n: int = 10,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Stage 4: Genos 模型分析。

        Parameters
        ----------
        stage3_csv : str, optional
            Stage 3 输出 CSV（不填则用 Stage 3 实例结果）
        top_n : int
            取前 N 个位点
        **kwargs : 传递给 Genos pipeline 的参数

        Returns
        -------
        pd.DataFrame: Genos peak 验证结果
        """
        if stage3_csv is None and hasattr(self, "_stage3_df"):
            stage3_df = self._stage3_df
            if stage3_df is None:
                raise ValueError("Stage 3 has not been run yet")
            csv_path = str(self.cfg.work_dir / "stage3_ranked.csv")
            stage3_df.to_csv(csv_path, index=False)
        else:
            csv_path = stage3_csv or str(self.cfg.work_dir / "3_clinvar_annotated.csv")

        result = stage4_genos_analysis(
            stage3_csv=csv_path,
            genome_fa=str(self.cfg.genome_fa),
            model_path=str(self.cfg.genos_model_path),
            output_dir=str(self.cfg.work_dir / "genos_result"),
            top_n=top_n,
            **kwargs,
        )
        return result

    # ── Full pipeline ──────────────────────────────────────────────────

    def run(
        self,
        symptoms: str,
        skip_stage2: bool = True,
        skip_stage4: bool = True,
    ) -> pd.DataFrame:
        """
        运行完整流水线。

        Parameters
        ----------
        symptoms : str
            患者症状描述
        skip_stage2 : bool
            若 True（默认），跳过 Stage 2
        skip_stage4 : bool
            若 True（默认），跳过 Stage 4

        Returns
        -------
        pd.DataFrame: Stage 3 排序结果
        """
        self.stage1_preprocess()

        if not skip_stage2 and self.cfg.gtex_tissue_dir:
            self.stage2_eqtl_annotation(symptoms=symptoms)

        result = self.stage3_score_and_rank(symptoms=symptoms)

        if not skip_stage4 and self.cfg.genos_model_path:
            self.stage4_genos_analysis(top_n=10)

        return result
