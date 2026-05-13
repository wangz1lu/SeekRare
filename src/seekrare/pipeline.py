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
    run_family_preprocess,
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
    dbsnp_vcf: Optional[Union[str, Path]] = None
    work_dir: Union[str, Path] = "seekrare_output"

    # Stage 2
    gtex_tissue_dir: Optional[Union[str, Path]] = None
    splicevardb_tsv: Optional[Union[str, Path]] = None

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
            dbsnp_vcf="/ref/dbsnp.vcf.gz",
        )
        p.stage1_preprocess()
        p.stage2_eqtl_annotation(symptoms="眼部病变")
        p.stage3_score_and_rank(symptoms="视网膜色素变性")
        p.stage4_genos_analysis(top_n=5)
    """

    def __init__(
        self,
        vcf_proband: Optional[str] = None,
        vcf_father: Optional[str] = None,
        vcf_mother: Optional[str] = None,
        ref_fasta: Optional[str] = None,
        gtf_file: Optional[str] = None,
        clinvar_vcf: Optional[str] = None,
        dbsnp_vcf: Optional[str] = None,
        work_dir: str = "seekrare_output",
        gtex_tissue_dir: Optional[str] = None,
        llm_provider: str = "openai",
        llm_model: str = "deepseek-v4-flash",
        api_key: Optional[str] = None,
        base_url: Optional[str] = None,
        genome_fa: Optional[str] = None,
        genos_model_path: Optional[str] = None,
        **kwargs,
    ):
        """
        支持两种调用方式：

        方式A — 关键字参数（推荐）:
            SeekRarePipeline(vcf_proband="child.vcf.gz", ref_fasta="/ref/GRCh38.fa", ...)

        方式B — dataclass/config字典:
            cfg = SeekRareConfig(vcf_proband="child.vcf.gz", ...)
            SeekRarePipeline(config=cfg)
        """
        # 如果 vcf_proband 不是 None，说明用的是方式A
        if vcf_proband is not None:
            config = SeekRareConfig(
                vcf_proband=vcf_proband,
                vcf_father=vcf_father,
                vcf_mother=vcf_mother,
                ref_fasta=ref_fasta,
                gtf_file=gtf_file,
                clinvar_vcf=clinvar_vcf,
                dbsnp_vcf=dbsnp_vcf,
                work_dir=work_dir,
                gtex_tissue_dir=gtex_tissue_dir,
                llm_provider=llm_provider,
                llm_model=llm_model,
                api_key=api_key,
                base_url=base_url,
                genome_fa=genome_fa,
                genos_model_path=genos_model_path,
            )
        else:
            # 方式B: config 字典
            cfg_dict = kwargs.get("config", {})
            if isinstance(cfg_dict, dict):
                config = SeekRareConfig(**cfg_dict)
            else:
                config = cfg_dict

        self.cfg = config
        self._stage1_df: Optional[pd.DataFrame] = None
        self._llm_interpretation: Optional[dict] = None

    # ── Stage 1 ──────────────────────────────────────────────────────────

    def stage1_preprocess(self) -> pd.DataFrame:
        """
        Stage 1: VCF → 基本注释 CSV。

        自动检测家系模式：
          - 若提供了 father_vcf + mother_vcf → 家系 trio 模式
            （de_novo / recessive / xlinked 分类，基于 bcftools merge）
          - 否则 → 单样本模式

        Returns
        -------
        pd.DataFrame: 基本注释完成的 variants（含 inheritance_mode 列）
        """
        cfg = self.cfg

        result = run_family_preprocess(
            work_dir=str(cfg.work_dir),
            child_vcf=str(cfg.vcf_proband),
            father_vcf=str(cfg.vcf_father) if cfg.vcf_father else None,
            mother_vcf=str(cfg.vcf_mother) if cfg.vcf_mother else None,
            ref_fasta=str(cfg.ref_fasta) if cfg.ref_fasta else None,
            gtf_file=str(cfg.gtf_file) if cfg.gtf_file else None,
            clinvar_vcf=str(cfg.clinvar_vcf) if cfg.clinvar_vcf else None,
            dbsnp_vcf=str(cfg.dbsnp_vcf) if cfg.dbsnp_vcf else None,
        )

        self._stage1_df = result["combined"]
        self._stage1_modes = result["modes"]
        return self._stage1_df

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

    def stage2_splicevardb_annotation(
        self,
        stage1_csv: Optional[str] = None,
        splicevardb_tsv: Optional[str] = None,
        output_csv: Optional[str] = None,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Stage 2: SpliceVARDB 剪接变异注释。

        Parameters
        ----------
        stage1_csv : str, optional   Stage 1 输出 CSV
        splicevardb_tsv : str, optional   SpliceVARDB TSV 文件
        output_csv : str, optional   输出 CSV（默认覆盖 Stage 1 CSV）
        """
        if stage1_csv is None:
            stage1_csv = str(self.cfg.work_dir / "3_clinvar_annotated.csv")
        if splicevardb_tsv is None:
            splicevardb_tsv = str(self.cfg.splicevardb_tsv) if self.cfg.splicevardb_tsv else ""
        if not splicevardb_tsv:
            logger.warning("[SpliceVARDB] splicevardb_tsv not provided, skipping")
            return pd.read_csv(stage1_csv) if Path(stage1_csv).exists() else pd.DataFrame()
        if output_csv is None:
            output_csv = stage1_csv  # overwrite by default

        df = stage2_splicevardb_annotation(
            stage1_csv=stage1_csv,
            splicevardb_tsv=splicevardb_tsv,
            output_csv=output_csv,
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


# ─────────────────────────────────────────────────────────────────────────────
# Standalone Stage 3 function
# ─────────────────────────────────────────────────────────────────────────────

def stage3_score_and_rank(
    csv_path: str,
    symptoms: str,
    llm_provider: str = "openai",
    llm_model: str = "deepseek-v4-flash",
    api_key: Optional[str] = None,
    base_url: Optional[str] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Stage 3: LLM 双重动态评分排序（独立函数）。

    Parameters
    ----------
    csv_path : str
        输入 CSV 路径（Stage 1 或 Stage 2 输出）
    symptoms : str
        患者症状描述
    llm_provider : str
    llm_model : str
    api_key : str
    base_url : str
    **kwargs : 其他参数传递给 Stage3Scorer

    Returns
    -------
    pd.DataFrame: 排序后的 candidates
    """
    from seekrare.scoring import Stage3Scorer

    scorer = Stage3Scorer(
        csv_path=csv_path,
        symptoms=symptoms,
        llm_provider=llm_provider,
        llm_model=llm_model,
        api_key=api_key,
        base_url=base_url,
        **kwargs,
    )
    return scorer.run()


# ─────────────────────────────────────────────────────────────────────────────
# Standalone Stage 4 functions
# ─────────────────────────────────────────────────────────────────────────────

def stage4_genos_analysis(
    sites: Union[str, list],
    stage3_csv: Optional[str] = None,
    genome_fa: Optional[str] = None,
    model_path: Optional[str] = None,
    output_dir: str = "seekrare_output/genos_result",
    **kwargs,
) -> pd.DataFrame:
    """
    Stage 4A: Genos 模型分析（独立函数）。

    Parameters
    ----------
    sites : str or list[tuple]
        位点选择："top:N" / "rows:R1,R2-R3" / [(chrom,pos,REF,ALT), ...]
    stage3_csv : str, optional
        Stage 3 CSV（top/rows 模式需要）
    genome_fa : str
        参考基因组 FASTA
    model_path : str
        Genos 模型路径
    output_dir : str
        输出目录
    **kwargs : 其他参数

    Returns
    -------
    pd.DataFrame: peak 验证结果
    """
    from seekrare.preprocess.stage4_genos import run_genos_analysis

    return run_genos_analysis(
        sites=sites,
        stage3_csv=stage3_csv,
        genome_fa=genome_fa,
        model_path=model_path,
        output_dir=output_dir,
        **kwargs,
    )


def stage4_alphafold_prediction(
    csv_path: str,
    ref_fasta: str,
    gtf_file: str,
    output_dir: str = "seekrare_output/alphafold_results",
    top_n: int = 10,
    **kwargs,
) -> pd.DataFrame:
    """
    Stage 4B: AlphaFold3 蛋白结构预测（独立函数）。

    Parameters
    ----------
    csv_path : str
        Stage 3 CSV 路径
    ref_fasta : str
        参考基因组 FASTA
    gtf_file : str
        GTF 注释文件
    output_dir : str
        输出目录
    top_n : int
        取前 N 个位点
    **kwargs : 其他参数

    Returns
    -------
    pd.DataFrame: AlphaFold 预测结果汇总
    """
    from seekrare.preprocess.stage4_alphafold import run_alphafold_prediction

    return run_alphafold_prediction(
        stage3_csv=csv_path,
        ref_fasta=ref_fasta,
        gtf_file=gtf_file,
        output_dir=output_dir,
        top_n=top_n,
        **kwargs,
    )
