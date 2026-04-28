"""
SeekRare — 三阶段罕见病诊断系统

Stage 1 (必须): VCF 家系预处理 → 基本注释
    bcftools norm/merge/filter → inheritance 分类
    → VCF→GT CSV → GTF gene annotation → ClinVar+OMIM+HPO annotation
    输出: variants.csv (基础列 + 必须注释列)

Stage 2 (可选): 高级注释（按需启用）
    eQTL annotation (GTEx)
    AlphaFold3 structure prediction
    Genos model annotation (STUB)
    输出: 在 Stage 1 基础上追加高级注释列

Stage 3 (LLM 分析):
    LLM symptom → HPO + weight vector
    Dual-dynamic scoring
    Ranking
    输出: candidate_variants.csv + 可选 model analysis report

Usage:
    from seekrare import SeekRarePipeline

    pipeline = SeekRarePipeline(
        vcf_proband="child.vcf.gz",
        vcf_father="father.vcf.gz",
        vcf_mother="mother.vcf.gz",
        ref_fasta="/ref/GRCh38.fa",
        gtf_file="/ref/genomic.gtf",
        clinvar_vcf="/ref/clinvar.vcf.gz",
        # Stage 2 (optional):
        gtex_tissue_dir="/ref/GTEx_eQTL_parquet/",
        # alphafold3_server="https://alphafold.ebi.ac.uk",
        # genos_api_key="...",
    )

    # Stage 1 only
    df = pipeline.stage1_preprocess()
    df.to_csv("variants_annotated.csv", index=False)

    # Stage 1 + 2
    df = pipeline.stage1_preprocess()
    df = pipeline.stage2_advanced_annotation(df)

    # Full pipeline (Stage 1 + 2 + 3)
    result = pipeline.run(
        symptoms="智力障碍，癫痫，全身肌张力低",
    )
    result.to_csv("candidates.csv", index=False)
"""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger

from seekrare.preprocessing import (
    vcf_to_gt_csv,
    annotate_by_gtf,
    merge_filter_clinvar,
    run_bcftools_preprocess,
    run_compound_het_filter,
)
from seekrare.annotation import (
    ClinVarLoader,
    HPOMatcher,
    GTExEQTLAnnotator,
    AlphaFold3Annotator,
    GenosAnnotationStub,
)
from seekrare.llm import LLMSymptomParser
from seekrare.scoring import DualDynamicScorer, rank_variants


# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class SeekRareConfig:
    """
    SeekRare 全局配置。

    Stage 1 参数（必须）:
        vcf_proband  : 先证者 VCF
        gtf_file     : NCBI genomic.gtf 基因注释
        clinvar_vcf  : ClinVar VCF 致病性注释

    Stage 1 参数（可选）:
        vcf_father / vcf_mother : 家系 VCF（用于 inheritance 分类）
        ref_fasta                : GRCh38 参考基因组（bcftools 需要）
        dbSNP_vcf               : dbSNP VCF（剔除 common 变异）

    Stage 2 参数（可选高级注释）:
        gtex_tissue_dir  : GTEx eQTL parquet 目录
        alphafold3_mode : "server" 或 "colabfold"
        alphafold3_url   : AlphaFold3 API URL
        genos_api_key    : Genos API key
        genos_url        : Genos API base URL

    Stage 3 参数:
        llm_provider  : "openai" / "anthropic" / "local"
        llm_model     : 模型名
        api_key       : API key
        base_url      : OpenAI-compatible base URL
        top_k         : 返回 top-K 候选数
    """

    # ── Stage 1 ───────────────────────────────────────────────────────────────
    vcf_proband: Union[str, Path]
    vcf_father: Optional[Union[str, Path]] = None
    vcf_mother: Optional[Union[str, Path]] = None
    ref_fasta: Optional[Union[str, Path]] = None
    gtf_file: Optional[Union[str, Path]] = None
    clinvar_vcf: Optional[Union[str, Path]] = None
    dbSNP_vcf: Optional[Union[str, Path]] = None

    # ── Stage 2 ─────────────────────────────────────────────────────────────
    gtex_tissue_dir: Optional[Union[str, Path]] = None
    alphafold3_mode: Optional[str] = None       # "server" or "colabfold"
    alphafold3_url: Optional[str] = None
    alphafold3_api_key: Optional[str] = None
    genos_api_key: Optional[str] = None
    genos_url: Optional[str] = None

    # ── Stage 3 ─────────────────────────────────────────────────────────────
    llm_provider: str = "openai"
    llm_model: str = "gpt-4o"
    api_key: Optional[str] = None
    base_url: Optional[str] = None
    llm_temperature: float = 0.0
    top_k: int = 50

    # ── Global ───────────────────────────────────────────────────────────────
    work_dir: Union[str, Path] = "seekrare_output"
    inheritance_modes: tuple[str, ...] = ("denovo", "recessive")

    def __post_init__(self):
        self.work_dir = Path(self.work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# Pipeline
# ─────────────────────────────────────────────────────────────────────────────

class SeekRarePipeline:
    """
    三阶段罕见病诊断流水线。

    Stage 1: VCF → 基本注释 CSV（必须）
    Stage 2: 高级注释（可选，eQTL / AlphaFold3 / Genos）
    Stage 3: LLM 分析（打分 + 排序）
    """

    def __init__(self, config: Union[SeekRareConfig, dict]):
        if isinstance(config, dict):
            config = SeekRareConfig(**config)
        self.cfg = config
        self._llm: Optional[LLMSymptomParser] = None
        self._hpo_matcher: Optional[HPOMatcher] = None
        self._gtf_annotator: Optional[ClinVarLoader] = None

    # ── Lazy initializers ────────────────────────────────────────────────────

    def _init_llm(self):
        if self._llm is None:
            self._llm = LLMSymptomParser(
                provider=self.cfg.llm_provider,
                model=self.cfg.llm_model,
                api_key=self.cfg.api_key,
                base_url=self.cfg.base_url,
                temperature=self.cfg.llm_temperature,
            )

    # ─────────────────────────────────────────────────────────────────────────
    # Stage 1: VCF Preprocessing + Basic Annotation
    # ─────────────────────────────────────────────────────────────────────────

    def stage1_preprocess(self) -> pd.DataFrame:
        """
        Stage 1: VCF 家系预处理 + 基本注释。

        完整流程:
        1. bcftools: norm → merge trio → quality filter → inheritance 分类 → 剔除 dbSNP common
        2. compound_het_filter: Python 精细复合杂合过滤
        3. VCF → GT CSV: CHROM / POS / REF / ALT / GT per sample
        4. GTF gene annotation: gene_name, feature_type
        5. ClinVar + OMIM + HPO annotation

        Returns
        -------
        pd.DataFrame
            基本注释完成的 variants DataFrame。
            基础列: CHROM, POS, REF, ALT, [GT_proband, GT_father, GT_mother]
            注释列: inheritance_type, gene_name, feature_type,
                   clinvar_sig, clinvar_mc, min_distance,
                   omim_diseases, omim_inheritance, omim_mim_number,
                   hpo_terms (from gene-phenotype mapping)
        """
        logger.info("=" * 60)
        logger.info("Stage 1: VCF Preprocessing & Basic Annotation")
        logger.info("=" * 60)

        cfg = self.cfg
        wd = cfg.work_dir

        # ── 1a. 家系 bcftools 预处理 ────────────────────────────────────────
        if cfg.vcf_father and cfg.vcf_mother and cfg.ref_fasta:
            logger.info("Step 1: bcftools 家系过滤")
            bcf_result = run_bcftools_preprocess(
                outdir=str(wd / "bcftools"),
                father_vcf=str(cfg.vcf_father),
                mother_vcf=str(cfg.vcf_mother),
                child_vcf=str(cfg.vcf_proband),
                ref_fasta=str(cfg.ref_fasta),
                dbsnp_vcf=str(cfg.dbSNP_vcf) if cfg.dbSNP_vcf else None,
            )
            # 用 inheritance 过滤后的 VCF 继续下游
            # 根据 inheritance_modes 决定输入文件
            inheritance_vcf = self._select_inheritance_vcf(bcf_result)
            logger.info(f"  Using inheritance VCF: {inheritance_vcf}")
        else:
            # 无家系 VCF，直接用先证者 VCF
            inheritance_vcf = str(cfg.vcf_proband)
            logger.info(f"Step 1: No trio VCFs, using proband VCF directly")

        # ── 1b. compound het 精细过滤 ─────────────────────────────────────────
        compound_het_csv = wd / "compound_het_candidates.csv"
        if "compound_het" in cfg.inheritance_modes and cfg.vcf_father and cfg.vcf_mother:
            logger.info("Step 1b: compound het 精细过滤")
            if bcf_result.get("father_het") and bcf_result.get("mother_het"):
                run_compound_het_filter(
                    father_het_vcf=bcf_result["father_het"],
                    mother_het_vcf=bcf_result["mother_het"],
                    out_csv=str(compound_het_csv),
                    min_qual=20,
                )

        # ── 2. VCF → GT CSV ───────────────────────────────────────────────────
        gt_csv = wd / "1_gt.csv"
        if not gt_csv.exists():
            vcf_to_gt_csv(str(cfg.vcf_proband), str(gt_csv))
            logger.info(f"Step 2: VCF → GT CSV: {gt_csv}")
        else:
            logger.info(f"  [跳过] {gt_csv} 已存在")

        df = pd.read_csv(gt_csv, dtype=str)
        df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
        logger.info(f"  {len(df)} variants from GT CSV")

        # ── 3. GTF gene annotation ─────────────────────────────────────────────
        if cfg.gtf_file:
            gtf_csv = wd / "2_gtf_annotated.csv"
            if not gtf_csv.exists():
                annotate_by_gtf(str(gt_csv), str(cfg.gtf_file), str(gtf_csv))
                logger.info(f"Step 3: GTF annotation → {gtf_csv}")
            else:
                logger.info(f"  [跳过] {gtf_csv} 已存在")

            df = pd.read_csv(str(gtf_csv), dtype=str) if gtf_csv.exists() else df
            df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
            logger.info(f"  After GTF: {len(df)} variants")

        # ── 4. ClinVar + OMIM + HPO annotation ─────────────────────────────────
        if cfg.clinvar_vcf:
            clinvar_csv = wd / "3_clinvar_annotated.csv"
            if not clinvar_csv.exists():
                merge_filter_clinvar(
                    str(gtf_csv) if gtf_csv.exists() else str(gt_csv),
                    str(cfg.clinvar_vcf),
                    str(clinvar_csv),
                )
                logger.info(f"Step 4: ClinVar+OMIM+HPO → {clinvar_csv}")
            else:
                logger.info(f"  [跳过] {clinvar_csv} 已存在")

            df = pd.read_csv(str(clinvar_csv), dtype=str) if clinvar_csv.exists() else df
            df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
            logger.info(f"  After ClinVar: {len(df)} variants")

        # ── 5. 添加 inheritance_type 列 ─────────────────────────────────────────
        if cfg.vcf_father and cfg.vcf_mother:
            df = self._add_inheritance_column(df, bcf_result)

        # ── 6. 加载 compound het 结果 ──────────────────────────────────────────
        if compound_het_csv.exists():
            df_compound = pd.read_csv(str(compound_het_csv), dtype=str)
            n_compound = len(df_compound["gene"].unique())
            logger.info(f"  Compound het candidates: {n_compound} genes from {compound_het_csv}")

        logger.info(f"Stage 1 完成: {len(df)} variants, 列数={len(df.columns)}")
        logger.info(f"  列名: {list(df.columns)}")

        self._stage1_df = df
        return df

    def _select_inheritance_vcf(self, bcf_result: dict) -> str:
        """根据 inheritance_modes 选择输入 VCF。"""
        modes = self.cfg.inheritance_modes

        # 优先级: denovo > recessive > compound_het > xlinked
        for mode in ("denovo", "recessive", "compound_het", "xlinked"):
            if mode in modes:
                key = mode if mode != "compound_het" else "father_het"
                vcf = bcf_result.get(key, "")
                if vcf and Path(vcf).exists():
                    return vcf

        # Fallback: denovo
        vcf = bcf_result.get("denovo", "")
        if vcf and Path(vcf).exists():
            return vcf

        # 最后 fallback: proband VCF
        return str(self.cfg.vcf_proband)

    def _add_inheritance_column(
        self,
        df: pd.DataFrame,
        bcf_result: dict,
    ) -> pd.DataFrame:
        """为每个变异标注其 inheritance type。"""
        df = df.copy()
        df["inheritance_type"] = "unknown"

        mode_vcfs = {
            "denovo": bcf_result.get("denovo"),
            "recessive": bcf_result.get("recessive"),
            "father_het": bcf_result.get("father_het"),
            "mother_het": bcf_result.get("mother_het"),
            "xlinked": bcf_result.get("xlinked"),
        }

        chrom_col = "CHROM" if "CHROM" in df.columns else "chrom"

        for mode, vcf_path in mode_vcfs.items():
            if not vcf_path or not Path(vcf_path).exists():
                continue

            try:
                import gzip
                opener = gzip.open if vcf_path.endswith(".gz") else open
                with opener(vcf_path, "rt") as f:
                    for line in f:
                        if line.startswith("#"):
                            continue
                        cols = line.strip().split("\t")
                        chrom = cols[0].lstrip("chr")
                        pos = int(cols[1])
                        mask = (df[chrom_col].astype(str).str.lstrip("chr") == chrom) & (df["POS"] == pos)
                        df.loc[mask, "inheritance_type"] = mode
            except Exception as e:
                logger.warning(f"Failed to annotate inheritance type '{mode}': {e}")

        n_annotated = (df["inheritance_type"] != "unknown").sum()
        logger.info(f"  Inheritance annotated: {n_annotated}/{len(df)} variants")
        return df

    # ─────────────────────────────────────────────────────────────────────────
    # Stage 2: Advanced Annotation (Optional)
    # ─────────────────────────────────────────────────────────────────────────

    def stage2_advanced_annotation(self, df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        Stage 2: 高级可选注释。

        按配置依次调用:
        - eQTL (GTEx): gtex_tissue_dir
        - AlphaFold3: alphafold3_mode
        - Genos: genos_api_key

        Parameters
        ----------
        df : pd.DataFrame, optional
            Stage 1 输出。若 None，则重新运行 Stage 1。

        Returns
        -------
        pd.DataFrame
            Stage 1 基础上追加高级注释列。
            新增列: eqtl_gene, eqtl_pval, eqtl_tissue (eQTL)
                   alphafold_predicted, alphafold_pdb_url (AlphaFold3)
                   genos_pathogenicity, genos_acmg_criteria (Genos, STUB)
        """
        logger.info("=" * 60)
        logger.info("Stage 2: Advanced Annotation (Optional)")
        logger.info("=" * 60)

        if df is None:
            df = getattr(self, "_stage1_df", None)
            if df is None:
                df = self.stage1_preprocess()

        cfg = self.cfg
        n_before = len(df.columns)

        # ── 2a. GTEx eQTL ─────────────────────────────────────────────────────
        if cfg.gtex_tissue_dir and Path(cfg.gtex_tissue_dir).exists():
            logger.info("Step 2a: GTEx eQTL annotation")
            annotator = GTExEQTLAnnotator(tissue_dir=str(cfg.gtex_tissue_dir))
            df = annotator.annotate_variants(df)
            logger.info(f"  eQTL 列已添加: {len(df.columns) - n_before} 列")
        else:
            logger.info("Step 2a: GTEx eQTL 跳过 (gtex_tissue_dir 未配置)")

        # ── 2b. AlphaFold3 ────────────────────────────────────────────────────
        if cfg.alphafold3_mode:
            logger.info(f"Step 2b: AlphaFold3 annotation (mode={cfg.alphafold3_mode})")
            annotator = AlphaFold3Annotator(
                mode=cfg.alphafold3_mode,
                base_url=cfg.alphafold3_url or "https://alphafold.ebi.ac.uk",
                api_key=cfg.alphafold3_api_key,
            )
            df = annotator.annotate_genes(df)
            logger.info(f"  AlphaFold3 列已添加: {len(df.columns) - n_before} 列")
        else:
            logger.info("Step 2b: AlphaFold3 跳过 (未配置)")

        # ── 2c. Genos ───────────────────────────────────────────────────────────
        if cfg.genos_api_key:
            logger.info("Step 2c: Genos annotation (STUB)")
            genos = GenosAnnotationStub(
                api_key=cfg.genos_api_key,
                base_url=cfg.genos_url or "https://api.genos.tech/v1",
            )
            df = genos.annotate_variants(df)
            logger.info(f"  Genos 列已添加 (STUB)")
        else:
            logger.info("Step 2c: Genos 跳过 (未配置)")

        n_added = len(df.columns) - n_before
        logger.info(f"Stage 2 完成: 新增 {n_added} 列，总列数={len(df.columns)}")
        logger.info(f"  所有列: {list(df.columns)}")

        self._stage2_df = df
        return df

    # ─────────────────────────────────────────────────────────────────────────
    # Stage 3: LLM-Powered Analysis
    # ─────────────────────────────────────────────────────────────────────────

    def stage3_analyze(
        self,
        df: Optional[pd.DataFrame] = None,
        symptoms: str = "",
        hpo_terms: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        """
        Stage 3: LLM 分析 → 打分 → 排序。

        Parameters
        ----------
        df : pd.DataFrame, optional
            Stage 2 输出。若 None，则重新运行 Stage 1+2。
        symptoms : str
            患者表型描述（自由文本）
        hpo_terms : list[str], optional
            已知 HPO terms

        Returns
        -------
        pd.DataFrame
            打分 + 排序后的 top-K 候选 variants
        """
        logger.info("=" * 60)
        logger.info("Stage 3: LLM-Powered Analysis")
        logger.info("=" * 60)

        if df is None:
            df = getattr(self, "_stage2_df", None)
            if df is None:
                self.stage1_preprocess()
                df = self.stage2_advanced_annotation()

        # ── LLM symptom interpretation ─────────────────────────────────────────
        self._init_llm()
        llm_out = self._llm.interpret(symptoms)
        self._llm_interpretation = llm_out

        logger.info(f"  LLM weight_vector: {llm_out.get('weight_vector', {})}")
        logger.info(f"  LLM relevant_hpos: {[h['hpo_id'] for h in llm_out.get('relevant_hpos', [])]}")

        # ── HPO matching ────────────────────────────────────────────────────────
        relevant_hpos = llm_out.get("relevant_hpos", [])
        if hpo_terms:
            hpo_ids = {h["hpo_id"] for h in relevant_hpos}
            for hpo_id in hpo_terms:
                if hpo_id not in hpo_ids:
                    relevant_hpos.append({"hpo_id": hpo_id, "score": 0.9})

        # ── Dual-dynamic scoring ────────────────────────────────────────────────
        weight_vector = llm_out.get("weight_vector", {})
        if not weight_vector:
            weight_vector = {
                "clinvar_sig": 0.35,
                "impact": 0.20,
                "gnomad_af": 0.25,
                "hpo_terms": 0.20,
            }

        scorer = DualDynamicScorer(weight_vector)
        df = scorer.score(df, relevant_hpos=relevant_hpos)

        # ── Ranking ─────────────────────────────────────────────────────────────
        df = rank_variants(df, top_k=self.cfg.top_k)

        logger.info(f"Stage 3 完成: top {len(df)} candidates")
        return df

    # ─────────────────────────────────────────────────────────────────────────
    # Full pipeline
    # ─────────────────────────────────────────────────────────────────────────

    def run(
        self,
        symptoms: str,
        hpo_terms: Optional[list[str]] = None,
        skip_stage2: bool = False,
    ) -> pd.DataFrame:
        """
        运行完整三阶段流水线。

        Parameters
        ----------
        symptoms : str
            患者表型描述
        hpo_terms : list[str], optional
            已知 HPO terms
        skip_stage2 : bool
            若 True，跳过 Stage 2（直接用 Stage 1 结果进 Stage 3）

        Returns
        -------
        pd.DataFrame
            最终排序后的候选 variants
        """
        self.stage1_preprocess()

        if not skip_stage2:
            self.stage2_advanced_annotation()

        df = self.stage3_analyze(
            symptoms=symptoms,
            hpo_terms=hpo_terms,
        )

        # 保存结果
        out_path = self.cfg.work_dir / "final_candidates.csv"
        df.to_csv(out_path, index=False)
        logger.info(f"结果已保存: {out_path}")

        return df
