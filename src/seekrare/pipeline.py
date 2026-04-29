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
    输出: candidate_variants.csv

Usage:
    from seekrare import SeekRarePipeline

    pipeline = SeekRarePipeline(
        vcf_proband="child.vcf.gz",
        vcf_father="father.vcf.gz",
        vcf_mother="mother.vcf.gz",
        # Stage 1 — 四个文件路径（用户自行下载）
        ref_fasta="/path/to/GRCh38_no_alt_analysis_set.fa",
        gtf_file="/path/to/genomic.gtf",
        clinvar_vcf="/path/to/clinvar.vcf.gz",
        dbSNP_vcf="/path/to/dbsnp.vcf.gz",
        # Stage 2 (optional)
        gtex_tissue_dir="/path/to/GTEx_eQTL_parquet/",
    )

    # Stage 1 only（最快验证）
    df = pipeline.stage1_preprocess()
    df.to_csv("variants_annotated.csv", index=False)

    # Full pipeline
    result = pipeline.run(
        symptoms="智力障碍，癫痫，全身肌张力低",
    )
    result.to_csv("candidates.csv", index=False)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger

from seekrare.preprocess import (
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
# Configuration — 文件路径全部由用户配置，无默认值
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class SeekRareConfig:
    """
    SeekRare 全局配置。

    所有文件路径均为必填（无默认值），由用户配置。
    Stage 1 的四个文件下载地址见:
        ref_fasta   : https://ftp.ebi.ac.uk/pub/genomes/GRCh38/
        gtf_file    : https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
        clinvar_vcf : https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/
        dbSNP_vcf  : https://ftp.ncbi.nlm.nih.gov/species/newise/2024/
    """

    # ── Stage 1 — VCF & 注释文件（用户自行下载）────────────────────────────
    vcf_proband: Union[str, Path]           # 先证者 VCF (.vcf.gz)
    vcf_father: Optional[Union[str, Path]] = None   # 父亲 VCF
    vcf_mother: Optional[Union[str, Path]] = None   # 母亲 VCF

    # 注释资源文件（用户下载后填入路径）
    ref_fasta: Optional[Union[str, Path]] = None  # GRCh38 参考基因组 FASTA + .fai
    gtf_file: Optional[Union[str, Path]] = None   # NCBI genomic.gtf
    clinvar_vcf: Optional[Union[str, Path]] = None # ClinVar VCF (.vcf.gz)
    dbSNP_vcf: Optional[Union[str, Path]] = None  # dbSNP VCF (.vcf.gz)

    # 家系过滤模式
    inheritance_modes: tuple[str, ...] = ("denovo", "recessive")

    # ── Stage 2 — 高级注释（可选）────────────────────────────────────────────
    gtex_tissue_dir: Optional[Union[str, Path]] = None   # GTEx eQTL parquet 目录
    alphafold3_mode: Optional[str] = None                # "server" 或 "colabfold"
    alphafold3_url: Optional[str] = None
    alphafold3_api_key: Optional[str] = None
    genos_api_key: Optional[str] = None
    genos_url: Optional[str] = None

    # ── Stage 3 — LLM ───────────────────────────────────────────────────────
    llm_provider: str = "openai"
    llm_model: str = "gpt-4o"
    api_key: Optional[str] = None
    base_url: Optional[str] = None
    llm_temperature: float = 0.0
    top_k: int = 50

    # ── 全局 ───────────────────────────────────────────────────────────────
    work_dir: Union[str, Path] = "seekrare_output"

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

    def __init__(self, config: Union[SeekRareConfig, dict, None] = None, /, **kwargs):
        """
        支持两种初始化方式:

        1. 关键字参数直接传:
               SeekRarePipeline(vcf_proband="child.vcf.gz", gtf_file="...", ...)
        2. config dict 传参:
               SeekRarePipeline({"vcf_proband": "child.vcf.gz", ...})
        """
        if config is not None and kwargs:
            raise TypeError("不能同时传入 config 和关键字参数")
        if config is not None:
            if isinstance(config, dict):
                cfg_dict = config
            else:
                cfg_dict = vars(config)
        else:
            cfg_dict = kwargs
        cfg = SeekRareConfig(**cfg_dict)
        self.cfg = cfg
        self._llm: Optional[LLMSymptomParser] = None
        self._llm_interpretation: Optional[dict] = None
        self._stage1_df: Optional[pd.DataFrame] = None
        self._stage2_df: Optional[pd.DataFrame] = None

    # ─────────────────────────────────────────────────────────────────────────
    # Stage 1: VCF Preprocessing + Basic Annotation
    # ─────────────────────────────────────────────────────────────────────────

    def stage1_preprocess(self) -> pd.DataFrame:
        """
        Stage 1: VCF 家系预处理 + 基本注释。

        流程:
        1. bcftools: norm → merge → filter → inheritance 分类 → 剔除 dbSNP common
        2. compound_het_filter: Python 精细复合杂合过滤
        3. VCF → GT CSV (CHROM/POS/REF/ALT/GT per sample)
        4. GTF gene annotation (gene_name, feature_type)
        5. ClinVar + OMIM + HPO annotation (clinvar_sig, clinvar_mc, min_distance, omim_*)

        Returns
        -------
        pd.DataFrame
            基本注释完成的 variants DataFrame。
            基础列: CHROM, POS, REF, ALT, [GT_proband, GT_father, GT_mother]
            注释列: inheritance_type, gene_name, feature_type,
                   clinvar_sig, clinvar_mc, min_distance,
                   omim_diseases, omim_inheritance, omim_mim_number
        """
        logger.info("=" * 60)
        logger.info("Stage 1: VCF Preprocessing & Basic Annotation")
        logger.info("=" * 60)

        cfg = self.cfg
        wd = cfg.work_dir

        # ── 1. 家系 bcftools 预处理 ─────────────────────────────────────────
        if cfg.vcf_father and cfg.vcf_mother and cfg.ref_fasta:
            _assert_file(cfg.vcf_father, "父亲 VCF")
            _assert_file(cfg.vcf_mother, "母亲 VCF")
            _assert_file(cfg.ref_fasta, "参考基因组 FASTA")

            logger.info("Step 1: bcftools 家系过滤")
            bcf_result = run_bcftools_preprocess(
                outdir=str(wd / "bcftools"),
                father_vcf=str(cfg.vcf_father),
                mother_vcf=str(cfg.vcf_mother),
                child_vcf=str(cfg.vcf_proband),
                ref_fasta=str(cfg.ref_fasta),
                dbsnp_vcf=str(cfg.dbSNP_vcf) if cfg.dbSNP_vcf else None,
            )
            inheritance_vcf = self._select_inheritance_vcf(bcf_result)
            logger.info(f"  使用 inheritance VCF: {inheritance_vcf}")
        else:
            inheritance_vcf = str(cfg.vcf_proband)
            bcf_result = {}
            logger.info("Step 1: 无家系 VCF / ref_fasta，直接使用先证者 VCF")

        # ── 2. compound het 精细过滤 ────────────────────────────────────────
        if "compound_het" in cfg.inheritance_modes and bcf_result.get("father_het"):
            compound_csv = wd / "compound_het_candidates.csv"
            if not compound_csv.exists():
                logger.info("Step 2: compound het 精细过滤")
                run_compound_het_filter(
                    father_het_vcf=bcf_result["father_het"],
                    mother_het_vcf=bcf_result.get("mother_het", ""),
                    out_csv=str(compound_csv),
                    min_qual=20,
                )

        # ── 3. VCF → GT CSV ─────────────────────────────────────────────────
        _assert_file(cfg.vcf_proband, "先证者 VCF")
        gt_csv = wd / "1_gt.csv"
        if not gt_csv.exists():
            vcf_to_gt_csv(str(cfg.vcf_proband), str(gt_csv))
        df = pd.read_csv(gt_csv, dtype=str)
        df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
        logger.info(f"Step 3: VCF → GT CSV: {len(df)} variants")

        # ── 4. GTF gene annotation ───────────────────────────────────────────
        if cfg.gtf_file:
            _assert_file(cfg.gtf_file, "GTF 文件")
            gtf_csv = wd / "2_gtf_annotated.csv"
            if not gtf_csv.exists():
                annotate_by_gtf(str(gt_csv), str(cfg.gtf_file), str(gtf_csv))
                logger.info(f"Step 4: GTF annotation → {gtf_csv}")
            df = pd.read_csv(str(gtf_csv), dtype=str) if gtf_csv.exists() else df
            df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")

        # ── 5. ClinVar + OMIM + HPO annotation ───────────────────────────────
        if cfg.clinvar_vcf:
            _assert_file(cfg.clinvar_vcf, "ClinVar VCF")
            clinvar_csv = wd / "3_clinvar_annotated.csv"
            if not clinvar_csv.exists():
                merge_filter_clinvar(
                    str(gtf_csv),
                    str(cfg.clinvar_vcf),
                    str(clinvar_csv),
                )
                logger.info(f"Step 5: ClinVar+OMIM+HPO → {clinvar_csv}")
            df = pd.read_csv(str(clinvar_csv), dtype=str) if clinvar_csv.exists() else df
            df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")

        # ── 6. 标注 inheritance_type 列 ─────────────────────────────────────
        if bcf_result:
            df = self._annotate_inheritance_type(df, bcf_result)

        logger.info(f"Stage 1 完成: {len(df)} variants, {len(df.columns)} 列")
        logger.info(f"  列: {list(df.columns)}")
        self._stage1_df = df
        return df

    def _select_inheritance_vcf(self, bcf_result: dict) -> str:
        """根据 inheritance_modes 选择输入 VCF。"""
        modes = self.cfg.inheritance_modes
        for mode in ("denovo", "recessive", "compound_het", "xlinked"):
            if mode in modes:
                key = mode if mode != "compound_het" else "father_het"
                vcf = bcf_result.get(key, "")
                if vcf and Path(vcf).exists():
                    return vcf
        return str(self.cfg.vcf_proband)

    def _annotate_inheritance_type(
        self, df: pd.DataFrame, bcf_result: dict
    ) -> pd.DataFrame:
        """为每个变异标注其 inheritance type。"""
        df = df.copy()
        df["inheritance_type"] = "unknown"
        chrom_col = "CHROM" if "CHROM" in df.columns else "chrom"

        mode_vcfs = {
            "denovo": bcf_result.get("denovo"),
            "recessive": bcf_result.get("recessive"),
            "compound_het": bcf_result.get("father_het"),
            "xlinked": bcf_result.get("xlinked"),
        }

        for mode, vcf_path in mode_vcfs.items():
            if not vcf_path or not Path(vcf_path).exists():
                continue
            try:
                import gzip
                opener = gzip.open if str(vcf_path).endswith(".gz") else open
                with opener(vcf_path, "rt") as f:
                    for line in f:
                        if line.startswith("#"):
                            continue
                        cols = line.strip().split("\t")
                        chrom = cols[0].lstrip("chr")
                        pos = int(cols[1])
                        mask = (
                            (df[chrom_col].astype(str).str.lstrip("chr") == chrom)
                            & (df["POS"] == pos)
                        )
                        df.loc[mask, "inheritance_type"] = mode
            except Exception as e:
                logger.warning(f"Failed to annotate inheritance '{mode}': {e}")

        n = (df["inheritance_type"] != "unknown").sum()
        logger.info(f"  Inheritance annotated: {n}/{len(df)} variants")
        return df

    # ─────────────────────────────────────────────────────────────────────────
    # Stage 2: Advanced Annotation (Optional)
    # ─────────────────────────────────────────────────────────────────────────

    def stage2_advanced_annotation(
        self, df: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """
        Stage 2: 高级可选注释。

        - GTEx eQTL (gtex_tissue_dir)
        - AlphaFold3 (alphafold3_mode)
        - Genos (genos_api_key)

        Parameters
        ----------
        df : pd.DataFrame, optional
            Stage 1 输出。若 None，则重新运行 Stage 1。

        Returns
        -------
        pd.DataFrame
            Stage 1 基础上追加高级注释列。
        """
        logger.info("=" * 60)
        logger.info("Stage 2: Advanced Annotation (Optional)")
        logger.info("=" * 60)

        if df is None:
            df = getattr(self, "_stage1_df", None) or self.stage1_preprocess()

        cfg = self.cfg
        n_before = len(df.columns)

        if cfg.gtex_tissue_dir:
            logger.info(f"Step 2a: GTEx eQTL — {cfg.gtex_tissue_dir}")
            annotator = GTExEQTLAnnotator(tissue_dir=str(cfg.gtex_tissue_dir))
            df = annotator.annotate_variants(df)

        if cfg.alphafold3_mode:
            logger.info(f"Step 2b: AlphaFold3 (mode={cfg.alphafold3_mode})")
            annotator = AlphaFold3Annotator(
                mode=cfg.alphafold3_mode,
                base_url=cfg.alphafold3_url or "https://alphafold.ebi.ac.uk",
                api_key=cfg.alphafold3_api_key,
            )
            df = annotator.annotate_genes(df)

        if cfg.genos_api_key:
            logger.info("Step 2c: Genos annotation (STUB)")
            genos = GenosAnnotationStub(
                api_key=cfg.genos_api_key,
                base_url=cfg.genos_url or "https://api.genos.tech/v1",
            )
            df = genos.annotate_variants(df)

        if len(df.columns) == n_before:
            logger.info("  Stage 2: 无高级注释启用（未配置可选模块）")

        logger.info(f"Stage 2 完成: {len(df.columns) - n_before} 列新增，总 {len(df.columns)} 列")
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
            已知 HPO term IDs

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
                self.stage2_advanced_annotation()
                df = self._stage2_df

        # ── LLM symptom → HPO + weight ───────────────────────────────────────
        self._llm = LLMSymptomParser(
            provider=self.cfg.llm_provider,
            model=self.cfg.llm_model,
            api_key=self.cfg.api_key,
            base_url=self.cfg.base_url,
            temperature=self.cfg.llm_temperature,
        )
        llm_out = self._llm.interpret(symptoms)
        self._llm_interpretation = llm_out
        logger.info(f"  weight_vector: {llm_out.get('weight_vector', {})}")
        logger.info(f"  relevant_hpos: {[h['hpo_id'] for h in llm_out.get('relevant_hpos', [])]}")

        # ── HPO matching ─────────────────────────────────────────────────────
        relevant_hpos = llm_out.get("relevant_hpos", [])
        if hpo_terms:
            hpo_ids = {h["hpo_id"] for h in relevant_hpos}
            for hpo_id in hpo_terms:
                if hpo_id not in hpo_ids:
                    relevant_hpos.append({"hpo_id": hpo_id, "score": 0.9})

        # ── Dual-dynamic scoring ──────────────────────────────────────────────
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

        # ── Ranking ───────────────────────────────────────────────────────────
        df = rank_variants(df, top_k=self.cfg.top_k)

        logger.info(f"Stage 3 完成: {len(df)} candidates")
        return df

    # ─────────────────────────────────────────────────────────────────────────
    # Full pipeline
    # ─────────────────────────────────────────────────────────────────────────

    def run(
        self,
        symptoms: str,
        hpo_terms: Optional[list[str]] = None,
        skip_stage2: bool = True,
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
            若 True（默认），跳过 Stage 2，直接用 Stage 1 结果进 Stage 3。
            启用 Stage 2：传入 skip_stage2=False 并配置 gtex_tissue_dir 等。

        Returns
        -------
        pd.DataFrame
            最终排序后的候选 variants
        """
        self.stage1_preprocess()

        if not skip_stage2:
            self.stage2_advanced_annotation()

        df = self.stage3_analyze(symptoms=symptoms, hpo_terms=hpo_terms)

        out_path = self.cfg.work_dir / "final_candidates.csv"
        df.to_csv(out_path, index=False)
        logger.info(f"结果已保存: {out_path}")

        return df


# ─────────────────────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────────────────────

def _assert_file(path: Union[str, Path], label: str) -> Path:
    """检查文件是否存在，不存在则报错。"""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"{label} 文件不存在: {path}\n"
            f"请下载后配置正确路径。\n"
            f"下载地址参考: https://github.com/wangz1lu/SeekRare#data-resources"
        )
    return p
