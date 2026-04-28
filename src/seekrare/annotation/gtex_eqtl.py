"""
GTEx eQTL annotation — Stage 2 optional annotation module.

对变异进行 GTEx 组织特异性 eQTL 注释。
使用 Parquet 格式的 GTEx eQTL 数据（每组织一个 .parquet 文件）。

GTEx 列说明:
- variant_id: chr_POS_REF_ALT 格式
- phenotype_id: 基因 symbol
- slope: eQTL effect size
- pval_nominal: eQTL p-value
- tss_distance: TSS distance
- ma_samples: samples with alt allele
- ma_count: count of alt allele

Usage:
    from seekrare.annotation import GTExEQTLAnnotator

    annotator = GTExEQTLAnnotator(
        tissue_dir="/path/to/GTEx_eQTL_parquet/",
    )
    df = annotator.annotate_variants(df)
    # 或输出为长格式 TSV（原始风格）
    annotator.annotate_to_tsv(vcf_path, output_path)
"""

from __future__ import annotations

import gzip
import os
import glob
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


class GTExEQTLAnnotator:
    """
    GTEx eQTL annotation.

    Parameters
    ----------
    tissue_dir : str or Path
        包含 GTEx .eQTLs.signif_pairs.parquet 文件的目录
        每个组织一个 parquet 文件
    min_pval : float
        eQTL p-value 阈值 (default: 1e-6)
    """

    def __init__(
        self,
        tissue_dir: Union[str, Path],
        min_pval: float = 1e-6,
    ):
        self.tissue_dir = Path(tissue_dir)
        self.min_pval = min_pval
        self._tissue_files: list[Path] = []

    def discover_tissues(self) -> list[str]:
        """发现所有可用组织."""
        self._tissue_files = sorted(self.tissue_dir.glob("*.eQTLs.signif_pairs.parquet"))
        tissues = []
        for f in self._tissue_files:
            # 文件名格式: Brain_Amygdala.v11.signif_pairs.parquet
            name = f.name.split(".v11")[0]
            tissues.append(name)
        logger.info(f"发现 {len(tissues)} 个 GTEx 组织: {tissues}")
        return tissues

    def _detect_suffix(self) -> str:
        """从第一个 parquet 文件检测 variant_id 后缀."""
        if not self._tissue_files:
            self.discover_tissues()

        f = self._tissue_files[0]
        gtex_sample = pd.read_parquet(f, columns=["variant_id"]).iloc[0]["variant_id"]
        # e.g. "chr1_12345_A_G_b38" → "_b38"
        parts = gtex_sample.rsplit("_", 1)
        suffix = "_" + parts[-1] if len(parts) > 1 else ""
        logger.info(f"GTEx variant_id suffix detected: '{suffix}'")
        return suffix

    def _build_vcf_set(self, vcf_path: str) -> set:
        """从 VCF 构建 variant_id 集合 (哈希查找，极快省内存)."""
        vcf_set = set()
        opener = gzip.open if vcf_path.endswith(".gz") else open

        with opener(vcf_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                chrom = cols[0].lstrip("chr")
                pos = cols[1]
                ref = cols[2].upper()
                alt = cols[3].upper()
                vid = f"chr{chrom}_{pos}_{ref}_{alt}"
                vcf_set.add(vid)

        logger.info(f"VCF 中共有 {len(vcf_set)} 个唯一变异位点")
        return vcf_set

    def annotate_to_tsv(
        self,
        vcf_path: str,
        output_path: str,
        tissue_dir: Optional[str] = None,
    ) -> str:
        """
        将 eQTL 注释结果写出为长格式 TSV（原始脚本风格）。

        参数
        ----
        vcf_path : str
            输入 VCF (.gz 或 plain)
        output_path : str
            输出 TSV 路径
        tissue_dir : str, optional
            覆盖 self.tissue_dir

        Returns
        -------
        str
            输出文件路径
        """
        if tissue_dir:
            self.tissue_dir = Path(tissue_dir)

        suffix = self._detect_suffix()
        vcf_set = self._build_vcf_set(vcf_path)

        tissue_files = sorted(Path(self.tissue_dir).glob("*.eQTLs.signif_pairs.parquet"))

        first_write = True
        total_matched = 0

        for f in tissue_files:
            t_name = f.name.split(".v11")[0]
            logger.info(f"  处理组织: {t_name}")

            gtex_chunk = pd.read_parquet(
                f,
                columns=["variant_id", "phenotype_id", "slope", "pval_nominal", "tss_distance"],
            )

            # 过滤：在 VCF 中 且 通过 p-value 阈值
            matched = gtex_chunk[
                gtex_chunk["variant_id"].isin(vcf_set)
                & (gtex_chunk["pval_nominal"] < self.min_pval)
            ].copy()
            matched["tissue"] = t_name

            if not matched.empty:
                mode = "w" if first_write else "a"
                header = first_write
                matched.to_csv(output_path, sep="\t", index=False, mode=mode, header=header)
                first_write = False
                total_matched += len(matched)
                logger.info(f"    匹配到 {len(matched)} 条记录")

            del gtex_chunk, matched

        logger.info(f"GTEx eQTL 注释完成: {total_matched} 条 → {output_path}")
        return output_path

    def annotate_variants(
        self,
        df: pd.DataFrame,
        tissues: Optional[list[str]] = None,
        top_n: int = 3,
    ) -> pd.DataFrame:
        """
        将 eQTL 信息作为新列添加到 variants DataFrame。

        对于每个变异，找其最显著的 eQTL（按 tissue 分组）。

        参数
        ----
        df : pd.DataFrame
            变异 DataFrame（需有 CHROM, POS, REF, ALT 列）
        tissues : list[str], optional
            指定要使用的组织列表（None = 全部）
        top_n : int
            每个变异保留的最显著 eQTL 数量 (default: 3)

        返回
        ----
        pd.DataFrame
            添加了 eQTL 列的原始 DataFrame
        """
        if not self._tissue_files:
            self.discover_tissues()

        suffix = self._detect_suffix()
        df = df.copy()

        # 构建 variant_id 列
        chrom_col = "CHROM" if "CHROM" in df.columns else "chrom"
        df["_variant_id"] = (
            "chr"
            + df[chrom_col].astype(str).str.lstrip("chr")
            + "_"
            + df["POS"].astype(str)
            + "_"
            + df["REF"].str.upper()
            + "_"
            + df["ALT"].str.upper()
            + suffix
        )

        # 收集所有组织的 eQTL 数据
        tissue_files = sorted(Path(self.tissue_dir).glob("*.eQTLs.signif_pairs.parquet"))
        if tissues:
            tissue_files = [f for f in tissue_files if f.name.split(".v11")[0] in tissues]

        all_eqtls = []
        for f in tissue_files:
            t_name = f.name.split(".v11")[0]
            chunk = pd.read_parquet(
                f,
                columns=["variant_id", "phenotype_id", "slope", "pval_nominal"],
            )
            chunk = chunk[chunk["pval_nominal"] < self.min_pval].copy()
            chunk["tissue"] = t_name
            all_eqtls.append(chunk)

        if not all_eqtls:
            logger.warning("未找到任何 eQTL 数据")
            df.drop(columns=["_variant_id"], inplace=True, errors="ignore")
            return df

        all_eqtls_df = pd.concat(all_eqtls, ignore_index=True)

        # 找出每个 variant 最显著的 eQTL
        all_eqtls_df = all_eqtls_df.sort_values("pval_nominal")
        best_per_variant = (
            all_eqtls_df.groupby("variant_id")
            .head(top_n)
            .groupby("variant_id")
            .agg(
                eqtl_gene=("phenotype_id", "first"),
                eqtl_pval=("pval_nominal", "first"),
                eqtl_slope=("slope", "first"),
                eqtl_tissue=("tissue", "first"),
                n_eqtl_tissues=("tissue", "nunique"),
            )
        ).reset_index()

        # 合并到 df
        n_before = len(df)
        df = df.merge(
            best_per_variant,
            left_on="_variant_id",
            right_on="variant_id",
            how="left",
        )
        df.drop(columns=["_variant_id", "variant_id"], inplace=True, errors="ignore")

        n_eqtl = df["eqtl_gene"].notna().sum()
        logger.info(f"GTEx eQTL: {n_eqtl}/{n_before} 变异有 eQTL 注释")

        return df
