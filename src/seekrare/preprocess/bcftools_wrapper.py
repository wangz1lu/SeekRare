"""
bcftools_wrapper.py — bcftools 家系过滤流水线的 Python 调用接口

提供 Python 级别的 API，可直接 import 使用，也可调用底层 bash 脚本。

Usage:
    from seekrare.preprocess import run_bcftools_preprocess

    result = run_bcftools_preprocess(
        outdir="/tmp/seekrare_preprocess",
        father_vcf="father.vcf.gz",
        mother_vcf="mother.vcf.gz",
        child_vcf="child.vcf.gz",
        ref_fasta="/ref/GRCh38.fa",
        dbsnp_vcf="/ref/dbsnp.vcf.gz",
        min_qual=30,
        min_dp=10,
        min_gq=20,
    )

    print(result["denovo"])       # path to de_novo VCF
    print(result["recessive"])    # path to recessive VCF
    print(result["father_het"])   # path to father het VCF (for compound het)
    print(result["mother_het"])   # path to mother het VCF (for compound het)
    print(result["xlinked"])      # path to X-linked VCF
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

from loguru import logger


@dataclass
class BcftoolsConfig:
    """
    bcftools 家系过滤配置。

    Attributes
    ----------
    outdir : str or Path
        输出目录
    father_vcf : str
        父亲 VCF (.vcf.gz)
    mother_vcf : str
        母亲 VCF (.vcf.gz)
    child_vcf : str
        子女/先证者 VCF (.vcf.gz)
    ref_fasta : str
        参考基因组 FASTA（GRCh38，带 .fai 索引）
    dbsnp_vcf : str, optional
        dbSNP VCF（用于剔除常见变异）
    min_qual : float
        最低 QUAL 阈值 (default: 30)
    min_dp : int
        最低测序深度 (default: 10)
    min_gq : int
        最低基因型质量 (default: 20)
    script_path : str
        bcftools_preprocess.sh 路径（默认用包内置）
    """
    outdir: Union[str, Path]
    father_vcf: str
    mother_vcf: str
    child_vcf: str
    ref_fasta: str
    dbsnp_vcf: Optional[str] = None
    min_qual: float = 30.0
    min_dp: int = 10
    min_gq: int = 20
    script_path: Optional[str] = None


def get_script_path() -> Path:
    """Get the bundled bcftools_preprocess.sh path."""
    import seekrare
    pkg_dir = Path(seekrare.__file__).parent
    script = pkg_dir.parent / "scripts" / "bcftools_preprocess.sh"

    # Fallback to source tree during dev
    if not script.exists():
        script = pkg_dir.parent.parent / "scripts" / "bcftools_preprocess.sh"

    if not script.exists():
        raise FileNotFoundError(
            f"bcftools_preprocess.sh not found at {script}. "
            "Set script_path explicitly or ensure package is properly installed."
        )
    return script


def run_bcftools_preprocess(
    outdir: Union[str, Path],
    father_vcf: str,
    mother_vcf: str,
    child_vcf: str,
    ref_fasta: str,
    dbsnp_vcf: Optional[str] = None,
    min_qual: float = 30.0,
    min_dp: int = 10,
    min_gq: int = 20,
    script_path: Optional[str] = None,
    dry_run: bool = False,
) -> dict:
    """
    Run bcftools trio preprocessing pipeline.

    Parameters
    ----------
    outdir : str or Path
        输出目录
    father_vcf, mother_vcf, child_vcf : str
        家系 VCF 文件路径
    ref_fasta : str
        参考基因组（GRCh38，带索引）
    dbsnp_vcf : str, optional
        dbSNP VCF（剔除常见变异）
    min_qual, min_dp, min_gq : float/int
        过滤阈值
    script_path : str, optional
        bcftools_preprocess.sh 路径
    dry_run : bool
        若 True，只打印命令而不执行

    Returns
    -------
    dict
        {
            "outdir": str,
            "denovo": str,           # de_novo.nocommon.vcf.gz
            "recessive": str,        # recessive.nocommon.vcf.gz
            "father_het": str,       # father_het.nocommon.vcf.gz
            "mother_het": str,       # mother_het.nocommon.vcf.gz
            "xlinked": str,          # xlinked.nocommon.vcf.gz
            "filtered_trio": str,    # 质控后的 trio VCF
        }
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Locate script
    if script_path is None:
        script = get_script_path()
    else:
        script = Path(script_path)

    if not script.exists():
        raise FileNotFoundError(f"Script not found: {script}")

    # Build env
    env = os.environ.copy()
    env["REF"] = str(Path(ref_fasta).resolve())
    env["MIN_QUAL"] = str(min_qual)
    env["MIN_DP"] = str(min_dp)
    env["MIN_GQ"] = str(min_gq)

    if dbsnp_vcf:
        env["DBSNP"] = str(Path(dbsnp_vcf).resolve())
    else:
        env["DBSNP"] = ""

    # Build command
    cmd = [
        "bash",
        str(script.resolve()),
        str(outdir.resolve()),
        str(Path(father_vcf).resolve()),
        str(Path(mother_vcf).resolve()),
        str(Path(child_vcf).resolve()),
    ]

    logger.info(f"Running bcftools preprocessing:")
    logger.info(f"  outdir: {outdir}")
    logger.info(f"  REF: {env['REF']}")
    logger.info(f"  DBSNP: {env.get('DBSNP', 'none')}")
    logger.info(f"  cmd: {' '.join(cmd)}")

    if dry_run:
        logger.info("[dry_run] Would execute bcftools pipeline")
        return _build_result_dict(outdir)

    # Run
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        env=env,
        timeout=3600,  # 1h max
    )

    if result.returncode != 0:
        logger.error(f"bcftools pipeline failed:\nSTDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}")
        raise RuntimeError(f"bcftools_preprocess.sh failed with code {result.returncode}")

    logger.info(f"bcftools pipeline completed:\n{result.stdout}")

    # Parse output to find file paths
    result_dict = _build_result_dict(outdir)

    # Save metadata
    meta_path = outdir / "preprocess_metadata.json"
    with open(meta_path, "w") as f:
        json.dump({
            "config": {
                "father_vcf": father_vcf,
                "mother_vcf": mother_vcf,
                "child_vcf": child_vcf,
                "ref_fasta": ref_fasta,
                "dbsnp_vcf": dbsnp_vcf,
                "min_qual": min_qual,
                "min_dp": min_dp,
                "min_gq": min_gq,
            },
            "outputs": result_dict,
        }, f, indent=2)

    logger.info(f"Metadata saved: {meta_path}")
    return result_dict


def _build_result_dict(outdir: Path) -> dict:
    """Build result dict with file paths."""
    final_dir = outdir / "final"
    inh_dir = outdir / "inheritance"
    filt_dir = outdir / "filtered"

    def path(subdir, name):
        p = outdir / subdir / name
        return str(p) if p.exists() else ""

    return {
        "outdir": str(outdir),
        "filtered_trio": path("filtered", "trio.filtered.vcf.gz"),
        "denovo": path("final", "denovo.nocommon.vcf.gz") or path("inheritance", "denovo.vcf.gz"),
        "recessive": path("final", "recessive.nocommon.vcf.gz") or path("inheritance", "recessive_hom.vcf.gz"),
        "father_het": path("final", "father_het.nocommon.vcf.gz") or path("inheritance", "father_het_only.vcf.gz"),
        "mother_het": path("final", "mother_het.nocommon.vcf.gz") or path("inheritance", "mother_het_only.vcf.gz"),
        "xlinked": path("final", "xlinked.nocommon.vcf.gz") or path("inheritance", "x_linked.vcf.gz"),
    }


def run_compound_het_filter(
    father_het_vcf: str,
    mother_het_vcf: str,
    out_csv: str,
    gtf_csv: Optional[str] = None,
    min_qual: int = 20,
    require_phase: bool = False,
    script_path: Optional[str] = None,
) -> str:
    """
    Run compound het Python filter on father_het and mother_het VCFs.

    Parameters
    ----------
    father_het_vcf, mother_het_vcf : str
        来自 run_bcftools_preprocess 的父亲/母亲 het VCF
    out_csv : str
        输出 CSV 路径
    gtf_csv : str, optional
        GTF 注释 CSV（用于补充基因信息）
    min_qual : int
        最低 QUAL 阈值
    require_phase : bool
        是否要求有相位
    script_path : str, optional
        compound_het_filter.py 路径

    Returns
    -------
    str
        输出 CSV 路径
    """
    if script_path is None:
        import seekrare
        pkg_dir = Path(seekrare.__file__).parent
        script = pkg_dir.parent / "scripts" / "compound_het_filter.py"
        if not script.exists():
            script = pkg_dir.parent.parent / "scripts" / "compound_het_filter.py"
    else:
        script = Path(script_path)

    cmd = [
        "python3", str(script.resolve()),
        father_het_vcf,
        mother_het_vcf,
        "--out", out_csv,
        "--min-qual", str(min_qual),
    ]
    if require_phase:
        cmd.append("--require-phase")
    if gtf_csv:
        cmd.extend(["--gtf", gtf_csv])

    logger.info(f"Running compound het filter: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        logger.error(f"compound_het_filter failed:\n{result.stderr}")
        raise RuntimeError(f"compound_het_filter.py failed: {result.stderr}")

    logger.info(f"Compound het filter output:\n{result.stdout}")
    return out_csv
