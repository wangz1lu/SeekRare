"""
SeekRare — 三阶段罕见病诊断系统

Stage 1: VCF 家系预处理 → 基本注释（必须）
Stage 2: 高级注释（可选：eQTL / AlphaFold3 / Genos）
Stage 3: LLM 分析 → 双动态打分 → 排序

快速开始:
    from seekrare import SeekRarePipeline
    pipeline = SeekRarePipeline(vcf_proband="child.vcf.gz", ...)
    result = pipeline.run(symptoms="智力障碍，癫痫...")
"""

from seekrare.pipeline import SeekRarePipeline, SeekRareConfig

__version__ = "0.1.0"
__all__ = ["SeekRarePipeline", "SeekRareConfig"]
