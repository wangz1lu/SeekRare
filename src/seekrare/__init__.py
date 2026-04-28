"""
SeekRare — 三阶段罕见病诊断系统

Three-stage pipeline:
    Stage 1: VCF 家系预处理 → 基本注释（必须）
    Stage 2: 高级注释（可选：eQTL / AlphaFold3 / Genos）
    Stage 3: LLM 分析 → 双动态打分 → 排序
"""

from seekrare.pipeline import SeekRarePipeline, SeekRareConfig

__version__ = "0.1.0"
__all__ = ["SeekRarePipeline", "SeekRareConfig"]
