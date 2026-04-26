"""
SeekRare — LLM-powered rare disease diagnosis system.

Two-stage pipeline:
  Stage 1: VCF → Multi-source annotation (VEP, ClinVar, CADD, SpliceAI, OMIM)
  Stage 2: LLM interpretation → Dual-dynamic scoring → Model analysis
"""

from seekrare.pipeline import SeekRarePipeline, SeekRareConfig

__version__ = "0.1.0"
__all__ = ["SeekRarePipeline", "SeekRareConfig"]
