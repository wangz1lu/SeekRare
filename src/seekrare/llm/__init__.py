"""
LLM integrations — Stage 3 modules.

symptom_parser.py  — 症状解析：文本 → HPO terms + weight vector
genos_client.py   — Genos 临床遗传学 LLM (STUB/预留)
alphafold_client.py — AlphaFold2 结构预测 (Stage 2 扩展)
"""

from seekrare.llm.symptom_parser import LLMSymptomParser

__all__ = ["LLMSymptomParser"]
