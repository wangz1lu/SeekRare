"""
LLM integrations — Stage 3 modules.

Stage 3: LLM-Powered Analysis
    symptom_parser.py  — 症状解析：症状文本 → HPO terms + weight vector
    genos_client.py   — Genos 临床遗传学专用 LLM (STUB/预留)
    alphafold_client.py — AlphaFold2 结构预测 (Stage 2 扩展)
"""

from seekrare.llm.symptom_parser import LLMSymptomParser

__all__ = ["LLMSymptomParser"]
