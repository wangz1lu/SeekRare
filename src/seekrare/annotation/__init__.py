"""
Variant annotation modules.

Stage 1 (基本注释 — Stage 1 流程中使用):
    clinvar_loader.py   — ClinVar VCF/CSV 注释（二分距离查找）
    hpo_matcher.py      — HPO ontology 语义匹配（Stage 3 打分中使用）

Stage 2 (高级注释 — Stage 2 流程中使用，需配置):
    gtex_eqtl.py        — GTEx eQTL 组织特异性注释
    alphafold3.py       — AlphaFold3 Server / ColabFold 结构预测
    genos_annotation.py — Genos 临床遗传学模型注释（STUB/预留）
"""

from seekrare.annotation.hpo_matcher import HPOMatcher, symptom_to_hpo
from seekrare.annotation.clinvar_loader import ClinVarLoader
from seekrare.annotation.gtex_eqtl import GTExEQTLAnnotator
from seekrare.annotation.alphafold3 import AlphaFold3Annotator, GenosAnnotationStub

__all__ = [
    "HPOMatcher", "symptom_to_hpo",
    "ClinVarLoader",
    "GTExEQTLAnnotator",
    "AlphaFold3Annotator",
    "GenosAnnotationStub",
]
