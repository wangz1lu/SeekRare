"""
Variant annotation modules.

Stage 1 (必须注释):
    clinvar_loader.py   — ClinVar VCF/CSV 注释（二分距离查找）
    hpo_matcher.py      — HPO ontology 语义匹配（症状→HPO terms）
    combiner.py         — 多源注释合并（基础工具）

Stage 2 (可选注释，预留扩展接口):
    gtex_eqtl.py        — GTEx eQTL 组织特异性注释
    alphafold3.py       — AlphaFold3 Server / ColabFold 结构预测
    genos_annotation.py — Genos 临床遗传学模型注释（STUB）
"""

from seekrare.annotation.hpo_matcher import HPOMatcher, symptom_to_hpo
from seekrare.annotation.clinvar_loader import ClinVarLoader
from seekrare.annotation.gtex_eqtl import GTExEQTLAnnotator
from seekrare.annotation.alphafold3 import AlphaFold3Annotator, GenosAnnotationStub

__all__ = [
    # Stage 1 — required
    "HPOMatcher", "symptom_to_hpo",
    "ClinVarLoader",
    # Stage 2 — optional
    "GTExEQTLAnnotator",
    "AlphaFold3Annotator",
    "GenosAnnotationStub",
]
