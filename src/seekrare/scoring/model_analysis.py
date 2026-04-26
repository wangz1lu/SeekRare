"""
Model Analysis — post-filtering LLM model analysis of candidate variants.

After SeekRare scoring + filtering, use an LLM to perform deep
pathogenicity analysis on the top candidates:

- Re-evaluate ACMG criteria
- Integrate all annotation sources (ClinVar, VEP, CADD, SpliceAI, OMIM)
- Assess genotype-phenotype consistency
- Rank candidates by 综合evidence
- Generate differential diagnosis list
- Suggest functional studies

Accepts any OpenAI-compatible API (OpenAI, Anthropic, local, Genos).

Usage:
    from seekrare.scoring import ModelAnalyzer
    analyzer = ModelAnalyzer(provider="openai", model="gpt-4o")
    report = analyzer.analyze(
        variants_df=filtered_df,
        patient_phenotype="...",
        hpo_terms=["HP:0001250", ...],
        top_k=10,
    )
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any, Optional, Union

import pandas as pd
from loguru import logger


@dataclass
class ModelAnalyzerConfig:
    """Configuration for model analysis."""
    provider: str = "openai"        # "openai", "anthropic", "local", "genos"
    model: str = "gpt-4o"
    api_key: Optional[str] = None
    base_url: Optional[str] = None   # For local/self-hosted
    temperature: float = 0.0
    max_tokens: int = 8192


SYSTEM_PROMPT = """You are a senior clinical molecular geneticist and genomic data scientist.
You analyze candidate variants from a rare disease NGS pipeline and provide
expert interpretation in a clinical genetics report.

Your expertise:
- ACMG/AMP variant classification guidelines
- Variant effect prediction (splice, missense, nonsense, frameshift)
- Population frequency interpretation (gnomAD, ExAC)
- ClinVar evidence synthesis
- Protein domain and structural impact (AlphaFold2)
- Splicing impact (SpliceAI, dbscSNV)
- Genotype-phenotype correlation
- OMIM/Orphanet disease matching
- Compound heterozygous and other complex inheritance models

You will receive:
1. Patient phenotype + HPO terms
2. A ranked list of candidate variants with full annotations
3. SeekRare scoring results

Your task:
1. Re-evaluate ACMG criteria for each candidate
2. Synthesize all evidence into a final pathogenicity classification
3. Assess how well each variant explains the patient's phenotype
4. Generate a ranked differential diagnosis
5. Recommend specific follow-up studies (segregation analysis, functional assays, etc.)
6. Flag any variants that need special caution (variants of uncertain significance
   that are likely benign, or pathogenic variants with no phenotypic fit)

Be precise, cite specific evidence, and distinguish definite conclusions
from reasonable hypotheses. Use Chinese/English bilingual for clinical terms."""

SYSTEM_PROMPT_CN = """你是一位资深临床分子遗传学家和基因组数据科学家。
你负责分析罕见病NGS流程中的候选变异，并提供专业的临床遗传学解读。

你的专长：
- ACMG/AMP 变异分类指南
- 变异效应预测（剪接、错义、无义、移码）
- 等位基因频率解读（gnomAD、ExAC）
- ClinVar 证据综合
- 蛋白结构域与结构影响（AlphaFold2）
- 剪接影响（SpliceAI、dbscSNV）
- 基因型-表型相关性
- OMIM/Orphanet 疾病匹配
- 复合杂合及其他复杂遗传模式

你会收到：
1. 患者表型 + HPO terms
2. 候选变异列表（含完整注释）
3. SeekRare 评分结果

你的任务：
1. 重新评估每个候选变异的 ACMG 标准
2. 综合所有证据给出最终致病性分类
3. 评估每个变异对患者表型的解释程度
4. 生成排序的鉴别诊断列表
5. 推荐具体的随访研究（分离分析、功能验证等）
6. 标注需要特别谨慎的变异

精确、引用具体证据，区分明确结论与合理假设。"""


class ModelAnalyzer:
    """
    Post-filtering LLM model analysis of candidate variants.

    Parameters
    ----------
    config : ModelAnalyzerConfig or dict
        Analyzer configuration
    """

    def __init__(self, config: Union[ModelAnalyzerConfig, dict]):
        if isinstance(config, dict):
            self.config = ModelAnalyzerConfig(**config)
        else:
            self.config = config

        self._client = None

    def _get_client(self):
        """Lazy-init LLM client based on provider."""
        if self._client is not None:
            return self._client

        provider = self.config.provider.lower()
        api_key = self.config.api_key or os.getenv("_API_KEY", "")
        base_url = self.config.base_url

        if provider in ("openai", "local", "genos"):
            try:
                from openai import OpenAI
            except ImportError:
                raise ImportError("openai package required: pip install openai")

            self._client = OpenAI(api_key=api_key, base_url=base_url)

        elif provider == "anthropic":
            try:
                import anthropic
            except ImportError:
                raise ImportError("anthropic package required: pip install anthropic")

            self._client = anthropic.Anthropic(api_key=api_key)

        else:
            raise ValueError(f"Unknown provider: {provider}")

        return self._client

    def _build_variant_summary(self, df: pd.DataFrame) -> str:
        """Build compact variant summary for LLM prompt."""
        cols = ["gene_name", "CHROM", "POS", "REF", "ALT",
                "clinvar_sig", "consequence", "cadd_phred",
                "gnomad_af", "impact", "sift_prediction",
                "polyphen_prediction", "spliceai_max_delta",
                "omim_diseases", "omim_inheritance",
                "seekrare_score"]

        available = [c for c in cols if c in df.columns]
        summary_lines = []

        for rank, (_, row) in enumerate(df.iterrows(), 1):
            parts = [f"#{rank} {row.get('gene_name', '?')} "]
            for col in available:
                val = row.get(col)
                if val is not None and str(val).strip():
                    parts.append(f"{col}={val}")
            summary_lines.append(" | ".join(parts))

        return "\n".join(summary_lines)

    def analyze(
        self,
        variants_df: pd.DataFrame,
        patient_phenotype: str,
        hpo_terms: Optional[list[str]] = None,
        top_k: int = 10,
        language: str = "en",
    ) -> dict:
        """
        Run model analysis on top candidate variants.

        Parameters
        ----------
        variants_df : pd.DataFrame
            Filtered/scored variants from SeekRare pipeline
        patient_phenotype : str
            Free-text patient phenotype
        hpo_terms : list[str], optional
            HPO term IDs relevant to patient
        top_k : int
            Number of top variants to analyze
        language : str
            "en" or "cn" for prompt language

        Returns
        -------
        dict
            {"candidates": [...], "differential_diagnosis": [...],
             "followup": [...], "warnings": [...], "raw_response": str}
        """
        top_df = variants_df.head(top_k)
        variant_summary = self._build_variant_summary(top_df)
        hpo_str = ", ".join(hpo_terms) if hpo_terms else "Not provided"
        system = SYSTEM_PROMPT_CN if language == "cn" else SYSTEM_PROMPT

        user_msg = self._build_prompt(
            variant_summary, patient_phenotype, hpo_str, language
        )

        logger.info(f"Model analysis: {len(top_df)} variants, {language=}")
        response = self._call_llm(system, user_msg)

        return self._parse_response(response, language)

    def _build_prompt(
        self, variant_summary: str,
        phenotype: str, hpo_str: str,
        language: str,
    ) -> str:
        if language == "cn":
            return f"""请分析以下候选变异：

患者表型：{phenotype}
HPO Terms：{hpo_str}

候选变异（按SeekRare评分排序）：
{variant_summary}

请用JSON格式返回，包含以下字段：
{{
  "candidates": [
    {{
      "rank": 1,
      "gene": "基因名",
      "variant": "c.HGVS",
      "acmg_class": "Pathogenic/Likely Pathogenic/VUS/Likely Benign/Benign",
      "acmg_criteria": ["PS1", "PM2", "PP3", ...],
      "pathogenicity_evidence": "主要证据总结",
      "phenotype_fit": "该变异对患者表型的解释程度",
      "confidence": "高/中/低",
      "notes": "任何需要补充的注意事项"
    }}
  ],
  "differential_diagnosis": ["鉴别诊断1", "鉴别诊断2", ...],
  "followup_studies": ["建议的随访研究1", ...],
  "warnings": ["需要特别注意的问题1", ...]
}}

确保JSON格式正确，可直接被Python json.loads解析。"""

        else:
            return f"""Analyze these candidate variants:

Patient Phenotype: {phenotype}
HPO Terms: {hpo_str}

Candidate Variants (ranked by SeekRare score):
{variant_summary}

Respond ONLY with valid JSON:
{{
  "candidates": [
    {{
      "rank": 1,
      "gene": "GENE_NAME",
      "variant": "c.HGVS",
      "acmg_class": "Pathogenic/Likely Pathogenic/VUS/Likely Benign/Benign",
      "acmg_criteria": ["PS1", "PM2", "PP3", ...],
      "pathogenicity_evidence": "Key supporting evidence",
      "phenotype_fit": "How well this variant explains patient phenotype",
      "confidence": "High/Medium/Low",
      "notes": "Additional caveats or notes"
    }}
  ],
  "differential_diagnosis": ["Differential diagnosis 1", ...],
  "followup_studies": ["Recommended functional studies", ...],
  "warnings": ["Caveats or concerns requiring attention", ...]
}}

Respond with ONLY JSON — no markdown, no explanation."""

    def _call_llm(self, system: str, user_msg: str) -> str:
        """Make LLM API call."""
        provider = self.config.provider.lower()

        if provider in ("openai", "local", "genos"):
            client = self._get_client()
            resp = client.chat.completions.create(
                model=self.config.model,
                messages=[
                    {"role": "system", "content": system},
                    {"role": "user", "content": user_msg},
                ],
                temperature=self.config.temperature,
                max_tokens=self.config.max_tokens,
            )
            return resp.choices[0].message.content

        elif provider == "anthropic":
            client = self._get_client()
            resp = client.messages.create(
                model=self.config.model,
                system=system,
                messages=[{"role": "user", "content": user_msg}],
                temperature=self.config.temperature,
                max_tokens=self.config.max_tokens,
            )
            return resp.content[0].text

        raise ValueError(f"Unknown provider: {provider}")

    def _parse_response(self, response: str, language: str) -> dict:
        """Parse LLM JSON response."""
        # Try to extract JSON from response
        try:
            # Handle markdown code blocks
            if "```json" in response:
                response = response.split("```json")[1].split("```")[0]
            elif "```" in response:
                response = response.split("```")[1].split("```")[0]
            elif response.strip().startswith("{"):
                response = response.strip()

            parsed = json.loads(response.strip())
            logger.info(f"Model analysis parsed successfully: {len(parsed.get('candidates', []))} candidates")
            parsed["raw_response"] = response
            return parsed

        except json.JSONDecodeError as e:
            logger.error(f"JSON parse failed: {e}")
            logger.debug(f"Response was: {response[:500]}")
            return {
                "candidates": [],
                "differential_diagnosis": [],
                "followup_studies": [],
                "warnings": [f"LLM response parse failed: {e}"],
                "raw_response": response,
            }

    def batch_analyze(
        self,
        variants_dfs: list[pd.DataFrame],
        phenotypes: list[str],
        **kwargs,
    ) -> list[dict]:
        """
        Batch analyze multiple patients.

        Parameters
        ----------
        variants_dfs : list[pd.DataFrame]
            One filtered variants DataFrame per patient
        phenotypes : list[str]
            Patient phenotype descriptions (parallel to variants_dfs)

        Returns
        -------
        list[dict]
            One analysis result dict per patient
        """
        results = []
        for i, (df, pheno) in enumerate(zip(variants_dfs, phenotypes)):
            logger.info(f"Batch analysis: patient {i+1}/{len(variants_dfs)}")
            result = self.analyze(df, pheno, **kwargs)
            result["patient_index"] = i
            results.append(result)

        return results
