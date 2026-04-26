"""
Genos LLM integration — clinical genetics domain-specific LLM.

Genos is a specialized LLM for clinical genetics that can:
- Interpret complex variant annotations
- Assess variant pathogenicity based on clinical context
- Integrate multiple annotation sources for 综合分析
- Generate clinical reports in Chinese/English

Usage:
    from seekrare.llm import GenosClient
    genos = GenosClient(api_key="...", base_url="https://genos.example.com/v1")
    response = genos.analyze_variant(variant_info, patient_context)
    report = genos.generate_clinical_report(candidate_variants, patient_phenotype)
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Any, Optional, Union
from loguru import logger


@dataclass
class GenosConfig:
    """Configuration for Genos LLM."""
    model: str = "genos-clinical-v1"
    temperature: float = 0.0
    max_tokens: int = 4096
    api_base: str = "https://api.genos.tech/v1"   # Placeholder - replace with actual


SYSTEM_PROMPT = """You are Genos, a specialized clinical genetics AI assistant.
You analyze genomic variants in the context of patient phenotypes to identify
the most likely disease-causing mutations.

Your capabilities:
- Variant pathogenicity assessment (ACMG guidelines)
- Integration of ClinVar, OMIM, HGMD, UniProt annotations
- Phenotype-genotype correlation analysis using HPO terms
- Protein structure impact prediction
- Literature evidence synthesis
- Chinese and English clinical report generation

Be precise, cite evidence, and distinguish between definite pathogenic
variants and variants requiring further functional validation.
"""


class GenosClient:
    """
    Genos LLM API client.

    Parameters
    ----------
    api_key : str
        API key for Genos service
    base_url : str, optional
        API base URL. Defaults to Genos official API.
        For self-hosted: "http://localhost:8000/v1"
    model : str
        Model name to use
    temperature : float
        Sampling temperature
    """

    def __init__(
        self,
        api_key: Optional[str] = None,
        base_url: Optional[str] = None,
        model: str = "genos-clinical-v1",
        temperature: float = 0.0,
    ):
        self.api_key = api_key or os.getenv("GENOS_API_KEY", "")
        self.base_url = (base_url or os.getenv("GENOS_API_BASE", "")).rstrip("/")
        self.model = model
        self.temperature = temperature

        if not self.api_key:
            logger.warning("No Genos API key provided. Set GENOS_API_KEY env var.")

    def _call(self, messages: list[dict], **kwargs) -> dict:
        """Make API call to Genos."""
        try:
            from openai import OpenAI
        except ImportError:
            raise ImportError("openai package required: pip install openai")

        client = OpenAI(api_key=self.api_key, base_url=self.base_url or None)

        params = {
            "model": self.model,
            "messages": messages,
            "temperature": kwargs.get("temperature", self.temperature),
            "max_tokens": kwargs.get("max_tokens", 4096),
        }

        if kwargs.get("response_format"):
            params["response_format"] = kwargs["response_format"]

        response = client.chat.completions.create(**params)
        return {"content": response.choices[0].message.content}

    def analyze_variant(
        self,
        variant_info: dict,
        patient_phenotype: str,
        hpo_terms: Optional[list[str]] = None,
    ) -> dict:
        """
        Analyze a single variant in clinical context.

        Parameters
        ----------
        variant_info : dict
            Variant annotations: gene, cDNA_change, protein_change, clinvar_sig,
            cadd_score, gnomad_af, impact, sift, polyphen, etc.
        patient_phenotype : str
            Free-text patient phenotype description
        hpo_terms : list[str], optional
            Relevant HPO term IDs

        Returns
        -------
        dict
            Analysis with: pathogenicity_assessment, evidence, acmg_criteria, recommendation
        """
        variant_str = json.dumps(variant_info, indent=2, ensure_ascii=False)
        hpo_str = ", ".join(hpo_terms) if hpo_terms else "Not provided"

        user_msg = f"""Analyze this variant in clinical context:

Patient Phenotype: {patient_phenotype}
HPO Terms: {hpo_str}

Variant Information:
{variant_str}

Provide a structured analysis including:
1. Pathogenicity assessment (Pathogenic/Likely Pathogenic/VUS/Likely Benign/Benign)
2. Key evidence supporting this assessment
3. ACMG criteria applied (PVS1, PS1, PM2, PP3, etc.)
4. Whether this variant likely explains the patient phenotype
5. Recommended follow-up (functional studies, segregation analysis, etc.)

Respond in JSON format."""
        messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": user_msg},
        ]

        result = self._call(messages, response_format={"type": "json_object"})
        try:
            return json.loads(result["content"])
        except Exception:
            logger.error(f"Genos JSON parse failed: {result['content'][:200]}")
            return {"raw": result["content"]}

    def batch_analyze(
        self,
        variants_df,
        patient_phenotype: str,
        hpo_terms: Optional[list[str]] = None,
        top_n: int = 20,
    ) -> list[dict]:
        """
        Batch analyze top candidate variants.

        Parameters
        ----------
        variants_df : pd.DataFrame
            Candidate variants DataFrame (from SeekRare pipeline)
        patient_phenotype : str
            Patient phenotype description
        hpo_terms : list[str], optional
        top_n : int
            Number of top variants to analyze

        Returns
        -------
        list[dict]
            List of analysis results per variant
        """
        top_df = variants_df.head(top_n)

        results = []
        for i, row in top_df.iterrows():
            variant_info = row.to_dict()
            logger.info(f"Analyzing variant {i+1}/{len(top_df)}: {row.get('gene_name', '?')}")

            analysis = self.analyze_variant(
                variant_info=variant_info,
                patient_phenotype=patient_phenotype,
                hpo_terms=hpo_terms,
            )
            analysis["variant_index"] = i
            analysis["gene"] = row.get("gene_name", "")
            results.append(analysis)

        return results

    def generate_clinical_report(
        self,
        candidate_variants: list[dict],
        patient_info: dict,
        phenotype: str,
        llm_interpretation: Optional[dict] = None,
    ) -> str:
        """
        Generate a comprehensive clinical report.

        Parameters
        ----------
        candidate_variants : list[dict]
            List of variant analysis results from batch_analyze()
        patient_info : dict
            Patient demographics: age, sex, family_history, etc.
        phenotype : str
            Patient phenotype description
        llm_interpretation : dict, optional
            LLM symptom interpretation output from SeekRare

        Returns
        -------
        str
            Formatted clinical report in Chinese/English
        """
        variants_str = json.dumps(candidate_variants[:10], indent=2, ensure_ascii=False)
        patient_str = json.dumps(patient_info, indent=2, ensure_ascii=False)

        user_msg = f"""Generate a clinical genetics report for this patient.

Patient Information:
{patient_str}

Phenotype: {phenotype}

Top Candidate Variants (analyzed):
{variants_str}

{f'LLM Phenotype Interpretation: {json.dumps(llm_interpretation, indent=2, ensure_ascii=False)}' if llm_interpretation else ''}

Generate a comprehensive clinical report in Chinese with:
1. Executive Summary (结论概要)
2. Patient demographics and presenting phenotype
3. Analysis methodology
4. Top candidate variants with detailed findings
5. ACMG pathogenicity classification for each candidate
6. Genotype-phenotype correlation analysis
7. Differential diagnosis considerations
8. Recommended genetic counseling and follow-up
9. Literature references (if applicable)

Use professional clinical genetics terminology.
Be conservative and distinguish evidence-based conclusions from speculation."""
        messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": user_msg},
        ]

        result = self._call(messages, max_tokens=8192)
        return result["content"]

    def assess_compound_heterozygous(
        self,
        gene: str,
        variants: list[dict],
        patient_phenotype: str,
    ) -> dict:
        """
        Assess compound heterozygous variant pairs in a gene.

        Parameters
        ----------
        gene : str
            Gene symbol
        variants : list[dict]
            List of variants in the gene (each with cDNA_change, protein_change, etc.)
        patient_phenotype : str
            Patient phenotype

        Returns
        -------
        dict
            Compound het assessment: combined_pathogenicity, pair_evidence, recommendation
        """
        vars_str = json.dumps(variants, indent=2, ensure_ascii=False)

        user_msg = f"""Assess compound heterozygous variant pairs in {gene}:

Patient Phenotype: {patient_phenotype}

Variants in {gene}:
{vars_str}

For compound heterozygous inheritance:
1. Identify which variant pairs are on opposite chromosomes (trans configuration)
2. Assess combined pathogenicity of each pair
3. Determine if any pair sufficiently explains the phenotype
4. Consider allele balance, conservation, and functional studies

Respond in JSON format."""
        messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": user_msg},
        ]

        result = self._call(messages, response_format={"type": "json_object"})
        try:
            return json.loads(result["content"])
        except Exception:
            return {"raw": result["content"]}
