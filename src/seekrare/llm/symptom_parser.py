"""
LLM Symptom Parser — the core of SeekRare's dual-dynamic scoring.

The LLM Symptom Parser takes a patient's free-text symptom description and outputs:
1. A list of relevant HPO terms with semantic relevance scores
2. A weight vector for annotation columns

These two outputs drive the dual-dynamic scoring engine.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Any, Optional

from loguru import logger

# Supported LLM providers
PROVIDER_TO_IMPORT = {
    "openai": ("openai", "OpenAI"),
    "anthropic": ("anthropic", "Anthropic"),
}


@dataclass
class LLMInterpretation:
    """Output from the LLM Symptom Parser."""
    relevant_hpos: list[dict]  # [{"hpo_id": "HP:0012345", "hpo_name": "...", "score": 0.95}, ...]
    weight_vector: dict[str, float]  # {"clinvar": 0.40, "hpo": 0.30, "cadd": 0.15, ...}
    reasoning: str  # LLM's reasoning for these weights


class LLMSymptomParser:
    """
    LLM-powered symptom interpreter.

    Parameters
    ----------
    provider : str
        "openai", "anthropic", or "local"
    model : str
        Model name (e.g., "gpt-4o", "claude-3-opus-20240229")
    api_key : str, optional
        API key. Falls back to environment variable.
    base_url : str, optional
        Base URL for OpenAI-compatible APIs (e.g., vLLM, Ollama)
    temperature : float
        Sampling temperature
    """

    SYSTEM_PROMPT = """You are a clinical geneticist assistant. Your task is to analyze a patient's symptom description and produce:

1. **Relevant HPO terms**: Extract the most relevant Human Phenotype Ontology (HPO) terms with semantic relevance scores (0.0-1.0).
2. **Annotation weight vector**: Assign importance weights to each variant annotation column. The weights must sum to 1.0.

Available annotation columns and their brief descriptions:
- clinvar_significance: ClinVar pathogenicity classification
- clinvar_stars: ClinVar expert panel review stars (0-4)
- cadd_score: CADD phred score (higher = more deleterious)
- gnomad_af: gnomAD allele frequency (lower = rarer)
- sift_score: SIFT prediction (0=damaging, 1=tolerated)
- polyphen_score: PolyPhen prediction (0=benign, 1=probably damaging)
- impact: Variant impact (HIGH/MODERATE/LOW)
- inheritance: Inferred inheritance pattern (AD/AR/denovo/unknown)
- hpo_terms: Associated HPO terms from annotation
- gene_constraint: Gene intolerance scores (pLI, oeLoF)

Important principles:
- For NEUROLOGICAL presentations: prioritize cadd_score, clinvar_significance, impact, hpo_terms
- For CARDIAC presentations: prioritize gnomad_af (keep very low), impact, gene_constraint
- For RARE symptoms: prioritize gnomad_af (lower is more important), hpo_terms
- For DEVELOPMENTAL presentations: prioritize hpo_terms, clinvar_significance, impact

Be precise. The weights you assign will directly determine which variants are ranked highest for this patient.
"""

    USER_PROMPT_TEMPLATE = """Analyze the following patient symptom description and produce HPO terms and weight vector.

Patient Symptoms:
{syptoms}

Please respond in JSON format:
{{
    "relevant_hpos": [
        {{"hpo_id": "HP:0001250", "hpo_name": "Seizure", "score": 0.95}},
        {{"hpo_id": "HP:0012469", "hpo_name": "Infantile spasms", "score": 0.88}},
        ...
    ],
    "weight_vector": {{
        "clinvar_significance": 0.35,
        "clinvar_stars": 0.10,
        "cadd_score": 0.20,
        "gnomad_af": 0.15,
        "hpo_terms": 0.15,
        "impact": 0.05
    }},
    "reasoning": "Brief explanation of why these weights were chosen..."
}}
"""

    def __init__(
        self,
        provider: str = "openai",
        model: str = "gpt-4o",
        api_key: Optional[str] = None,
        base_url: Optional[str] = None,
        temperature: float = 0.0,
    ):
        self.provider = provider
        self.model = model
        self.temperature = temperature

        # Resolve API key
        if api_key:
            self.api_key = api_key
        elif provider == "openai":
            self.api_key = os.getenv("OPENAI_API_KEY") or os.getenv("ANTHROPIC_API_KEY")
        elif provider == "anthropic":
            self.api_key = os.getenv("ANTHROPIC_API_KEY")
        else:
            self.api_key = None

        self.base_url = base_url
        self._client = self._init_client()

    def _init_client(self):
        """Initialize the LLM client based on provider."""
        if self.provider == "openai":
            try:
                from openai import OpenAI
            except ImportError:
                raise ImportError("Please install openai: pip install openai")
            return OpenAI(api_key=self.api_key, base_url=self.base_url)

        elif self.provider == "anthropic":
            try:
                from anthropic import Anthropic
            except ImportError:
                raise ImportError("Please install anthropic: pip install anthropic")
            return Anthropic(api_key=self.api_key)

        elif self.provider == "local":
            # Local model via OpenAI-compatible API
            from openai import OpenAI
            return OpenAI(api_key="not-needed", base_url=self.base_url or "http://localhost:8000/v1")

        else:
            raise ValueError(f"Unknown provider: {provider}")

    def interpret(self, symptoms: str) -> dict[str, Any]:
        """
        Interpret patient symptoms and return LLM interpretation.

        Parameters
        ----------
        symptoms : str
            Free-text patient symptom description

        Returns
        -------
        dict
            {
                "relevant_hpos": [...],
                "weight_vector": {...},
                "reasoning": "..."
            }
        """
        logger.info(f"LLM interpreting symptoms (provider={self.provider}, model={self.model})")

        user_prompt = self.USER_PROMPT_TEMPLATE.format(symptoms=symptoms)

        if self.provider == "openai" or self.provider == "local":
            response = self._client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": self.SYSTEM_PROMPT},
                    {"role": "user", "content": user_prompt},
                ],
                response_format={"type": "json_object"},
                temperature=self.temperature,
            )
            content = response.choices[0].message.content

        elif self.provider == "anthropic":
            response = self._client.messages.create(
                model=self.model,
                max_tokens=2048,
                system=self.SYSTEM_PROMPT,
                messages=[{"role": "user", "content": user_prompt}],
            )
            content = response.content[0].text

        else:
            raise ValueError(f"Unknown provider: {self.provider}")

        try:
            parsed = json.loads(content)
            logger.info(f"LLM output: {len(parsed.get('relevant_hpos', []))} HPO terms, "
                        f"weights sum={sum(parsed.get('weight_vector', {}).values()):.2f}")
            return parsed
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse LLM JSON output: {e}\nContent: {content[:500]}")
            # Return safe fallback
            return {
                "relevant_hpos": [],
                "weight_vector": {
                    "clinvar_significance": 0.30,
                    "clinvar_stars": 0.15,
                    "cadd_score": 0.25,
                    "gnomad_af": 0.10,
                    "hpo_terms": 0.15,
                    "impact": 0.05,
                },
                "reasoning": "Fallback due to JSON parse error",
            }
