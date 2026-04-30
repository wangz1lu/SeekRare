"""
stage3_annotate.py — Stage 3: LLM 打分 + 排序

工作流:
1. stage3_prep.summarize_stage1() — 解析标签，统计 unique 值
2. LLM 接收 prompt → 输出 JSON 打分
3. 本地 Python 计算每行分数
4. 排序输出 top-K
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger

from seekrare.scoring.stage3_prep import (
    summarize_stage1,
    build_llm_prompt,
    parse_clndisdb_tags,
    parse_clndn_tags,
)


class Stage3Scorer:
    """
    Stage 3: LLM 驱动的变异排序。

    参数
    ----
    csv_path : str
        Stage 1 输出 CSV (3_clinvar_annotated.csv)
    symptoms : str
        患者症状描述（自由文本）
    top_k : int
        返回 top-K 候选，默认 50
    llm_provider : str
        "openai" / "anthropic"
    llm_model : str
        模型名
    api_key : str, optional
    base_url : str, optional
    """

    def __init__(
        self,
        csv_path: str,
        symptoms: str,
        top_k: int = 50,
        llm_provider: str = "openai",
        llm_model: str = "gpt-4o",
        api_key: Optional[str] = None,
        base_url: Optional[str] = None,
    ):
        self.csv_path = str(csv_path)
        self.symptoms = symptoms
        self.top_k = top_k
        self.llm_provider = llm_provider
        self.llm_model = llm_model
        self.api_key = api_key
        self.base_url = base_url
        self._scores: Optional[dict] = None

    # ── LLM 调用 ──────────────────────────────────────────────────────────

    def _call_llm(self, prompt: str) -> dict:
        """调用 LLM 获取 JSON 格式打分。"""
        if self.llm_provider == "openai":
            from openai import OpenAI
            client = OpenAI(api_key=self.api_key or "", base_url=self.base_url)
            resp = client.chat.completions.create(
                model=self.llm_model,
                messages=[
                    {"role": "system",
                     "content": "You are a clinical genetics scoring assistant. "
                               "Return ONLY valid JSON matching the requested schema."},
                    {"role": "user", "content": prompt},
                ],
                temperature=0.0,
                max_tokens=4096,
            )
            content = resp.choices[0].message.content

        elif self.llm_provider == "anthropic":
            import anthropic
            client = anthropic.Anthropic(api_key=self.api_key or "")
            resp = client.messages.create(
                model=self.llm_model,
                system="You are a clinical genetics scoring assistant. Return ONLY valid JSON.",
                messages=[{"role": "user", "content": prompt}],
                temperature=0.0,
                max_tokens=4096,
            )
            content = resp.content[0].text

        else:
            raise ValueError(f"Unknown provider: {self.llm_provider}")

        # 提取 JSON
        content = content.strip()
        if content.startswith("```json"):
            content = content[7:]
        if content.startswith("```"):
            content = re.sub(r"^```[a-z]*\n?", "", content, count=1)
        content = content.strip().rstrip("```").rstrip()

        return json.loads(content)

    # ── 打分计算 ────────────────────────────────────────────────────────────

    @staticmethod
    def _score_feature_type(val: str, score_map: dict) -> float:
        if pd.isna(val) or str(val).strip() == "":
            return 0.0
        return score_map.get(str(val).strip(), 0.0)

    @staticmethod
    def _score_clnsig(val: str, score_map: dict) -> float:
        if pd.isna(val) or str(val).strip() == "":
            return 0.0
        return score_map.get(str(val).strip(), 0.0)

    @staticmethod
    def _score_clnrevstat(val: str, score_map: dict) -> float:
        if pd.isna(val) or str(val).strip() == "":
            return 0.0
        return score_map.get(str(val).strip(), 0.0)

    @staticmethod
    def _score_tags(val: str, score_map: dict, parser_fn) -> float:
        """
        解析 CLNDISDB / CLNDN 标签，取该行所有标签中最高分。
        """
        tags = parser_fn(val)
        if not tags:
            return 0.0
        scores = [score_map.get(t, 0.0) for t in tags]
        return max(scores) if scores else 0.0

    def score_row(self, row: pd.Series, weights: dict, score_maps: dict) -> float:
        """计算单行加权总分。"""
        w = weights
        m = score_maps

        s_ft = self._score_feature_type(row.get("feature_type", ""), m["feature_type"])
        s_sig = self._score_clnsig(row.get("CLNSIG", ""), m["clnsig"])
        s_rs = self._score_clnrevstat(row.get("CLNREVSTAT", ""), m["clnrevstat"])
        s_db = self._score_tags(row.get("CLNDISDB", ""), m["clndisdb_tags"], parse_clndisdb_tags)
        s_dn = self._score_tags(row.get("CLNDN", ""), m["clndn_tags"], parse_clndn_tags)

        return (
            s_ft * w.get("feature_type", 0.0)
            + s_sig * w.get("CLNSIG", 0.0)
            + s_rs * w.get("CLNREVSTAT", 0.0)
            + s_db * w.get("CLNDISDB", 0.0)
            + s_dn * w.get("CLNDN", 0.0)
        )

    # ── 主流程 ─────────────────────────────────────────────────────────────

    def run(self) -> pd.DataFrame:
        """
        执行 Stage 3 全流程。

        Returns
        -------
        pd.DataFrame
            排序后的 top-K 候选变异（含 seekrare_score 列）
        """
        logger.info("Stage 3: LLM Scoring & Ranking")
        logger.info(f"  Input: {self.csv_path}")
        logger.info(f"  Symptoms: {self.symptoms}")

        # ── 1. 数据摘要 ─────────────────────────────────────────────────────
        summary = summarize_stage1(self.csv_path)
        logger.info(f"  Summary: {summary['n_total']:,} rows, {summary['n_with_clinvar']:,} with ClinVar")
        logger.info(f"  CLNDISDB tags: {len(summary['clndisdb_tag_counts'])}")
        logger.info(f"  CLNDN tags: {len(summary['clndn_tag_counts'])}")

        # ── 2. LLM 打分 prompt ───────────────────────────────────────────────
        prompt = build_llm_prompt(self.symptoms, summary)
        logger.info("  Sending prompt to LLM...")

        raw_scores = self._call_llm(prompt)
        logger.info(f"  LLM response keys: {list(raw_scores.keys())}")

        weights = raw_scores.get("col_weights", {})
        score_maps = {
            "feature_type": raw_scores.get("feature_type_scores", {}),
            "clnsig": raw_scores.get("clnsig_scores", {}),
            "clnrevstat": raw_scores.get("clnrevstat_scores", {}),
            "clndisdb_tags": raw_scores.get("clndisdb_tag_scores", {}),
            "clndn_tags": raw_scores.get("clndn_tag_scores", {}),
        }

        self._scores = {"weights": weights, "maps": score_maps}

        # ── 3. 本地打分 ─────────────────────────────────────────────────────
        logger.info("  Scoring all rows locally...")
        df = pd.read_csv(self.csv_path, dtype=str, low_memory=False)
        df["seekrare_score"] = df.apply(
            lambda row: self.score_row(row, weights, score_maps),
            axis=1,
        )

        # ── 4. 排序 ──────────────────────────────────────────────────────────
        df = df.sort_values("seekrare_score", ascending=False)
        df["rank"] = range(1, len(df) + 1)

        # ── 5. 输出 top-K ───────────────────────────────────────────────────
        result = df.head(self.top_k).copy()
        logger.info(f"  Stage 3 完成: top {len(result)} candidates")
        logger.info(f"  Top 3: {result[['CPRA','gene_name','seekrare_score']].head(3).to_string()}")

        return result

    def save(self, df: pd.DataFrame, output_path: str):
        """保存结果。"""
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, index=False)
        logger.info(f"  Stage 3 结果已保存: {output_path}")
