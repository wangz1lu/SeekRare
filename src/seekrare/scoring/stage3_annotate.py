"""
stage3_annotate.py — Stage 3: LLM 动态打分 + 排序

工作流:
1. summarize_stage1() — 解析标签，统计 unique 值
2. LLM 接收 prompt → 输出 JSON 打分（动态列的 per-value scores + 权重）
3. 本地 Python 计算每行分数（静态列直接查表，动态列 LLM 打分）
4. 排序输出 top-K

动态评分列（LLM 给出 per-value score）:
  gene_name, HPO, OMIM, Orphanet, inheritance_mode, MC

固定评分列（代码内置映射表）:
  feature_type:  CDS=1, exon=0.9, gene=0.7, start_codon=0.8, stop_codon=0.8, transcript=0.5
  significance:  Pathogenic=1, Likely_Pathogenic=0.85, Uncertain_significance=0.5,
                Conflicting_classifications_of_pathogenicity=0.4, Likely_benign=0.1, Benign=0
                ( "/" 分隔时取最坏情况)
  clinvarstar:   0~5 直接映射数字
  eqtl_tissue:   有内容=0.5, 空=0
  splicevardb:    Splice-altering=1, Low-frequency=0.6, Conflicting=0.5, Normal=0.2, 其他=0
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


# ─────────────────────────────────────────────────────────────────────────────
# 固定打分表（内置，LLM 不参与）
# ─────────────────────────────────────────────────────────────────────────────

FEATURE_TYPE_MAP = {
    "CDS": 1.0,
    "exon": 0.9,
    "gene": 0.7,
    "start_codon": 0.8,
    "stop_codon": 0.8,
    "transcript": 0.5,
}


SIGNIFICANCE_WORST = {
    # 分数越低越良性，分数越高越致病；"/" 分隔取最低分（最良性 = 最坏情况）
    "Benign": -1.0,
    "Likely_benign": -0.5,
    "Conflicting_classifications_of_pathogenicity": 0.4,
    "Uncertain_significance": 0.5,
    "Likely_Pathogenic": 0.85,
    "Pathogenic": 1.0,
    "other": 0.5,
}


def _score_significance(val: str) -> float:
    """
    处理 "/" 分隔的多值 significance。
    逻辑：按 significance 方向分两段，分别取极值：
      - 若包含 Benign → 取最低分（最良性）
      - 若包含 Pathogenic → 取最高分（最致病）
      - 否则取最高分（最坏情况）
    无 ClinVar 注释 → -0.5 惩罚。
    """
    if pd.isna(val) or str(val).strip() == "":
        return -0.5
    parts = str(val).strip().split("/")
    scores = []
    for p in parts:
        p = p.strip()
        scores.append(SIGNIFICANCE_WORST.get(p, 0.5))

    # 若有任何 Benign 类 → 取最低（最良性 = 最坏结果）
    benign_keywords = ["benign"]
    has_benign = any(
        any(bk in str(p).lower() for bk in benign_keywords) for p in parts
    )
    if has_benign:
        return min(scores)   # Benign/Likely_benign → -1.0

    # 否则取最高（Pathogenic/Likely_Pathogenic → 1.0）
    return max(scores) if scores else -0.5


CLINVARSTAR_MAP = {
    "0": 0.0, "1": 1.0, "2": 2.0, "3": 3.0, "4": 4.0, "5": 5.0,
    0: 0.0, 1: 1.0, 2: 2.0, 3: 3.0, 4: 4.0, 5: 5.0,
}


def _score_clinvarstar(val) -> float:
    if pd.isna(val):
        return 0.0
    return CLINVARSTAR_MAP.get(str(val).strip(), 0.0)


SPLICEVARDB_MAP = {
    "Splice-altering": 1.0,
    "Low-frequency": 0.6,
    "Conflicting": 0.5,
    "Normal": 0.2,
}


def _score_splicevardb(val: str) -> float:
    if pd.isna(val) or str(val).strip() == "":
        return 0.0
    return SPLICEVARDB_MAP.get(str(val).strip(), 0.0)


def _score_eqtl_tissue(val: str) -> float:
    """有内容（组织名）就给分"""
    if pd.isna(val) or str(val).strip() == "":
        return 0.0
    return 0.5


# ─────────────────────────────────────────────────────────────────────────────
# Stage3Scorer
# ─────────────────────────────────────────────────────────────────────────────

class Stage3Scorer:
    """
    Stage 3: LLM 驱动的变异排序。

    参数
    ----
    csv_path : str
        Stage 1/2 输出 CSV
    symptoms : str
        患者症状描述（自由文本）
        患者性别 "male" / "female"（影响 inheritance_mode 打分）
    top_k : int
        返回 top-K 候选，默认 50
    llm_provider : str
        "openai" / "anthropic"
    llm_model : str
        模型名
    api_key : str, optional
    base_url : str, optional
    """

    # 固定列（内置映射表）
    STATIC_COLS = [
        "feature_type", "significance", "clinvarstar",
        "eqtl_tissue", "splicevardb",
    ]

    # 动态列（LLM 给出 per-value scores）
    DYNAMIC_COLS = [
        "gene_name", "HPO", "OMIM", "Orphanet", "diseasename",
        "inheritance_mode", "MC",
    ]

    def __init__(
        self,
        csv_path: str,
        symptoms: str,
        top_k: int = 50,
        llm_provider: str = "openai",
        llm_model: str = "deepseek-v4-flash",
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
                max_tokens=16384,
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
                max_tokens=16384,
            )
            content = resp.content[0].text

        else:
            raise ValueError(f"Unknown provider: {self.llm_provider}")

        content = content.strip()
        if content.startswith("```json"):
            content = content[7:]
        if content.startswith("```"):
            content = re.sub(r"^```[a-z]*\n?", "", content, count=1)
        content = content.strip().rstrip("```").rstrip()

        return json.loads(content)

    # ── 打分计算 ────────────────────────────────────────────────────────────

    @staticmethod
    def _score_gene_name(val: str, score_map: dict) -> float:
        if pd.isna(val) or str(val).strip() == "":
            return 0.0
        return score_map.get(str(val).strip(), 0.0)

    @staticmethod
    def _score_hpo(val: str, score_map: dict) -> float:
        """解析 HPO 标签（HP:xxxx），取最高分"""
        if pd.isna(val) or str(val).strip() == "":
            return 0.0
        tags = str(val).strip().split(";")
        scores = [score_map.get(t.strip(), 0.0) for t in tags if t.strip()]
        return max(scores) if scores else 0.0

    @staticmethod
    def _score_omim(val: str, score_map: dict) -> float:
        """解析 OMIM 标签，取最高分"""
        if pd.isna(val) or str(val).strip() == "":
            return 0.0
        tags = str(val).strip().split(";")
        scores = [score_map.get(t.strip(), 0.0) for t in tags if t.strip()]
        return max(scores) if scores else 0.0

    @staticmethod
    def _score_orphanet(val: str, score_map: dict) -> float:
        """解析 Orphanet 标签，取最高分"""
        if pd.isna(val) or str(val).strip() == "":
            return 0.0
        tags = str(val).strip().split(";")
        scores = [score_map.get(t.strip(), 0.0) for t in tags if t.strip()]
        return max(scores) if scores else 0.0

    @staticmethod
    def _score_diseasename(val: str, score_map: dict) -> float:
        """
        解析 disease_name 列（用 | 或 ; 分隔多个疾病名），取最高分。
        """
        if pd.isna(val) or str(val).strip() == "":
            return 0.0
        parts = re.split(r"[|;]", str(val).strip())
        scores = [score_map.get(p.strip(), 0.0) for p in parts if p.strip()]
        return max(scores) if scores else 0.0

    @staticmethod
    def _score_inheritance_mode(val: str, score_map: dict) -> float:
        if pd.isna(val) or str(val).strip() == "":
            return 0.0
        return score_map.get(str(val).strip(), 0.0)

    @staticmethod
    def _score_mc(val: str, score_map: dict) -> float:
        """MC 取值可能是 SO:xxxx 格式，直接查 score_map"""
        if pd.isna(val) or str(val).strip() == "":
            return 0.0
        return score_map.get(str(val).strip(), 0.0)

    def score_row(self, row: pd.Series, weights: dict, dynamic_maps: dict) -> float:
        """计算单行加权总分（静态列 + 动态列）"""
        w = weights
        dm = dynamic_maps

        # ── 静态列（内置映射表）──────────────────────────────
        s_ft = FEATURE_TYPE_MAP.get(str(row.get("feature_type", "")).strip(), 0.0)
        s_sig = _score_significance(row.get("significance", ""))
        s_star = _score_clinvarstar(row.get("clinvarstar", ""))
        s_eqtl = _score_eqtl_tissue(row.get("eqtl_tissue", ""))
        s_sp = _score_splicevardb(row.get("splicevardb", ""))

        # ── 动态列（LLM 打分）───────────────────────────────
        s_gene = self._score_gene_name(row.get("gene_name", ""), dm.get("gene_name", {}))
        s_hpo = self._score_hpo(row.get("HPO", ""), dm.get("HPO", {}))
        s_omim = self._score_omim(row.get("OMIM", ""), dm.get("OMIM", {}))
        s_orphan = self._score_orphanet(row.get("Orphanet", ""), dm.get("Orphanet", {}))
        s_disease = self._score_diseasename(row.get("diseasename", ""), dm.get("diseasename", {}))
        s_inh = self._score_inheritance_mode(row.get("inheritance_mode", ""), dm.get("inheritance_mode", {}))
        s_mc = self._score_mc(row.get("MC", ""), dm.get("MC", {}))

        # ── gene × disease 一致性奖励 ──────────────────────────
        # 若 gene_name 和 diseasename 都与症状高度相关（>0.6），给额外奖励
        consistency_bonus = 0.0
        if s_gene > 0.6 and s_disease > 0.6:
            consistency_bonus = (s_gene * s_disease) * 0.1  # 最高 +0.1

        return (
            s_ft    * w.get("feature_type", 0.0)
            + s_sig * w.get("significance", 0.0)
            + s_star * w.get("clinvarstar", 0.0)
            + s_eqtl * w.get("eqtl_tissue", 0.0)
            + s_sp   * w.get("splicevardb", 0.0)
            + s_gene  * w.get("gene_name", 0.0)
            + s_hpo   * w.get("HPO", 0.0)
            + s_omim  * w.get("OMIM", 0.0)
            + s_orphan * w.get("Orphanet", 0.0)
            + s_disease * w.get("diseasename", 0.0)
            + s_inh   * w.get("inheritance_mode", 0.0)
            + s_mc    * w.get("MC", 0.0)
            + consistency_bonus
        )

    # ── 主流程 ─────────────────────────────────────────────────────────────

    def run(self) -> pd.DataFrame:
        """执行 Stage 3 全流程。"""
        logger.info("Stage 3: LLM Dynamic Scoring & Ranking")
        logger.info(f"  Input: {self.csv_path}")
        logger.info(f"  Symptoms: {self.symptoms}")

        # ── 1. 数据摘要（收集动态列 unique 值）────────────────────────────
        summary = summarize_stage1(self.csv_path)
        logger.info(f"  Summary: {summary['n_total']:,} rows")

        # ── 2. LLM 打分 prompt ─────────────────────────────────────────────
        prompt = build_llm_prompt_new(self.symptoms, summary)
        logger.info("  Sending prompt to LLM...")

        raw = self._call_llm(prompt)
        logger.info(f"  LLM response keys: {list(raw.keys())}")

        weights = raw.get("col_weights", {})

        # 动态列 per-value score maps
        dynamic_maps = {
            "gene_name":        raw.get("gene_name_scores", {}),
            "HPO":              raw.get("HPO_scores", {}),
            "OMIM":             raw.get("OMIM_scores", {}),
            "Orphanet":         raw.get("Orphanet_scores", {}),
            "diseasename":      raw.get("diseasename_scores", {}),
            "inheritance_mode": raw.get("inheritance_mode_scores", {}),
            "MC":               raw.get("MC_scores", {}),
        }

        self._scores = {"weights": weights, "maps": dynamic_maps}

        # ── 3. 本地打分 ─────────────────────────────────────────────────────
        logger.info("  Scoring all rows locally...")
        df = pd.read_csv(self.csv_path, dtype=str, low_memory=False)
        df["seekrare_score"] = df.apply(
            lambda row: self.score_row(row, weights, dynamic_maps),
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
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, index=False)
        logger.info(f"  Stage 3 结果已保存: {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 新版 prompt builder（适配新列体系）
# ─────────────────────────────────────────────────────────────────────────────

def build_llm_prompt_new(symptoms: str, summary: dict) -> str:
    """
    生成新版 LLM prompt：
      - 静态列说明（内置打分，无需 LLM 输出）
      - 动态列 unique 值（LLM 输出 per-value scores）
      - col_weights（动态列权重，静态列权重固定）
    """
    # 收集各动态列 unique tags
    gene_names = list(summary.get("gene_name_counts", {}).keys())[:100]
    hpo_tags = list(summary.get("HPO_tag_counts", {}).keys())[:100]
    omim_tags = list(summary.get("OMIM_tag_counts", {}).keys())[:100]
    orphanet_tags = list(summary.get("Orphanet_tag_counts", {}).keys())[:100]
    inh_modes = list(summary.get("inheritance_mode_counts", {}).keys())
    mc_values = list(summary.get("MC_counts", {}).keys())[:50]

    lines = [
        f"患者症状: {symptoms}",
        f"",
        f"变异数据统计（共 {summary['n_total']:,} 行）:",
        f"",
        f"【静态评分列】（代码内置映射表，LLM 不需输出这些的分数）:",
        f"  feature_type:  CDS=1.0, exon=0.9, gene=0.7, start_codon=0.8, stop_codon=0.8, transcript=0.5",
        f"  significance: '/' 分隔时取最坏情况。Pathogenic=1.0, Likely_Pathogenic=0.85,",
        f"                Uncertain_significance=0.5, Conflicting=0.4, Likely_benign=0.1, Benign=0.0",
        f"  clinvarstar:  直接用星级数字 0~5",
        f"  eqtl_tissue:  有内容=0.5，无内容=0",
        f"  splicevardb:  Splice-altering=1.0, Low-frequency=0.6, Conflicting=0.5, Normal=0.2",
        f"",
        f"【动态评分列】（LLM 必须输出每列的 per-value scores 和权重）:",
        f"",
        f"【gene_name unique 值】({len(gene_names)} 个，展示部分):",
    ]
    for g in gene_names[:30]:
        lines.append(f"  {g}")

    lines.extend([
        f"",
        f"【HPO 标签分布】(HP:xxxx 格式，展示部分):",
    ])
    hpo_counts = summary.get("HPO_tag_counts", {})
    for tag, cnt in list(hpo_counts.items())[:50]:
        lines.append(f"  [{cnt}] {tag}")

    lines.extend([
        f"",
        f"【OMIM 标签分布】(展示部分):",
    ])
    omim_counts = summary.get("OMIM_tag_counts", {})
    for tag, cnt in list(omim_counts.items())[:50]:
        lines.append(f"  [{cnt}] {tag}")

    lines.extend([
        f"",
        f"【Orphanet 标签分布】(展示部分):",
    ])
    orphanet_counts = summary.get("Orphanet_tag_counts", {})
    for tag, cnt in list(orphanet_counts.items())[:50]:
        lines.append(f"  [{cnt}] {tag}")

    lines.extend([
        f"",
        f"【diseasename 唯一值】（部分展示）:",
    ])
    disease_names = list(summary.get("diseasename_counts", {}).keys())[:50]
    for dn in disease_names[:30]:
        lines.append(f"  {dn}")

    lines.extend([
        f"",
        f"【inheritance_mode 分布】:",
    ])
    inh_counts = summary.get("inheritance_mode_counts", {})
    for k, v in inh_counts.items():
        lines.append(f"  {k}: {v}")

    lines.extend([
        f"",
        f"【MC 值分布】(SO:xxxx 格式，展示部分):",
    ])
    mc_counts = summary.get("MC_counts", {})
    for k, v in list(mc_counts.items())[:50]:
        lines.append(f"  [{v}] {k}")

    lines.extend([
        f"",
        f"请根据患者症状，为以上动态列的每个取值给出 0~1 的相关性分数（1=最相关，0=不相关）:",
        f"",
        f"  1. gene_name_scores: {{基因名: 分数, ...}}",
        f"  2. HPO_scores: {{HP:xxxx: 分数, ...}}",
        f"  3. OMIM_scores: {{OMIM:xxxxx: 分数, ...}}",
        f"  4. Orphanet_scores: {{Orphanet:xxxxx: 分数, ...}}",
        f"  5. inheritance_mode_scores: {{denovo/recessive/xlinked: 分数}}",
        f"     (完全根据疾病/症状特征判断：如先天性发育异常多为de_novo，地中海贫血多为recessive，DMD等多为xlinked)",
        f"  6. MC_scores: {{SO:xxxx: 分数, ...}}",
        f"  7. col_weights: 各列权重（归一化为 sum=1.0，只给动态列权重，静态列权重固定）",
        f"",
        f"返回 JSON 格式：",
        f'''{{
  "gene_name_scores": {{"SASS6": 0.9, "DBT": 0.8, ...}},
  "HPO_scores": {{"HP:0004321": 0.95, "HP:0001513": 0.7, ...}},
  "OMIM_scores": {{"OMIM:616126": 0.9, ...}},
  "Orphanet_scores": {{"Orphanet:319563": 0.8, ...}},
  "diseasename_scores": {{"Thalassemia": 0.95, "Sickle cell disease": 0.9, ...}},
  "inheritance_mode_scores": {{"de_novo": 0.9, "recessive": 0.7, "xlinked": 0.3}},
  "MC_scores": {{"SO:0001583": 0.7, "SO:0001627": 0.9, ...}},
  "col_weights": {{
    "gene_name": 0.15, "HPO": 0.15, "OMIM": 0.15,
    "Orphanet": 0.08, "diseasename": 0.12, "inheritance_mode": 0.10, "MC": 0.05,
    "feature_type": 0.10, "significance": 0.10, "clinvarstar": 0.025,
    "eqtl_tissue": 0.025, "splicevardb": 0.025
  }}
}}''',
    ])

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# 兼容旧接口
# ─────────────────────────────────────────────────────────────────────────────

def stage3_score_and_rank(
    csv_path: str,
    symptoms: str,
    output_csv: Optional[str] = None,
    top_k: int = 50,
    llm_provider: str = "openai",
    llm_model: str = "deepseek-v4-flash",
    api_key: Optional[str] = None,
    base_url: Optional[str] = None,
) -> pd.DataFrame:
    """
    Stage 3 独立函数接口。

    静态列直接用内置映射表打分，动态列由 LLM 给出 per-value scores。
    """
    scorer = Stage3Scorer(
        csv_path=csv_path,
        symptoms=symptoms,
        top_k=top_k,
        llm_provider=llm_provider,
        llm_model=llm_model,
        api_key=api_key,
        base_url=base_url,
    )
    result = scorer.run()
    if output_csv:
        scorer.save(result, output_csv)
    return result