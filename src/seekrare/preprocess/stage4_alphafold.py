"""
stage4_alphafold.py — Stage 4 可选模块：AlphaFold3 蛋白结构预测

从 Stage 3 排序结果中选取变异位点，提取其对应的突变蛋白序列，
提交给 AlphaFold3 Server 或 ColabFold 进行结构预测。

核心流程（基于 vcf_to_protein_fasta 逻辑）:
  1. 从 Stage 3 CSV 读取 (CHROM, POS, REF, ALT, gene_name)
  2. 用 GTF 找到该基因的 CDS 转录本区间
  3. 从参考基因组提取 CDS 序列，应用变异
  4. 翻译为突变蛋白序列（AA FASTA）
  5. 调用 AlphaFold3 预测结构

Usage:
    from seekrare.preprocess.stage4_alphafold import run_alphafold_prediction
    results = run_alphafold_prediction(
        stage3_csv="seekrare_output/stage3_ranked.csv",
        ref_fasta="/path/to/GRCh38.fa",
        gtf_file="/path/to/genomic.gtf",
        top_n=5,
        alphafold_mode="server",
        output_dir="alphafold_results",
    )
"""

from __future__ import annotations

import concurrent.futures
import logging
import os
import re
import subprocess
import sys
import urllib.request
import urllib.parse
import urllib.error
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
from loguru import logger

# ─────────────────────────────────────────────────────────────────────────────
# 密码子表
# ─────────────────────────────────────────────────────────────────────────────

CODON2AA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# ─────────────────────────────────────────────────────────────────────────────
# 数据结构
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class TranscriptRecord:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    cds_parts: List[Tuple[int, int]] = field(default_factory=list)


@dataclass
class VcfRecord:
    chrom: str
    pos: int
    ref: str
    alt: str


# ─────────────────────────────────────────────────────────────────────────────
# GTF 解析
# ─────────────────────────────────────────────────────────────────────────────

def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"


def parse_gtf_attributes(raw: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in raw.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def parse_gtf(gtf_path: str) -> Tuple[Dict[str, str], Dict[str, TranscriptRecord]]:
    """
    解析 GTF 文件，返回:
      gene_to_name: gene_id → gene_name
      transcripts: transcript_id → TranscriptRecord
    """
    gene_to_name: Dict[str, str] = {}
    transcripts: Dict[str, TranscriptRecord] = {}

    with open(gtf_path, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            if not raw_line or raw_line.startswith("#"):
                continue
            cols = raw_line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            chrom_raw, _, feature_type, start, end, _, strand, _, attr_text = cols
            chrom = normalize_chrom(chrom_raw)
            attrs = parse_gtf_attributes(attr_text)
            start_i = int(start)
            end_i = int(end)

            if feature_type == "gene":
                gene_id = attrs.get("gene_id", "")
                gene_name = attrs.get("gene_name", "") or attrs.get("gene", "") or gene_id
                if gene_id:
                    gene_to_name[gene_id] = gene_name
                continue

            if feature_type in ("transcript", "mRNA", "mrna"):
                tx_id = attrs.get("transcript_id", "")
                parent = attrs.get("gene_id", "")
                if tx_id and parent:
                    transcripts[tx_id] = TranscriptRecord(
                        transcript_id=tx_id,
                        gene_id=parent,
                        chrom=chrom,
                        strand=strand,
                    )
                continue

            if feature_type == "CDS":
                parent_raw = attrs.get("transcript_id", "")
                if not parent_raw:
                    continue
                gene_id_raw = attrs.get("gene_id", "")
                for tx_id in parent_raw.split(","):
                    tx_id = tx_id.strip()
                    if not tx_id:
                        continue
                    if tx_id not in transcripts:
                        transcripts[tx_id] = TranscriptRecord(
                            transcript_id=tx_id,
                            gene_id=gene_id_raw or "UNKNOWN",
                            chrom=chrom,
                            strand=strand,
                        )
                    transcripts[tx_id].cds_parts.append((start_i, end_i))

    return gene_to_name, transcripts


def build_gene_index(
    transcripts: Dict[str, TranscriptRecord],
    gene_to_name: Dict[str, str],
) -> Dict[str, List[Tuple[int, int, str, str]]]:
    """
    按染色体构建基因区间索引: chrom → [(start, end, gene_id, gene_name)]
    """
    index: Dict[str, List[Tuple[int, int, str, str]]] = defaultdict(list)

    # 按 gene_id 分组
    gene_tx: Dict[str, List[str]] = defaultdict(list)
    for tx in transcripts.values():
        gene_tx[tx.gene_id].append(tx.transcript_id)

    for gene_id, tx_ids in gene_tx.items():
        if not tx_ids:
            continue
        gene_name = gene_to_name.get(gene_id, gene_id)

        # 取该基因所有转录本 CDS 的最大范围
        all_parts: List[Tuple[int, int]] = []
        for tx_id in tx_ids:
            tx = transcripts.get(tx_id)
            if tx:
                all_parts.extend(tx.cds_parts)

        if not all_parts:
            continue

        chrom = transcripts[tx_ids[0]].chrom
        start = min(p[0] for p in all_parts)
        end = max(p[1] for p in all_parts)
        index[chrom].append((start, end, gene_id, gene_name))

    for chrom in index:
        index[chrom].sort(key=lambda x: x[0])

    return index


def find_genes_by_pos(
    intervals: List[Tuple[int, int, str, str]], pos: int
) -> List[Tuple[str, str]]:
    return [(gene_id, gene_name) for start, end, gene_id, gene_name in intervals if start <= pos <= end]


# ─────────────────────────────────────────────────────────────────────────────
# 基因组序列操作
# ─────────────────────────────────────────────────────────────────────────────

def load_fasta(ref_path: str) -> Dict[str, str]:
    genome: Dict[str, List[str]] = defaultdict(list)
    current: Optional[str] = None

    with open(ref_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = normalize_chrom(line[1:].split()[0])
                if current not in genome:
                    genome[current] = []
                continue
            if current is None:
                raise ValueError("FASTA format error: sequence before header.")
            genome[current].append(line.upper())

    return {k: "".join(v) for k, v in genome.items()}


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def extract_ref_seq(genome: Dict[str, str], chrom: str, start: int, end: int) -> Optional[str]:
    seq = genome.get(chrom)
    if seq is None:
        return None
    if start < 1 or end > len(seq):
        return None
    return seq[start - 1:end].upper()


def build_transcript_cds(
    transcript: TranscriptRecord, genome: Dict[str, str]
) -> Tuple[str, List[Tuple[int, int, int]]]:
    """从基因组提取 CDS 序列，构建 (cds_idx, genome_pos, offset) 映射"""
    if transcript.strand == "-":
        parts = sorted(transcript.cds_parts, key=lambda x: x[0], reverse=True)
    else:
        parts = sorted(transcript.cds_parts, key=lambda x: x[0])

    assembled: List[str] = []
    mapping: List[Tuple[int, int, int]] = []
    cds_cursor = 0

    for start, end in parts:
        segment = extract_ref_seq(genome, transcript.chrom, start, end)
        if segment is None:
            continue
        if transcript.strand == "-":
            segment = reverse_complement(segment)
            g_positions = list(range(end, start - 1, -1))
        else:
            g_positions = list(range(start, end + 1))

        assembled.append(segment)
        for idx, gpos in enumerate(g_positions):
            mapping.append((cds_cursor + idx, gpos, idx))
        cds_cursor += len(segment)

    return "".join(assembled), mapping


def apply_variant_to_cds(
    cds_ref: str,
    cds_mapping: List[Tuple[int, int, int]],
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    strand: str,
) -> Optional[Tuple[str, int]]:
    """将 VCF 变异应用到 CDS 序列，返回 (突变后 CDS, 变异在 CDS 起始位置)"""
    gpos_to_cds: Dict[int, int] = {gpos: cds_idx for cds_idx, gpos, _ in cds_mapping}

    start_cds = gpos_to_cds.get(pos)
    if start_cds is None:
        return None

    ref_len = len(ref)
    replaced_positions = [pos + i for i in range(ref_len)]
    if any(gpos not in gpos_to_cds for gpos in replaced_positions):
        return None

    cds_positions = [gpos_to_cds[gpos] for gpos in replaced_positions]
    min_pos = min(cds_positions)
    max_pos = max(cds_positions)

    if (max_pos - min_pos + 1) != ref_len:
        return None

    ref_on_cds = cds_ref[min_pos:max_pos + 1]

    # 检查 REF 匹配
    if strand == "-":
        alt_seq = reverse_complement(alt)
        ref_expected = reverse_complement(ref)
    else:
        alt_seq = alt
        ref_expected = ref

    if ref_on_cds != ref_expected:
        logging.warning(
            "REF mismatch at %s:%d %s>%s, expected CDS ref=%s got=%s",
            chrom, pos, ref, alt, ref_expected, ref_on_cds
        )
        return None

    mutated = cds_ref[:min_pos] + alt_seq + cds_ref[max_pos + 1:]
    return mutated, min_pos


def translate_cds(cds_seq: str) -> str:
    aa: List[str] = []
    usable_len = (len(cds_seq) // 3) * 3
    for i in range(0, usable_len, 3):
        codon = cds_seq[i:i + 3].upper()
        aa.append(CODON2AA.get(codon, "X"))
    return "".join(aa)


def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_fasta(path: str, records: List[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            fh.write(f"{wrap_fasta(seq)}\n")


# ─────────────────────────────────────────────────────────────────────────────
# AlphaFold3 预测
# ─────────────────────────────────────────────────────────────────────────────

class AlphaFold3Client:
    """
    AlphaFold3 结构预测客户端。

    支持两种模式:
      server  — EBI AlphaFold Server (https://alphafold.ebi.ac.uk)
      colabfold — 自部署 ColabFold server

    Parameters
    ----------
    mode : str
        "server" 或 "colabfold"
    base_url : str
        API endpoint
    api_key : str, optional
        ColabFold 需要 API key
    """

    SERVER_URL = "https://alphafold.ebi.ac.uk"

    def __init__(
        self,
        mode: str = "server",
        base_url: Optional[str] = None,
        api_key: Optional[str] = None,
        timeout: int = 120,
    ):
        self.mode = mode
        self.base_url = base_url or self.SERVER_URL
        self.api_key = api_key
        self.timeout = timeout

    def submit_sequence(
        self,
        sequence: str,
        gene_name: str,
        job_name: Optional[str] = None,
    ) -> str:
        """
        提交单个蛋白序列到 AlphaFold3 Server。

        Parameters
        ----------
        sequence : str
            蛋白氨基酸序列
        gene_name : str
            基因名（用于生成 job 名称）
        job_name : str, optional
            自定义 job 名称

        Returns
        -------
        str: job_id 或预测结果
        """
        if self.mode == "server":
            return self._submit_server(sequence, gene_name, job_name)
        else:
            return self._submit_colabfold(sequence, gene_name, job_name)

    def _submit_server(
        self, sequence: str, gene_name: str, job_name: Optional[str] = None
    ) -> str:
        """EBI AlphaFold Server（实际上这里是 AlphaFold2/3 预测入口）"""
        job_name = job_name or f"{gene_name}_variant"

        # AlphaFold Server 的 search API
        url = f"{self.base_url}/search/text/sequence"
        data = urllib.parse.urlencode({
            "sequence": sequence,
            "name": job_name,
        }).encode()

        req = urllib.request.Request(
            url,
            data=data,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            method="POST",
        )

        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                result = resp.read().decode()
                # 解析返回的 job_id（Server 返回 JSON）
                import json
                result_obj = json.loads(result)
                return result_obj.get("jobId", result_obj.get("id", str(result_obj)))
        except Exception as e:
            raise RuntimeError(f"AlphaFold Server submission failed: {e}")

    def _submit_colabfold(
        self, sequence: str, gene_name: str, job_name: Optional[str] = None
    ) -> str:
        """ColabFold self-hosted API"""
        if not self.api_key:
            raise ValueError("ColabFold mode requires api_key")

        job_name = job_name or f"{gene_name}_variant"

        payload = {
            "sequence": sequence,
            "name": job_name,
            "model_preset": "auto",
        }

        data = json.dumps(payload).encode()
        req = urllib.request.Request(
            f"{self.base_url}/predict",
            data=data,
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {self.api_key}",
            },
            method="POST",
        )

        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                result = json.loads(resp.read().decode())
                return result.get("job_id", str(result))
        except Exception as e:
            raise RuntimeError(f"ColabFold submission failed: {e}")

    def poll_result(self, job_id: str, poll_interval: int = 30, max_wait: int = 3600) -> dict:
        """
        轮询 AlphaFold Server 预测结果。

        Parameters
        ----------
        job_id : str
            提交时返回的 job_id
        poll_interval : int
            轮询间隔（秒）
        max_wait : int
            最大等待时间（秒）

        Returns
        -------
        dict: 包含 pdb_url, pdb_content 等字段
        """
        url = f"{self.base_url}/result/{job_id}/details"

        import time
        elapsed = 0
        while elapsed < max_wait:
            try:
                req = urllib.request.Request(url, method="GET")
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    result = json.loads(resp.read().decode())
                    if result.get("status") in ("COMPLETED", "finished", "success"):
                        return result
                    elif result.get("status") in ("FAILED", "error"):
                        raise RuntimeError(f"AlphaFold job failed: {result}")
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    pass  # 还没准备好
                else:
                    raise

            time.sleep(poll_interval)
            elapsed += poll_interval

        raise TimeoutError(f"AlphaFold job {job_id} timed out after {max_wait}s")

    def download_pdb(self, pdb_url: str, output_path: str) -> None:
        """下载 PDB 文件到本地"""
        req = urllib.request.Request(pdb_url, method="GET")
        with urllib.request.urlopen(req, timeout=self.timeout) as resp:
            with open(output_path, "wb") as out:
                out.write(resp.read())
        logger.info(f"  PDB saved → {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 核心：VCF 位点 → 蛋白序列
# ─────────────────────────────────────────────────────────────────────────────

def variant_to_protein(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    gene_name: str,
    transcripts: Dict[str, TranscriptRecord],
    gene_to_name: Dict[str, str],
    gene_index: Dict[str, List[Tuple[int, int, str, str]]],
    genome: Dict[str, str],
) -> Optional[Tuple[str, str, str, str]]:
    """
    将一个 (CHROM, POS, REF, ALT, gene_name) 转换为突变蛋白序列。

    Returns
    -------
    (transcript_id, ref_aa, mut_aa, cds_mut_start) or None
    """
    # 找到这个位置涉及的基因
    candidates = find_genes_by_pos(gene_index.get(chrom, []), pos)

    # 优先匹配 gene_name
    matched_tx = None
    for gid, gname in candidates:
        if gname.upper() == gene_name.upper():
            for tx_id, tx in transcripts.items():
                if tx.gene_id == gid and tx.cds_parts:
                    matched_tx = tx
                    break
            if matched_tx:
                break

    if not matched_tx:
        # Fallback: 用第一个匹配的基因
        for gid, gname in candidates:
            for tx_id, tx in transcripts.items():
                if tx.gene_id == gid and tx.cds_parts:
                    matched_tx = tx
                    break
            if matched_tx:
                break

    if not matched_tx:
        return None

    tx = matched_tx

    # 构建 CDS 序列和映射
    cds_ref, cds_mapping = build_transcript_cds(tx, genome)
    if not cds_ref:
        return None

    ref_aa = translate_cds(cds_ref)

    # 应用变异
    result = apply_variant_to_cds(cds_ref, cds_mapping, chrom, pos, ref, alt, tx.strand)
    if result is None:
        return None

    mut_cds, cds_mut_start = result
    mut_aa = translate_cds(mut_cds)

    if mut_aa == ref_aa:
        return None  # 同义突变

    return tx.transcript_id, ref_aa, mut_aa, str(cds_mut_start)


def read_sites_from_csv(csv_path: str, top_n: Optional[int] = None) -> List[dict]:
    """
    从 Stage 3 CSV 读取位点列表。

    Returns
    -------
    list of dict: {chrom, pos, ref, alt, gene_name}
    """
    df = pd.read_csv(csv_path)
    if top_n:
        df = df.head(top_n)

    sites = []
    for _, row in df.iterrows():
        try:
            chrom = str(row.get("CHROM", row.get("chrom", ""))).strip()
            pos = int(float(row.get("POS", row.get("pos", 0))))
            ref = str(row.get("REF", row.get("ref", ""))).strip()
            alt = str(row.get("ALT", row.get("alt", ""))).strip()
            gene_name = str(row.get("gene_name", row.get("gene", ""))).strip()

            if chrom and pos and ref and alt:
                sites.append({
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gene_name": gene_name or "UNKNOWN",
                })
        except Exception:
            continue

    return sites


def build_transcript_cache(
    transcripts: Dict[str, TranscriptRecord],
    genome: Dict[str, str],
) -> Dict[str, Tuple[str, str]]:
    """
    为所有有 CDS 的转录本构建 (cds_ref, ref_aa) 缓存。

    Returns
    -------
    transcript_id → (cds_seq, ref_aa)
    """
    cache: Dict[str, Tuple[str, str]] = {}
    for tx_id, tx in transcripts.items():
        if not tx.cds_parts:
            continue
        cds_ref, _ = build_transcript_cds(tx, genome)
        if cds_ref:
            cache[tx_id] = (cds_ref, translate_cds(cds_ref))
    return cache


# ─────────────────────────────────────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────────────────────────────────────

def run_alphafold_prediction(
    stage3_csv: str,
    ref_fasta: str,
    gtf_file: str,
    output_dir: str,
    top_n: int = 10,
    alphafold_mode: str = "server",
    alphafold_url: Optional[str] = None,
    alphafold_api_key: Optional[str] = None,
    submit_only: bool = False,
    poll_interval: int = 30,
    max_wait: int = 3600,
) -> pd.DataFrame:
    """
    Stage 4 AlphaFold3 蛋白结构预测。

    从 Stage 3 CSV 取 top-N 位点，提取突变蛋白序列，
    提交 AlphaFold3 预测，返回结果汇总。

    Parameters
    ----------
    stage3_csv : str
        Stage 3 排序后的 CSV 路径
    ref_fasta : str
        参考基因组 FASTA 路径
    gtf_file : str
        GTF 注释文件路径
    output_dir : str
        输出目录
    top_n : int
        取前 N 个位点
    alphafold_mode : str
        "server"（EBI AlphaFold）或 "colabfold"（自部署）
    alphafold_url : str, optional
        自定义 AlphaFold API URL
    alphafold_api_key : str, optional
        ColabFold API key
    submit_only : bool
        若 True，只提交不等待结果
    poll_interval : int
        轮询结果间隔（秒）
    max_wait : int
        最大等待时间（秒）

    Returns
    -------
    pd.DataFrame: 预测结果汇总表
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Stage 4: AlphaFold3 Protein Structure Prediction")
    logger.info("=" * 60)

    # 1. 加载资源
    logger.info(f"加载参考基因组: {ref_fasta}")
    genome = load_fasta(ref_fasta)
    logger.info(f"  %d 条染色体", len(genome))

    logger.info(f"解析 GTF: {gtf_file}")
    gene_to_name, transcripts = parse_gtf(gtf_file)
    logger.info(f"  %d 基因, %d 转录本", len(gene_to_name), len(transcripts))

    gene_index = build_gene_index(transcripts, gene_to_name)
    tx_cache = build_transcript_cache(transcripts, genome)
    logger.info(f"  %d 转录本有 CDS 序列", len(tx_cache))

    # 2. 读取位点
    sites = read_sites_from_csv(stage3_csv, top_n=top_n)
    logger.info(f"从 Stage 3 读取 %d 个位点", len(sites))

    # 3. AlphaFold 客户端
    af_client = AlphaFold3Client(
        mode=alphafold_mode,
        base_url=alphafold_url,
        api_key=alphafold_api_key,
    )

    # 4. 对每个位点提取蛋白序列并提交
    fasta_records: List[Tuple[str, str]] = []
    results: List[dict] = []

    for site in sites:
        chrom = normalize_chrom(site["chrom"])
        pos = site["pos"]
        ref = site["ref"].upper()
        alt = site["alt"].upper()
        gene_name = site["gene_name"]

        info = variant_to_protein(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            gene_name=gene_name,
            transcripts=transcripts,
            gene_to_name=gene_to_name,
            gene_index=gene_index,
            genome=genome,
        )

        if info is None:
            logger.warning(f"  跳过 (无法提取蛋白): {chrom}:{pos} {ref}>{alt} @ {gene_name}")
            results.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene_name": gene_name,
                "status": "no_protein",
                "transcript_id": None,
                "ref_aa": None,
                "mut_aa": None,
                "af_job_id": None,
                "af_pdb_url": None,
                "is_reasonable_structure": None,
            })
            continue

        tx_id, ref_aa, mut_aa, cds_start = info

        # 写入 FASTA
        header = (
            f"{chrom}:{pos}:{ref}>{alt}|{gene_name}|tx={tx_id}"
            f"|CDSpos={cds_start}|len={len(mut_aa)}"
        )
        fasta_records.append((header, mut_aa))

        # 提交 AlphaFold
        job_id = None
        try:
            job_id = af_client.submit_sequence(
                sequence=mut_aa,
                gene_name=f"{gene_name}_{chrom}_{pos}",
                job_name=f"{gene_name}_{chrom}_{pos}_{ref}_{alt}",
            )
            logger.info(f"  提交: {chrom}:{pos} {ref}>{alt} @ {gene_name} → job={job_id}")
        except Exception as e:
            logger.error(f"  提交失败: {chrom}:{pos} {ref}>{alt} → {e}")

        results.append({
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "gene_name": gene_name,
            "status": "submitted" if job_id else "submit_failed",
            "transcript_id": tx_id,
            "ref_aa": ref_aa,
            "mut_aa": mut_aa,
            "af_job_id": job_id,
            "af_pdb_url": None,
            "is_reasonable_structure": None,
        })

    # 保存蛋白 FASTA
    fasta_path = output_dir / "mutant_proteins.fa"
    write_fasta(str(fasta_path), fasta_records)
    logger.info(f"  蛋白 FASTA → {fasta_path}")

    # 保存汇总 CSV
    result_df = pd.DataFrame(results)
    csv_path = output_dir / "alphafold_summary.csv"
    result_df.to_csv(str(csv_path), index=False)
    logger.info(f"  结果汇总 → {csv_path}")

    # 如果需要轮询结果
    if not submit_only:
        for i, row in result_df[result_df["status"] == "submitted"].iterrows():
            job_id = row["af_job_id"]
            if not job_id:
                continue
            try:
                result = af_client.poll_result(job_id, poll_interval, max_wait)
                result_df.at[i, "af_pdb_url"] = result.get("pdb_url")
                result_df.at[i, "is_reasonable_structure"] = result.get("is_reasonable_structure")
            except Exception as e:
                logger.error(f"  轮询失败 job={job_id}: {e}")

        result_df.to_csv(str(csv_path), index=False)

    logger.info(f"Stage 4 AlphaFold 完成: {output_dir}")
    return result_df


# ─────────────────────────────────────────────────────────────────────────────
# 快捷函数（供 pipeline.py 调用）
# ─────────────────────────────────────────────────────────────────────────────

def stage4_alphafold_prediction(
    stage3_csv: str,
    ref_fasta: str,
    gtf_file: str,
    output_dir: str,
    top_n: int = 10,
    **kwargs,
) -> pd.DataFrame:
    """
    Stage 4 AlphaFold wrapper（供 pipeline.py 调用）。
    """
    return run_alphafold_prediction(
        stage3_csv=stage3_csv,
        ref_fasta=ref_fasta,
        gtf_file=gtf_file,
        output_dir=output_dir,
        top_n=top_n,
        **kwargs,
    )