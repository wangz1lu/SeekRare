"""
AlphaFold3 Server annotation — Stage 2 optional annotation module.

对候选基因的蛋白质序列进行结构预测，评估变异对蛋白结构的影响。
支持 AlphaFold3 public server API 和自部署版本。

AlphaFold3 Server: https://alphafold.ebi.ac.uk
ColabFold (self-hosted): https://github.com/sokrypton/ColabFold

Usage:
    from seekrare.annotation import AlphaFold3Annotator

    annotator = AlphaFold3Annotator(
        mode="server",           # "server" (EBI) or "colabfold"
        base_url="https://alphafold.ebi.ac.uk",
        output_dir="alphafold_results",
    )
    result = annotator.predict_sequence(sequence, gene_name)
    annotator.download_pdb(result["pdb_url"], "output.pdb")
"""

from __future__ import annotations

import json
import time
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


class AlphaFold3Annotator:
    """
    AlphaFold3 结构预测 annotation。

    Mode "server":  使用 EBI AlphaFold2 server (免费，无需 API key)
                     AlphaFold3 目前仅通过 AlphaFold Server 提供：
                     https://alphafold.ebi.ac.uk/search/text/sequence

    Mode "colabfold": 自部署 ColabFold (需要 API key 和 self-hosted server)

    Parameters
    ----------
    mode : str
        "server" (EBI AlphaFold2) 或 "colabfold" (自部署)
    base_url : str
        API base URL
    output_dir : str or Path
        PDB 输出目录
    api_key : str, optional
        API key（ColabFold 自部署需要）
    """

    def __init__(
        self,
        mode: str = "server",
        base_url: str = "https://alphafold.ebi.ac.uk",
        output_dir: Union[str, Path] = "alphafold_results",
        api_key: Optional[str] = None,
    ):
        self.mode = mode.lower()
        self.base_url = base_url.rstrip("/")
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.api_key = api_key

    def predict_sequence(
        self,
        sequence: str,
        gene_name: str = "unknown",
        wait: bool = True,
        poll_interval: int = 30,
        max_wait: int = 3600,
    ) -> dict:
        """
        提交蛋白质序列进行结构预测。

        Returns
        -------
        dict
            {"job_id": str, "status": str, "pdb_url": str, "ae_url": str}
        """
        logger.info(f"AlphaFold3: submitting {gene_name} ({len(sequence)} AA)")

        if self.mode == "server":
            return self._predict_server(sequence, gene_name, wait, poll_interval, max_wait)
        else:
            return self._predict_colabfold(sequence, gene_name)

    def _predict_server(
        self,
        sequence: str,
        gene_name: str,
        wait: bool,
        poll_interval: int,
        max_wait: int,
    ) -> dict:
        """AlphaFold2 EBI server submission."""
        data = urllib.parse.urlencode({"query": sequence}).encode()
        req = urllib.request.Request(
            f"{self.base_url}/search/text/sequence",
            data=data,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            method="POST",
        )

        with urllib.request.urlopen(req, timeout=60) as resp:
            result = json.loads(resp.read().decode())
            job_id = result.get("jobId", "")
            logger.info(f"  Job submitted: {job_id}")

        if not wait:
            return {"job_id": job_id, "status": "PENDING"}

        # Poll until complete
        status_url = f"{self.base_url}/search/job/{job_id}"
        elapsed = 0

        while elapsed < max_wait:
            time.sleep(poll_interval)
            elapsed += poll_interval

            try:
                with urllib.request.urlopen(status_url, timeout=30) as sresp:
                    status_data = json.loads(sresp.read().decode())
                    status = status_data.get("jobStatus", "RUNNING")

                logger.info(f"  Status: {status} ({elapsed}s)")

                if status == "SUCCESS":
                    pdb_url = status_data.get("pdbUrl", "")
                    ae_url = status_data.get("aeUrl", "")
                    return {"job_id": job_id, "status": status, "pdb_url": pdb_url, "ae_url": ae_url}
                elif status == "ERROR":
                    return {"job_id": job_id, "status": "ERROR"}

            except Exception as e:
                logger.warning(f"  Status check failed: {e}")

        return {"job_id": job_id, "status": "TIMEOUT"}

    def _predict_colabfold(self, sequence: str, gene_name: str) -> dict:
        """ColabFold self-hosted API."""
        import requests

        payload = {"query": sequence, "model_preset": "auto", "num_recycles": 3}
        resp = requests.post(
            f"{self.base_url}/predict",
            json=payload,
            headers={"Authorization": f"Bearer {self.api_key}"},
            timeout=30,
        )
        resp.raise_for_status()
        return resp.json()

    def download_pdb(self, pdb_url: str, output_path: Union[str, Path]) -> Path:
        """下载 PDB 文件."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(pdb_url, output_path)
        logger.info(f"AlphaFold3: PDB saved → {output_path}")
        return output_path

    def predict_and_download(
        self,
        sequence: str,
        gene_name: str,
        wait: bool = True,
    ) -> Optional[Path]:
        """一步完成: 提交 + 等待 + 下载."""
        result = self.predict_sequence(sequence, gene_name, wait=wait)
        if result.get("status") == "SUCCESS" and result.get("pdb_url"):
            out_path = self.output_dir / f"{gene_name}_af2.pdb"
            return self.download_pdb(result["pdb_url"], out_path)
        return None

    def annotate_genes(
        self,
        df: pd.DataFrame,
        protein_sequence_col: str = "protein_sequence",
        gene_col: str = "gene_name",
    ) -> pd.DataFrame:
        """
        对 DataFrame 中的基因进行 AlphaFold 注释（批量模式）。

        注意: 实际结构预测需要调用外部 API，这里仅记录待预测的基因列表。

        Parameters
        ----------
        df : pd.DataFrame
            变异 DataFrame
        protein_sequence_col : str
            蛋白序列列名（如果有）
        gene_col : str
            基因名列

        Returns
        -------
        pd.DataFrame
            添加了 alphafold_* 列的 DataFrame
        """
        df = df.copy()

        # 记录需要预测的基因
        genes_to_predict = df[gene_col].dropna().unique().tolist()

        df["alphafold_predicted"] = False
        df["alphafold_pdb_url"] = None
        df["alphafold_status"] = "not_submitted"

        logger.info(f"AlphaFold3: {len(genes_to_predict)} 个基因待预测结构")
        logger.info("AlphaFold3: 建议使用 predict_and_download() 批量提交")
        logger.info("AlphaFold3: AlphaFold2 server 目前免费，但有速率限制")

        return df


# ── Stage 2 extension stub: Genos annotation ──────────────────────────────────


class GenosAnnotationStub:
    """
    Genos 模型 annotation — Stage 2 optional annotation module.

    对变异进行 Genos 模型（临床遗传学专用 LLM）的自动化注释。
    该模块为预留接口，具体功能取决于 Genos 模型 API 的实际规格。

    Genos 模型能力（预留）:
    - 变异致病性评分
    - ACMG criteria 自动评估
    - 蛋白结构影响预测
    - 文献证据综合

    Usage (预留):
        from seekrare.annotation import GenosAnnotation

        genos = GenosAnnotation(api_key="...")
        df = genos.annotate_variants(df, patient_phenotype="...")
    """

    def __init__(
        self,
        api_key: Optional[str] = None,
        base_url: str = "https://api.genos.tech/v1",
        model: str = "genos-clinical-v1",
    ):
        self.api_key = api_key
        self.base_url = base_url
        self.model = model
        logger.warning(
            "GenosAnnotation is a STUB. "
            "API endpoint and spec to be confirmed. "
            "Do not use in production until finalized."
        )

    def annotate_variants(
        self,
        df: pd.DataFrame,
        patient_phenotype: str = "",
        **kwargs,
    ) -> pd.DataFrame:
        """
        对变异 DataFrame 进行 Genos 模型注释（预留）.

        Parameters
        ----------
        df : pd.DataFrame
            变异 DataFrame
        patient_phenotype : str
            患者表型描述

        Returns
        -------
        pd.DataFrame
            添加了 genos_* 列的 DataFrame（目前仅为 stub）
        """
        df = df.copy()
        df["genos_pathogenicity"] = None
        df["genos_acmg_criteria"] = None
        df["genos_evidence"] = None
        df["genos_confidence"] = None
        logger.warning("GenosAnnotation.annotate_variants() called but is a STUB — no actual annotation performed")
        return df
