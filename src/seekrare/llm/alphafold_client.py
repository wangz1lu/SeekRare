"""
AlphaFold Server integration.

Calls AlphaFold2 server API to predict protein structures for candidate genes.
Useful for validating candidate variants affecting protein structure.

Usage:
    from seekrare.alphafold import AlphaFoldClient
    af = AlphaFoldClient()
    result = af.predict_sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH...")
    af.download_pdb(result["job_id"], "output.pdb")
"""

from __future__ import annotations

import os
import time
import urllib.request
import urllib.parse
import json
from pathlib import Path
from typing import Optional, Union
from loguru import logger


class AlphaFoldClient:
    """
    AlphaFold2 Server API client.

    Uses the public AlphaFold server:
    https://alphafold.ebi.ac.uk/search/text/sequence

    Note: For high-throughput use, consider AlphaFold ColabFold
    (https://github.com/sokrypton/ColabFold) which can be self-hosted.
    """

    BASE_URL = "https://alphafold.ebi.ac.uk"

    def __init__(self, output_dir: Union[str, Path] = "alphafold_results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def predict_sequence(
        self,
        sequence: str,
        gene_name: str = "unknown",
        wait: bool = True,
        poll_interval: int = 30,
        max_wait: int = 3600,
    ) -> dict:
        """
        Submit a protein sequence for AlphaFold prediction.

        Parameters
        ----------
        sequence : str
            Protein sequence (amino acids)
        gene_name : str
            Gene name for naming output files
        wait : bool
            If True, poll until job is complete
        poll_interval : int
            Seconds between status checks
        max_wait : int
            Maximum seconds to wait

        Returns
        -------
        dict
            {"job_id": str, "status": str, "pdb_url": str, "ae_url": str}
        """
        logger.info(f"Submitting AlphaFold prediction for {gene_name} ({len(sequence)} AA)")

        # Submit job
        data = urllib.parse.urlencode({"query": sequence}).encode()
        req = urllib.request.Request(
            f"{self.BASE_URL}/search/text/sequence",
            data=data,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            method="POST",
        )

        with urllib.request.urlopen(req, timeout=60) as resp:
            result = json.loads(resp.read().decode())
            job_id = result.get("jobId", "")
            logger.info(f"  AlphaFold job submitted: {job_id}")

        status_url = f"{self.BASE_URL}/search/job/{job_id}"
        elapsed = 0

        if wait:
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
                        logger.info(f"  AlphaFold complete! PDB: {pdb_url}")
                        return {
                            "job_id": job_id,
                            "status": status,
                            "pdb_url": pdb_url,
                            "ae_url": ae_url,
                        }
                    elif status == "ERROR":
                        logger.error("  AlphaFold job failed")
                        return {"job_id": job_id, "status": status}

                except Exception as e:
                    logger.warning(f"  Status check failed: {e}")

        return {"job_id": job_id, "status": "PENDING"}

    def download_pdb(self, pdb_url: str, output_path: Union[str, Path]) -> Path:
        """Download PDB file from AlphaFold result URL."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Downloading PDB: {pdb_url}")
        urllib.request.urlretrieve(pdb_url, output_path)
        logger.info(f"  Saved: {output_path}")
        return output_path

    def predict_and_download(
        self,
        sequence: str,
        gene_name: str,
        wait: bool = True,
    ) -> Optional[Path]:
        """Submit, wait, and download PDB in one call."""
        result = self.predict_sequence(sequence, gene_name, wait=wait)

        if result.get("status") == "SUCCESS" and result.get("pdb_url"):
            out_path = self.output_dir / f"{gene_name}_af2.pdb"
            return self.download_pdb(result["pdb_url"], out_path)

        logger.warning(f"AlphaFold did not complete for {gene_name}")
        return None


# ── ColabFold (self-hosted alternative) ────────────────────────────────────────

class ColabFoldClient:
    """
    ColabFold API client for self-hosted high-throughput structure prediction.

    Uses the ColabFold API (https://github.com/sokrypton/ColabFold).
    Requires a running ColabFold server.
    """

    def __init__(self, base_url: str = "http://localhost:8000", output_dir: str = "colabfold_results"):
        self.base_url = base_url.rstrip("/")
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def predict(
        self,
        sequence: str,
        gene_name: str = "unknown",
        num_recycles: int = 3,
        model_preset: str = "auto",
    ) -> dict:
        """Submit sequence to ColabFold server."""
        import requests

        payload = {
            "sequence": sequence,
            "model_preset": model_preset,
            "num_recycles": num_recycles,
        }

        logger.info(f"Submitting ColabFold for {gene_name}")
        resp = requests.post(f"{self.base_url}/predict", json=payload, timeout=30)
        resp.raise_for_status()
        result = resp.json()

        job_id = result.get("job_id", "")
        logger.info(f"  ColabFold job: {job_id}")

        return result

    def download_result(self, job_id: str, output_path: Union[str, Path]) -> Path:
        """Download ColabFold result PDB."""
        import requests

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        resp = requests.get(f"{self.base_url}/result/{job_id}/pdb", timeout=60)
        resp.raise_for_status()

        with open(output_path, "w") as f:
            f.write(resp.text)

        logger.info(f"  Saved: {output_path}")
        return output_path
