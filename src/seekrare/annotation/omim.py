"""
OMIM (Online Mendelian Inheritance in Man) disease-gene annotation.

Loads OMIM gene map /morbidmap and annotates variants with:
- Associated diseases per gene
- Inheritance pattern (AD/AR/XLR/Y-linked)
- MIM number
- Phenotype description

Usage:
    from seekrare.annotation import OMIMLoader
    omim = OMIMLoader("/ref/omim/genemap2.txt")
    annotated = omim.annotate_genes(df)   # df needs gene_name column
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import pandas as pd
from loguru import logger


class OMIMLoader:
    """
    OMIM gene-disease annotation loader.

    Supports:
    - genemap2.txt (gene map, tab-separated)
    - mim2gene.txt (MIM number ↔ gene mapping)
    - morbidmap.txt (phenotype map, inherited disorders)

    Parameters
    ----------
    omim_dir : str or Path
        Directory containing OMIM data files
    """

    def __init__(self, omim_dir: Union[str, Path]):
        self.omim_dir = Path(omim_dir)
        self._gene_map: dict[str, list[dict]] = {}
        self._mim_to_gene: dict[str, str] = {}
        self._phenotype_map: dict[str, list[dict]] = {}

    def load(self):
        """Load all OMIM data files."""
        self._load_genemap()
        self._load_mim2gene()
        self._load_morbidmap()
        n_genes = len(self._gene_map)
        n_phen = len(self._phenotype_map)
        logger.info(f"OMIM loaded: {n_genes} genes, {n_phen} phenotypes")

    def _load_genemap(self):
        """Load genemap2.txt."""
        gm_path = self.omim_dir / "genemap2.txt"
        if not gm_path.exists():
            logger.warning(f"genemap2.txt not found at {gm_path}")
            return

        try:
            df = pd.read_csv(gm_path, sep="\t", low_memory=False, dtype=str)
        except Exception as e:
            logger.error(f"Failed to load genemap2.txt: {e}")
            return

        for _, row in df.iterrows():
            approved_sym = str(row.get("Approved Symbol", "")).strip()
            if not approved_sym or approved_sym == "nan":
                continue

            mim_num = str(row.get("MIM Number", "")).strip()
            cyto = str(row.get("Cyto Location", "")).strip()
            gene_types = str(row.get("GeneIDS", "")).strip()
            comments = str(row.get("Comments", "")).strip()

            entries = self._gene_map.setdefault(approved_sym, [])
            entries.append({
                "mim_number": mim_num if mim_num != "nan" else "",
                "cyto_location": cyto if cyto != "nan" else "",
                "gene_ids": gene_types if gene_types != "nan" else "",
                "comments": comments if comments != "nan" else "",
                "type": "gene",
            })

    def _load_mim2gene(self):
        """Load mim2gene.txt."""
        m2g_path = self.omim_dir / "mim2gene.txt"
        if not m2g_path.exists():
            return

        try:
            df = pd.read_csv(m2g_path, sep="\t", low_memory=False, dtype=str)
        except Exception:
            return

        for _, row in df.iterrows():
            mim = str(row.get("MIM Number", "")).strip()
            gene_sym = str(row.get("Approved Symbol", "")).strip()
            if gene_sym and gene_sym != "nan":
                self._mim_to_gene[mim] = gene_sym

    def _load_morbidmap(self):
        """Load morbidmap.txt (phenotype → inheritance)."""
        mm_path = self.omim_dir / "morbidmap.txt"
        if not mm_path.exists():
            logger.warning(f"morbidmap.txt not found at {mm_path}")
            return

        try:
            df = pd.read_csv(mm_path, sep="\t", low_memory=False, dtype=str)
        except Exception as e:
            logger.error(f"Failed to load morbidmap.txt: {e}")
            return

        for _, row in df.iterrows():
            pheno = str(row.get("Phenotype", "")).strip()
            genes_str = str(row.get("MIM Number", "")).strip()
            inheritance = str(row.get("Inheritance", "")).strip()

            if not genes_str or genes_str == "nan":
                continue

            # morbidmap format: "MIM Number" is actually gene MIM
            gene_mim = genes_str.split(",")[0].strip()
            gene_sym = self._mim_to_gene.get(gene_mim, "")

            key = gene_sym.upper() if gene_sym else gene_mim
            if key:
                self._phenotype_map.setdefault(key, [])
                self._phenotype_map[key].append({
                    "phenotype": pheno if pheno != "nan" else "",
                    "gene_mim": gene_mim,
                    "inheritance": inheritance if inheritance != "nan" else "",
                })

    def annotate_genes(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Add OMIM disease annotations to a DataFrame with gene_name column.

        Parameters
        ----------
        df : pd.DataFrame
            Must have a gene_name (or gene_symbol) column

        Returns
        -------
        pd.DataFrame
            df with new OMIM columns
        """
        if not self._gene_map and not self._phenotype_map:
            self.load()

        df = df.copy()
        gene_col = "gene_name" if "gene_name" in df.columns else "gene_symbol"

        for col in ["omim_diseases", "omim_inheritance", "omim_mim_number"]:
            if col not in df.columns:
                df[col] = None

        for idx, row in df.iterrows():
            gene = str(row.get(gene_col, "")).strip().upper()
            if not gene or gene in ("NAN", ""):
                continue

            diseases = []
            inheritances = set()
            mim_nums = []

            # From gene map
            for entry in self._gene_map.get(gene, []):
                if entry["mim_number"]:
                    mim_nums.append(entry["mim_number"])

            # From phenotype map
            for entry in self._phenotype_map.get(gene, []):
                if entry["phenotype"]:
                    diseases.append(entry["phenotype"])
                if entry["inheritance"]:
                    inheritances.add(entry["inheritance"])

            if diseases:
                df.at[idx, "omim_diseases"] = "; ".join(diseases)
            if inheritances:
                df.at[idx, "omim_inheritance"] = "|".join(sorted(inheritances))
            if mim_nums:
                df.at[idx, "omim_mim_number"] = "|".join(mim_nums[:3])

        n_annotated = df["omim_diseases"].notna().sum()
        logger.info(f"OMIM: {n_annotated}/{len(df)} variants annotated")
        return df

    @staticmethod
    def infer_inheritance_from_phenotype(phenotype: str) -> str:
        """Infer inheritance pattern from phenotype string."""
        p = phenotype.lower()
        if any(k in p for k in ["autosomal dominant", "ad", "haploinsufficiency"]):
            return "AD"
        elif any(k in p for k in ["autosomal recessive", "ar", "homozygous"]):
            return "AR"
        elif "x-linked" in p:
            return "XL"
        return "Unknown"
