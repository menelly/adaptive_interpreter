#!/usr/bin/env python3
"""MSigDB Hallmarks loader for CumBurSum.

Parses a .gmt file (tab-separated: name, url, gene1, gene2, ...) and
exposes the same interface as ReactomeLoader:
    - pathway_to_genes: {hallmark_id: {'name': str, 'genes': set[str], 'size': int}}
    - gene_to_pathways: {gene: set[hallmark_id]}

Hallmarks are broader and curated (50 sets, ~36-200 genes each) — good
orthogonal signal to Reactome's reaction-granular pathways.
"""
from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Set

logger = logging.getLogger(__name__)

DEFAULT_GMT = Path(__file__).parent / "data" / "h.all.v2023.2.Hs.symbols.gmt"


class HallmarksLoader:
    """Load MSigDB Hallmark gene sets."""

    def __init__(
        self,
        gmt_path: Optional[Path] = None,
        exclude_mt: bool = True,
    ):
        self.gmt_path = Path(gmt_path) if gmt_path else DEFAULT_GMT
        self.exclude_mt = exclude_mt
        self.pathway_to_genes: Dict[str, Dict] = {}
        self.gene_to_pathways: Dict[str, Set[str]] = defaultdict(set)
        self._loaded = False

    def load(self) -> Dict[str, Dict]:
        if not self.gmt_path.exists():
            raise FileNotFoundError(f"Hallmarks .gmt file not found: {self.gmt_path}")
        pw: Dict[str, Dict] = {}
        with self.gmt_path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                name, _url, *genes = parts
                # Pretty name: strip HALLMARK_ prefix, replace underscores.
                pretty = name.replace("HALLMARK_", "").replace("_", " ").title()
                pretty = f"Hallmark: {pretty}"
                gene_set: Set[str] = set()
                for g in genes:
                    g = g.strip()
                    if not g:
                        continue
                    if self.exclude_mt and g.startswith("MT-"):
                        continue
                    gene_set.add(g)
                pw[name] = {
                    "name": pretty,
                    "genes": gene_set,
                    "size": len(gene_set),
                }
        self.pathway_to_genes = pw
        self.gene_to_pathways = defaultdict(set)
        for pid, rec in pw.items():
            for g in rec["genes"]:
                self.gene_to_pathways[g].add(pid)
        self._loaded = True
        logger.info(
            "Hallmarks loaded: %d gene sets, %d unique genes",
            len(pw), len(self.gene_to_pathways),
        )
        return pw

    def pathways_for_gene(self, gene: str) -> Set[str]:
        if not self._loaded:
            self.load()
        return self.gene_to_pathways.get(gene, set())

    def name(self, pid: str) -> str:
        if not self._loaded:
            self.load()
        rec = self.pathway_to_genes.get(pid)
        return rec["name"] if rec else pid


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    L = HallmarksLoader()
    pw = L.load()
    print(f"{len(pw)} Hallmarks loaded.")
    for pid in ["HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MYOGENESIS",
                "HALLMARK_COAGULATION", "HALLMARK_NOTCH_SIGNALING"]:
        rec = pw.get(pid)
        if rec:
            print(f"  {rec['name']}: {rec['size']} genes")
