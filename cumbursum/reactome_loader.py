#!/usr/bin/env python3
"""Reactome pathway loader for CumBurSum.

Parses Reactome flat files and returns gene-symbol -> pathway mappings,
using a UniProt idmapping file (shared with the parent CASCADE pipeline).

Inputs (expected in cumbursum/data/, downloadable via data/download.sh):
    - UniProt2Reactome_All_Levels.txt  (uniprot \t pathway_id \t url \t name \t evidence \t species)
    - ReactomePathways.txt             (pathway_id \t name \t species)
    - ReactomePathwaysRelation.txt     (parent_id \t child_id)  [for hierarchy rollup]

Output:
    {pathway_id: {'name': str, 'genes': set[str], 'size': int}}

Filters to Homo sapiens only. Optionally skips MT-* genes (recommended
for reference-alignment pipelines where mtDNA is aligned to a different
reference than nuclear DNA).
"""
from __future__ import annotations

import gzip
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Set

logger = logging.getLogger(__name__)

HUMAN_SPECIES = "Homo sapiens"
DEFAULT_DATA_DIR = Path(__file__).parent / "data"


def _resolve_uniprot_idmap() -> Path:
    """Resolve the UniProt idmapping file path, in order:

      1. `CUMBURSUM_UNIPROT_IDMAP` environment variable
      2. `cumbursum/data/HUMAN_9606_idmapping.dat.gz`
      3. `~/conservation_data/HUMAN_9606_idmapping.dat.gz` (CASCADE default)

    Downloaders: see `cumbursum/data/README.md`.
    """
    import os
    env = os.environ.get("CUMBURSUM_UNIPROT_IDMAP")
    if env:
        return Path(env)
    local = DEFAULT_DATA_DIR / "HUMAN_9606_idmapping.dat.gz"
    if local.exists():
        return local
    return Path.home() / "conservation_data" / "HUMAN_9606_idmapping.dat.gz"


DEFAULT_UNIPROT_IDMAP = _resolve_uniprot_idmap()

# Swiss-Prot canonical ID pattern (6-char reviewed IDs).
_SWISSPROT_CANONICAL = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


def _load_uniprot_to_gene(
    idmap_path: Path = DEFAULT_UNIPROT_IDMAP,
) -> Dict[str, str]:
    """Minimal in-pathways loader for UniProt -> Gene_Name.

    Inlined here (instead of importing analyzers.uniprot_mapper) to avoid
    the analyzers package's transitive imports. Format is the standard
    UniProt idmapping 3-col TSV (uniprot_id \t db_type \t db_id).
    """
    if not idmap_path.exists():
        raise FileNotFoundError(f"UniProt idmapping file not found: {idmap_path}")

    up_to_gene: Dict[str, str] = {}
    with gzip.open(idmap_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            uniprot_id, db_type, db_id = parts
            if db_type != "Gene_Name":
                continue
            up_to_gene[uniprot_id] = db_id
    return up_to_gene


class ReactomeLoader:
    """Load Reactome human pathways keyed by gene symbol."""

    def __init__(
        self,
        data_dir: Optional[Path] = None,
        uniprot_idmap: Optional[Path] = None,
        exclude_mt: bool = True,
        rollup_hierarchy: bool = False,
    ):
        self.data_dir = Path(data_dir) if data_dir else DEFAULT_DATA_DIR
        self.uniprot_idmap = Path(uniprot_idmap) if uniprot_idmap else DEFAULT_UNIPROT_IDMAP
        self.exclude_mt = exclude_mt
        self.rollup_hierarchy = rollup_hierarchy

        # Output structures (populated by load())
        self.pathway_to_genes: Dict[str, Dict] = {}
        self.gene_to_pathways: Dict[str, Set[str]] = defaultdict(set)
        self.parent_of: Dict[str, Set[str]] = defaultdict(set)  # child -> parents
        self.children_of: Dict[str, Set[str]] = defaultdict(set)  # parent -> children
        self._loaded = False

    def load(self) -> Dict[str, Dict]:
        """Populate pathway_to_genes and gene_to_pathways. Returns pathway_to_genes."""
        uniprot_file = self.data_dir / "UniProt2Reactome_All_Levels.txt"
        pathways_file = self.data_dir / "ReactomePathways.txt"
        if not uniprot_file.exists():
            raise FileNotFoundError(f"Reactome UniProt mapping not found: {uniprot_file}")
        if not pathways_file.exists():
            raise FileNotFoundError(f"Reactome pathways not found: {pathways_file}")

        up_to_gene = _load_uniprot_to_gene(self.uniprot_idmap)

        # Parse the pathway-name table (fallback names, plus species check).
        pathway_names: Dict[str, str] = {}
        with pathways_file.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                pid, name, species = parts[0], parts[1], parts[2]
                if species != HUMAN_SPECIES:
                    continue
                pathway_names[pid] = name.strip()

        # Parse the UniProt -> pathway mapping, filtered to human.
        pathway_genes: Dict[str, Set[str]] = defaultdict(set)
        n_rows = 0
        n_human = 0
        n_mapped = 0
        n_unmapped_uniprot: Set[str] = set()

        with uniprot_file.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                n_rows += 1
                uniprot_id, pid, _url, name, _evidence, species = parts[:6]
                if species != HUMAN_SPECIES:
                    continue
                n_human += 1

                # Strip any isoform suffix ("P12345-2" -> "P12345").
                base_id = uniprot_id.split("-")[0]
                gene = up_to_gene.get(base_id) or up_to_gene.get(uniprot_id)
                if not gene:
                    n_unmapped_uniprot.add(base_id)
                    continue

                if self.exclude_mt and gene.startswith("MT-"):
                    continue

                pathway_genes[pid].add(gene)
                # Prefer the name from the mapping row (always present) but fall
                # back on the table if missing.
                if pid not in pathway_names and name:
                    pathway_names[pid] = name.strip()
                n_mapped += 1

        # Optional: hierarchical rollup — parent pathways inherit all
        # descendants' genes. This lets "Mitochondrial translation" (parent)
        # pool hits that otherwise sit in "Mitochondrial translation
        # elongation" (child) alone. Load the relation file if requested.
        if self.rollup_hierarchy:
            relation_file = self.data_dir / "ReactomePathwaysRelation.txt"
            if relation_file.exists():
                human_relations = 0
                with relation_file.open("r", encoding="utf-8", errors="replace") as f:
                    for line in f:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 2:
                            continue
                        parent, child = parts[0], parts[1]
                        # Only keep human relations.
                        if not (parent.startswith("R-HSA-") and child.startswith("R-HSA-")):
                            continue
                        self.parent_of[child].add(parent)
                        self.children_of[parent].add(child)
                        human_relations += 1
                logger.info("Loaded %d human pathway parent-child relations.", human_relations)

                # Compute transitive closure: for each pathway, gather all
                # descendants' genes. Iterate until stable (should converge
                # in a few passes — Reactome depth is typically ≤7).
                for _ in range(10):
                    changed = False
                    for parent, children in self.children_of.items():
                        parent_genes = pathway_genes.setdefault(parent, set())
                        before = len(parent_genes)
                        for child in children:
                            child_genes = pathway_genes.get(child, set())
                            parent_genes.update(child_genes)
                        if len(parent_genes) > before:
                            changed = True
                    if not changed:
                        break
            else:
                logger.warning(
                    "Hierarchy rollup requested but %s not found. Skipping.",
                    relation_file,
                )

        # Build the final structure. Keep only pathways with >=1 gene.
        self.pathway_to_genes = {
            pid: {
                "name": pathway_names.get(pid, pid),
                "genes": genes,
                "size": len(genes),
            }
            for pid, genes in pathway_genes.items()
            if genes
        }

        # Gene-to-pathways inverted index.
        self.gene_to_pathways = defaultdict(set)
        for pid, rec in self.pathway_to_genes.items():
            for g in rec["genes"]:
                self.gene_to_pathways[g].add(pid)

        self._loaded = True
        logger.info(
            "Reactome loaded: %d pathways, %d unique genes, %d/%d human rows mapped (%d unknown UniProts)",
            len(self.pathway_to_genes),
            len(self.gene_to_pathways),
            n_mapped,
            n_human,
            len(n_unmapped_uniprot),
        )
        return self.pathway_to_genes

    # Convenience -----------------------------------------------------------

    def pathways_for_gene(self, gene: str) -> Set[str]:
        if not self._loaded:
            self.load()
        return self.gene_to_pathways.get(gene, set())

    def genes_in_pathway(self, pathway_id: str) -> Set[str]:
        if not self._loaded:
            self.load()
        rec = self.pathway_to_genes.get(pathway_id)
        return rec["genes"] if rec else set()

    def name(self, pathway_id: str) -> str:
        if not self._loaded:
            self.load()
        rec = self.pathway_to_genes.get(pathway_id)
        return rec["name"] if rec else pathway_id


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    loader = ReactomeLoader()
    pwdb = loader.load()
    print(f"{len(pwdb):,} human Reactome pathways loaded.")
    # Quick sanity print: some known-relevant pathway sizes.
    probes = [
        "Collagen biosynthesis and modifying enzymes",
        "Extracellular matrix organization",
        "The citric acid (TCA) cycle and respiratory electron transport",
    ]
    name_to_id = {rec["name"]: pid for pid, rec in pwdb.items()}
    for p in probes:
        pid = name_to_id.get(p)
        if pid:
            print(f"  {p}: {pid} ({pwdb[pid]['size']} genes)")
        else:
            print(f"  {p}: NOT FOUND")
