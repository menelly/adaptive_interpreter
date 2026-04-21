#!/usr/bin/env python3
"""Cross-database consensus aggregator for CumBurSum.

Premise: a biological theme is *real* if it surfaces as elevated in more
than one independent pathway database. The 50 MSigDB Hallmarks and the
2,800+ Reactome pathways (with hierarchical rollup) were curated by
different teams using different criteria — so when they BOTH flag a theme,
the signal is reinforced against database-specific artifact.

This module:
  1. Runs burden + permutation in each database independently
  2. Identifies each database's top-ranked pathways
  3. Groups them by the genes driving them
  4. Flags "consensus themes" — gene-sets that drive top pathways in >1 DB

Output: a list of ConsensusTheme objects suitable for the report header.
"""
from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple

from cumbursum.burden_scorer import MECHANISMS, PathwayBurden
from cumbursum.null_distribution import BurdenSignificance

logger = logging.getLogger(__name__)


@dataclass
class DatabaseHit:
    db_name: str       # "Reactome" or "Hallmarks"
    pathway_id: str
    pathway_name: str
    pathway_size: int
    genes_hit: Set[str]
    burden: float
    z: float
    p_value: float
    mechanism: str


@dataclass
class ConsensusTheme:
    """A biological theme surfacing in ≥1 database.

    Named by the best-fit pathway (shortest name among the contributing
    pathways with best p). `supporting_dbs` tells you how many DBs backed
    it; `driver_genes` is the union of variant-hit genes across all
    contributing pathways.
    """
    theme_name: str
    driver_genes: Set[str] = field(default_factory=set)
    supporting_hits: List[DatabaseHit] = field(default_factory=list)
    convergence_score: float = 0.0  # sum of driver adj_scores (set post-clustering)

    @property
    def n_dbs(self) -> int:
        return len({h.db_name for h in self.supporting_hits})

    @property
    def best_p(self) -> float:
        return min(h.p_value for h in self.supporting_hits) if self.supporting_hits else 1.0

    @property
    def is_consensus(self) -> bool:
        return self.n_dbs >= 2


def top_hits_from_run(
    db_name: str,
    burdens: Dict[str, PathwayBurden],
    sig: Dict[str, Dict[str, BurdenSignificance]],
    top_n: int = 25,
    p_ceiling: float = 0.05,
    mechanism: str = "total",
) -> List[DatabaseHit]:
    """Extract top-ranked pathways from one database's run."""
    hits: List[DatabaseHit] = []
    for pid, rec in sig.items():
        s = rec.get(mechanism)
        if not s:
            continue
        if s.p_value > p_ceiling:
            continue
        b = burdens[pid]
        hits.append(DatabaseHit(
            db_name=db_name,
            pathway_id=pid,
            pathway_name=b.pathway_name,
            pathway_size=b.pathway_size,
            genes_hit=set(b.genes_hit),
            burden=s.observed,
            z=s.z,
            p_value=s.p_value,
            mechanism=mechanism,
        ))
    hits.sort(key=lambda h: h.p_value)
    return hits[:top_n]


def group_hits_into_themes(
    all_hits: List[DatabaseHit],
    min_overlap: float = 0.5,
) -> List[ConsensusTheme]:
    """Group hits into themes by driving-gene overlap.

    Two hits belong to the same theme if their driver-gene sets have
    Jaccard overlap >= min_overlap.

    This is a greedy single-link clustering, which is imperfect but
    adequate for v0 — and keeps the theme-grouping transparent for the
    paper. (More sophisticated: agglomerative with ward linkage, but not
    justified for a set this small.)
    """
    themes: List[ConsensusTheme] = []
    for hit in all_hits:
        placed = False
        for theme in themes:
            shared = theme.driver_genes & hit.genes_hit
            union = theme.driver_genes | hit.genes_hit
            jaccard = len(shared) / len(union) if union else 0.0
            if jaccard >= min_overlap:
                theme.driver_genes.update(hit.genes_hit)
                theme.supporting_hits.append(hit)
                placed = True
                break
        if not placed:
            themes.append(ConsensusTheme(
                theme_name=hit.pathway_name,
                driver_genes=set(hit.genes_hit),
                supporting_hits=[hit],
            ))

    # After clustering, pick each theme's best-fit name.
    # Priority order:
    #   1. If a Hallmark hit exists in the theme, use its name (Hallmarks are
    #      broad, curated, biologically interpretable).
    #   2. Else, use the Reactome hit with the best p-value AND smallest
    #      pathway size (most specific).
    for theme in themes:
        hallmark_hits = [h for h in theme.supporting_hits if h.db_name == "Hallmarks"]
        if hallmark_hits:
            best_hallmark = min(hallmark_hits, key=lambda h: h.p_value)
            theme.theme_name = best_hallmark.pathway_name
        else:
            best_p = min(h.p_value for h in theme.supporting_hits)
            candidates = [h for h in theme.supporting_hits if h.p_value == best_p]
            best = min(candidates, key=lambda h: h.pathway_size)
            theme.theme_name = best.pathway_name

    # Sort: consensus themes first, then by best p.
    themes.sort(key=lambda t: (not t.is_consensus, t.best_p))
    return themes


def annotate_convergence_scores(
    themes: List[ConsensusTheme],
    all_burdens_by_db: Dict[str, Dict[str, PathwayBurden]],
) -> None:
    """Compute convergence_score per theme = sum of max adj_score for each
    driver gene (max across any pathway, any DB, where that gene contributes).
    This lets a theme built from 3 LP variants clearly outscore a theme
    built from 1 VUS — independent of pathway-size-corrected p-values."""
    # Build gene -> max adj_score lookup across all observed burdens
    gene_to_max_score: Dict[str, float] = {}
    for db_burdens in all_burdens_by_db.values():
        for b in db_burdens.values():
            for v in b.variants:
                prev = gene_to_max_score.get(v.gene, 0.0)
                if v.adj_score > prev:
                    gene_to_max_score[v.gene] = v.adj_score
    for theme in themes:
        theme.convergence_score = sum(
            gene_to_max_score.get(g, 0.0) for g in theme.driver_genes
        )


def build_consensus(
    reactome_burdens: Dict[str, PathwayBurden],
    reactome_sig: Dict[str, Dict[str, BurdenSignificance]],
    hallmark_burdens: Dict[str, PathwayBurden],
    hallmark_sig: Dict[str, Dict[str, BurdenSignificance]],
    top_n_per_db: int = 25,
    p_ceiling_reactome: float = 0.01,  # strict: ~3000 tests
    p_ceiling_hallmarks: float = 0.075,  # moderate: ~50 tests
    mechanism: str = "total",
    min_gene_overlap: float = 0.5,
) -> List[ConsensusTheme]:
    """High-level entry point. Returns themes ranked by consensus and p.

    Different p_ceilings per database reflect the multiple-testing burden:
      - Reactome (~3000 tested pathways): 0.01 raw ≈ ~30 expected false
        positives before BH. Tight to keep theme noise down.
      - Hallmarks (~50 tests): 0.075 raw ≈ ~4 expected false positives.
        Broader because Hallmarks is a much smaller, curated-broad
        hypothesis space with independent signal.
    """
    hits = []
    hits.extend(top_hits_from_run(
        "Reactome", reactome_burdens, reactome_sig,
        top_n=top_n_per_db, p_ceiling=p_ceiling_reactome, mechanism=mechanism,
    ))
    hits.extend(top_hits_from_run(
        "Hallmarks", hallmark_burdens, hallmark_sig,
        top_n=top_n_per_db, p_ceiling=p_ceiling_hallmarks, mechanism=mechanism,
    ))
    logger.info(
        "Consensus input: %d Reactome (p<%s) + %d Hallmarks (p<%s) hits on %s burden",
        sum(1 for h in hits if h.db_name == "Reactome"), p_ceiling_reactome,
        sum(1 for h in hits if h.db_name == "Hallmarks"), p_ceiling_hallmarks,
        mechanism,
    )
    themes = group_hits_into_themes(hits, min_overlap=min_gene_overlap)
    annotate_convergence_scores(themes, {
        "Reactome": reactome_burdens,
        "Hallmarks": hallmark_burdens,
    })
    # Re-sort: consensus first, then by (descending convergence, ascending p)
    themes.sort(key=lambda t: (not t.is_consensus, -t.convergence_score, t.best_p))
    n_consensus = sum(1 for t in themes if t.is_consensus)
    logger.info(
        "Consensus output: %d themes total, %d consensus (≥2 DBs)",
        len(themes), n_consensus,
    )
    return themes
