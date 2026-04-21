#!/usr/bin/env python3
"""Within-individual permutation null for CumBurSum.

Asks: "Is this pathway more loaded than a random pathway of the same size
would be for THIS person's variants?"

Method:
  For each pathway P of size K (|P|=K genes in pathway):
    For i in 1..N_perm:
      Draw K random genes from the universe of Reactome-covered genes.
      Compute burden score using the person's variants restricted to those
      K genes (only genes that have variants contribute; this is the same
      operation as scoring a real pathway).
    Empirical p = (# perms with burden >= observed + 1) / (N + 1).

Why this design:
  - Random sampling from the Reactome gene universe (not from the person's
    variant-hit genes) means every pathway is tested against "a pathway
    of equivalent size drawn from the same annotation universe." The
    person's variant distribution is the constant; pathway membership
    is what's shuffled.
  - Size-matched: sampling K genes per perm controls for pathway size bias.
  - Per-mechanism: each mechanism (total/dn/lof/gof) gets its own p-value.
  - The "+1" in both numerator and denominator avoids p=0 and encodes
    "we've only run N perms, so the smallest defensible p is 1/(N+1)".

Benjamini-Hochberg FDR correction across the pathway set.
"""
from __future__ import annotations

import logging
import random
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Set, Tuple

from cumbursum.burden_scorer import MECHANISMS, PathwayBurden, Variant

logger = logging.getLogger(__name__)

DEFAULT_N_PERM = 1000


@dataclass
class BurdenSignificance:
    """Significance results for one pathway, one mechanism."""
    pathway_id: str
    mechanism: str
    observed: float
    mean_null: float
    std_null: float
    z: float
    p_value: float
    p_adj: float = 1.0  # populated by benjamini_hochberg


def _sample_burden_for_random_geneset(
    gene_sample: Set[str],
    gene_to_variants: Dict[str, List[Variant]],
    mechanism: str,
) -> float:
    """Sum the given mechanism's adj_* over all variants whose gene is in the sample."""
    attr = {"total": "adj_score", "dn": "adj_dn", "lof": "adj_lof", "gof": "adj_gof"}[mechanism]
    total = 0.0
    hits = gene_sample & gene_to_variants.keys()
    for g in hits:
        for v in gene_to_variants[g]:
            total += getattr(v, attr)
    return total


def permutation_significance(
    burdens: Dict[str, PathwayBurden],
    variants: Iterable[Variant],
    gene_universe: List[str],
    mechanisms: Tuple[str, ...] = MECHANISMS,
    n_perm: int = DEFAULT_N_PERM,
    random_seed: Optional[int] = 42,
) -> Dict[str, Dict[str, BurdenSignificance]]:
    """Compute per-pathway, per-mechanism empirical p-values.

    Returns {pathway_id: {mechanism: BurdenSignificance}}.

    gene_universe should be the full list of genes known to Reactome (or
    more narrowly, those with any pathway annotation) — this is what we
    sample random "pseudo-pathways" from.
    """
    if random_seed is not None:
        random.seed(random_seed)

    # Build the per-gene variant index once.
    gene_to_variants: Dict[str, List[Variant]] = {}
    for v in variants:
        gene_to_variants.setdefault(v.gene, []).append(v)

    universe = list(gene_universe)
    universe_size = len(universe)
    logger.info("Permutation null: %d pathways × %d mechanisms × %d perms (universe=%d)",
                len(burdens), len(mechanisms), n_perm, universe_size)

    results: Dict[str, Dict[str, BurdenSignificance]] = {}

    # Group pathways by size so we only draw each size-specific null once.
    pathways_by_size: Dict[int, List[PathwayBurden]] = {}
    for b in burdens.values():
        pathways_by_size.setdefault(b.pathway_size, []).append(b)

    # Compute nulls once per (size, mechanism), reuse across all pathways of that size.
    for size, group in pathways_by_size.items():
        null_by_mech: Dict[str, List[float]] = {m: [] for m in mechanisms}
        for _ in range(n_perm):
            sample = set(random.sample(universe, k=min(size, universe_size)))
            for m in mechanisms:
                null_by_mech[m].append(
                    _sample_burden_for_random_geneset(sample, gene_to_variants, m)
                )

        # Precompute mean/std for each mechanism at this size.
        null_stats = {}
        for m, vals in null_by_mech.items():
            mean = sum(vals) / len(vals)
            var = sum((x - mean) ** 2 for x in vals) / max(len(vals) - 1, 1)
            std = var ** 0.5
            null_stats[m] = (mean, std, sorted(vals))

        for b in group:
            rec: Dict[str, BurdenSignificance] = {}
            for m in mechanisms:
                observed = {"total": b.total, "dn": b.dn, "lof": b.lof, "gof": b.gof}[m]
                mean, std, sorted_vals = null_stats[m]
                # Empirical p (right-tail, with +1 correction).
                n_ge = sum(1 for x in sorted_vals if x >= observed)
                p = (n_ge + 1) / (n_perm + 1)
                z = (observed - mean) / std if std > 0 else 0.0
                rec[m] = BurdenSignificance(
                    pathway_id=b.pathway_id,
                    mechanism=m,
                    observed=observed,
                    mean_null=mean,
                    std_null=std,
                    z=z,
                    p_value=p,
                )
            results[b.pathway_id] = rec

    return results


def benjamini_hochberg(
    sig_results: Dict[str, Dict[str, BurdenSignificance]],
    mechanisms: Tuple[str, ...] = MECHANISMS,
) -> None:
    """In-place BH-FDR adjustment, per mechanism, across all pathways."""
    for m in mechanisms:
        # Collect (pid, sig) for this mechanism.
        items = [(pid, rec[m]) for pid, rec in sig_results.items() if m in rec]
        # Sort ascending by raw p.
        items.sort(key=lambda it: it[1].p_value)
        n = len(items)
        if n == 0:
            continue
        # Standard BH: q_i = p_i * n / rank
        prev_q = 1.0
        # Walk from largest rank down to enforce monotonicity.
        adj_by_pid: Dict[str, float] = {}
        for rank in range(n, 0, -1):
            pid, sig = items[rank - 1]
            q = sig.p_value * n / rank
            if q < prev_q:
                prev_q = q
            adj_by_pid[pid] = min(prev_q, 1.0)
        for pid, q in adj_by_pid.items():
            sig_results[pid][m].p_adj = q


if __name__ == "__main__":
    import argparse
    from pathlib import Path
    from cumbursum.reactome_loader import ReactomeLoader
    from cumbursum.burden_scorer import (
        compute_pathway_burdens, load_cascade_variants, rank_burdens,
    )

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cascade", required=True)
    ap.add_argument("--min-score", type=float, default=0.0)
    ap.add_argument("--n-perm", type=int, default=DEFAULT_N_PERM)
    ap.add_argument("--mechanism", choices=list(MECHANISMS), default="total")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    loader = ReactomeLoader()
    pwdb = loader.load()
    variants = load_cascade_variants(Path(args.cascade), min_score=args.min_score)
    burdens = compute_pathway_burdens(variants, pwdb)
    universe = sorted(loader.gene_to_pathways.keys())
    sig = permutation_significance(
        burdens, variants, universe,
        mechanisms=MECHANISMS, n_perm=args.n_perm,
    )
    benjamini_hochberg(sig)

    # Print top by BH-adjusted p on the chosen mechanism.
    rows = []
    for pid, rec in sig.items():
        r = rec[args.mechanism]
        b = burdens[pid]
        rows.append((r.p_adj, r.p_value, r.z, r.observed, b))
    rows.sort(key=lambda t: (t[0], t[1]))
    print(f"\nTop {args.top} pathways by BH-adjusted p ({args.mechanism} burden):\n")
    print(f"{'p_adj':>8}  {'p':>7}  {'z':>6}  {'obs':>7}  {'hits':>4}/{'size':<4}  pathway")
    for p_adj, p_raw, z, obs, b in rows[:args.top]:
        print(f"{p_adj:>8.4f}  {p_raw:>7.4f}  {z:>6.2f}  {obs:>7.3f}  {b.n_genes_hit:>4}/{b.pathway_size:<4}  {b.pathway_name}")
