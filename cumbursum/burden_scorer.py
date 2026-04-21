#!/usr/bin/env python3
"""Pathway burden scorer for CumBurSum.

Given a CASCADE results TSV (per-variant mechanism-adjusted scores) and a
Reactome pathway database (pathway_id -> {name, genes}), compute per-pathway
burden for each mechanism:

    total_burden = Σ adj_score   over variants whose gene is in the pathway
    dn_burden    = Σ adj_dn      ''
    lof_burden   = Σ adj_lof     ''
    gof_burden   = Σ adj_gof     ''

Design notes:
  - Burden is a simple sum. The *significance* (permutation) is what makes
    a high burden meaningful; raw sums are just counts × intensity.
  - A variant's contribution to a pathway is the variant's own adj_* score,
    regardless of how many pathways the gene belongs to. Shared genes DO
    contribute to multiple pathways — this is intentional (shared insult),
    but means we need overlap-aware FDR later.
  - MT-* genes dropped before scoring (Dante alignment gotcha; see handoff).
  - Non-numeric adj_* values (empty strings) coerce to 0.
"""
from __future__ import annotations

import csv
import json
import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Pattern, Set, Tuple

logger = logging.getLogger(__name__)

MECHANISMS = ("total", "dn", "lof", "gof")
SCORE_COLUMNS = {
    "total": "adj_score",
    "dn": "adj_dn",
    "lof": "adj_lof",
    "gof": "adj_gof",
}
MIN_PATHWAY_SIZE = 5    # smaller pathways are too unstable for permutation stats
MAX_PATHWAY_SIZE = 500  # huge "root" pathways aren't informative

DEFAULT_NOISE_PATTERNS_PATH = (
    Path(__file__).parent.parent / "utils" / "clinical_noise_patterns.json"
)
# Burden-scorer additional patterns: same spirit as the shared file but
# slightly broader (mucins-with-suffix, PRAMEF-variants, etc.) because
# pathway burden is more sensitive to these than per-variant reporting.
_EXTRA_NOISE_PATTERNS = [
    r"^MUC\d",        # MUC3A, MUC17, MUC19, MUC20 — polymorphic mucins (tandem-repeat noise)
    r"^KRT\d",        # keratin paralogs (large family, many unreliable calls)
    r"^USP17L",       # USP17L1..40 — paralog cluster
    r"^NBPF",         # NBPF1..26 neuroblastoma-breakpoint paralogs
]


def _compile_noise_patterns(
    patterns_path: Optional[Path] = None,
    include_extras: bool = True,
) -> List[Pattern]:
    patterns: List[str] = []
    path = Path(patterns_path) if patterns_path else DEFAULT_NOISE_PATTERNS_PATH
    if path.exists():
        try:
            data = json.loads(path.read_text())
            patterns.extend(p["pattern"] for p in data.get("patterns", []) if "pattern" in p)
        except Exception as e:
            logger.warning("Failed to parse noise patterns file %s: %s", path, e)
    if include_extras:
        patterns.extend(_EXTRA_NOISE_PATTERNS)
    return [re.compile(p) for p in patterns]


def _is_noise_gene(gene: str, compiled_patterns: List[Pattern]) -> bool:
    return any(p.search(gene) for p in compiled_patterns)


@dataclass
class Variant:
    """A single CASCADE-scored variant, reduced to what burden needs."""
    gene: str
    adj_score: float
    adj_dn: float
    adj_lof: float
    adj_gof: float
    hgvs: str = ""
    consequence: str = ""
    gnomad_freq: float = 0.0
    classification: str = ""  # CASCADE final_classification / adj_classification, if available


@dataclass
class PathwayBurden:
    pathway_id: str
    pathway_name: str
    pathway_size: int
    genes_hit: Set[str] = field(default_factory=set)
    variants: List[Variant] = field(default_factory=list)
    total: float = 0.0
    dn: float = 0.0
    lof: float = 0.0
    gof: float = 0.0

    @property
    def n_variants(self) -> int:
        return len(self.variants)

    @property
    def n_genes_hit(self) -> int:
        return len(self.genes_hit)

    def as_row(self) -> Dict:
        return {
            "pathway_id": self.pathway_id,
            "pathway_name": self.pathway_name,
            "pathway_size": self.pathway_size,
            "n_genes_hit": self.n_genes_hit,
            "n_variants": self.n_variants,
            "total_burden": round(self.total, 4),
            "dn_burden": round(self.dn, 4),
            "lof_burden": round(self.lof, 4),
            "gof_burden": round(self.gof, 4),
            "genes_hit": ",".join(sorted(self.genes_hit)),
        }


def _to_float(x: str) -> float:
    if x is None:
        return 0.0
    s = str(x).strip()
    if not s or s.lower() in ("nan", "na", "none"):
        return 0.0
    try:
        return float(s)
    except ValueError:
        return 0.0


def load_cascade_variants(
    tsv_path: Path,
    exclude_mt: bool = True,
    exclude_noise_families: bool = True,
    noise_patterns_path: Optional[Path] = None,
    min_score: float = 0.0,
) -> List[Variant]:
    """Load a CASCADE results TSV into Variant objects.

    min_score: drop variants whose adj_score is below this. 0.0 keeps all
    (including benigns). Use ~0.3 to focus on meaningful signal.
    exclude_noise_families: drop variants in known-noise gene families
    (mucins, olfactory receptors, HLA, etc.) — these are otherwise tandem-
    repeat / highly polymorphic regions that over-inflate pathway burden.
    """
    tsv_path = Path(tsv_path)
    noise_patterns = _compile_noise_patterns(noise_patterns_path) if exclude_noise_families else []

    variants: List[Variant] = []
    skipped_mt = skipped_noise = skipped_low = 0
    noise_hits_by_family: Dict[str, int] = defaultdict(int)

    with tsv_path.open("r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = (row.get("gene") or "").strip()
            if not gene:
                continue
            if exclude_mt and gene.startswith("MT-"):
                skipped_mt += 1
                continue
            if noise_patterns and _is_noise_gene(gene, noise_patterns):
                skipped_noise += 1
                noise_hits_by_family[gene] += 1
                continue
            adj_score = _to_float(row.get("adj_score"))
            if adj_score < min_score:
                skipped_low += 1
                continue
            variants.append(Variant(
                gene=gene,
                adj_score=adj_score,
                adj_dn=_to_float(row.get("adj_dn")),
                adj_lof=_to_float(row.get("adj_lof")),
                adj_gof=_to_float(row.get("adj_gof")),
                hgvs=(row.get("hgvs") or "").strip(),
                consequence=(row.get("molecular_consequence") or "").strip(),
                gnomad_freq=_to_float(row.get("gnomad_freq")),
                classification=(row.get("adj_classification") or
                                row.get("final_classification") or "").strip(),
            ))
    logger.info(
        "Loaded %d variants from %s (exclude_mt=%s, exclude_noise=%s, min_score=%s; "
        "skipped MT=%d, noise=%d, low=%d)",
        len(variants), tsv_path.name, exclude_mt, exclude_noise_families, min_score,
        skipped_mt, skipped_noise, skipped_low,
    )
    if noise_hits_by_family:
        top_noise = sorted(noise_hits_by_family.items(), key=lambda t: -t[1])[:5]
        logger.info("  Top noise-family genes filtered: %s", top_noise)
    return variants


def compute_pathway_burdens(
    variants: Iterable[Variant],
    pathway_to_genes: Dict[str, Dict],
    min_pathway_size: int = MIN_PATHWAY_SIZE,
    max_pathway_size: int = MAX_PATHWAY_SIZE,
) -> Dict[str, PathwayBurden]:
    """For each pathway with >=1 variant hit, aggregate burden.

    Returns {pathway_id: PathwayBurden}. Pathways with no hits are omitted.
    Pathways outside the size window are also omitted.
    """
    # Build gene -> [variants] index so each pathway is one set-intersection.
    gene_to_variants: Dict[str, List[Variant]] = defaultdict(list)
    for v in variants:
        gene_to_variants[v.gene].append(v)

    burdens: Dict[str, PathwayBurden] = {}
    for pid, rec in pathway_to_genes.items():
        pw_size = rec["size"]
        if pw_size < min_pathway_size or pw_size > max_pathway_size:
            continue
        genes = rec["genes"]
        hit_genes = genes & gene_to_variants.keys()
        if not hit_genes:
            continue

        b = PathwayBurden(
            pathway_id=pid,
            pathway_name=rec["name"],
            pathway_size=pw_size,
            genes_hit=set(hit_genes),
        )
        for g in hit_genes:
            for v in gene_to_variants[g]:
                b.variants.append(v)
                b.total += v.adj_score
                b.dn += v.adj_dn
                b.lof += v.adj_lof
                b.gof += v.adj_gof
        burdens[pid] = b

    logger.info("Computed burden for %d pathways (size %d-%d)",
                len(burdens), min_pathway_size, max_pathway_size)
    return burdens


def rank_burdens(
    burdens: Dict[str, PathwayBurden],
    mechanism: str = "total",
    top_n: Optional[int] = None,
) -> List[PathwayBurden]:
    """Return burdens sorted descending by the given mechanism's score."""
    if mechanism not in MECHANISMS:
        raise ValueError(f"mechanism must be one of {MECHANISMS}")
    key = {"total": "total", "dn": "dn", "lof": "lof", "gof": "gof"}[mechanism]
    ranked = sorted(burdens.values(), key=lambda b: getattr(b, key), reverse=True)
    return ranked[:top_n] if top_n else ranked


if __name__ == "__main__":
    import argparse
    from cumbursum.reactome_loader import ReactomeLoader

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cascade", required=True, help="CASCADE results TSV")
    ap.add_argument("--min-score", type=float, default=0.0)
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    loader = ReactomeLoader()
    pwdb = loader.load()
    variants = load_cascade_variants(Path(args.cascade), min_score=args.min_score)
    burdens = compute_pathway_burdens(variants, pwdb)
    print(f"\nTop {args.top} pathways by total burden:\n")
    print(f"{'burden':>8}  {'dn':>6}  {'lof':>6}  {'hits':>4}/{'size':<4}  pathway")
    for b in rank_burdens(burdens, "total", top_n=args.top):
        print(f"{b.total:>8.3f}  {b.dn:>6.3f}  {b.lof:>6.3f}  {b.n_genes_hit:>4}/{b.pathway_size:<4}  {b.pathway_name}")
