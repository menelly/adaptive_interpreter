#!/usr/bin/env python3
"""Per-gene constraint weights from gnomAD v4.1.

Reads gnomad.v4.1.constraint_metrics.tsv and exposes a gene -> weight
lookup that biases burden scoring toward genes that are *intolerant to
variation* (high LOF pLI, low missense oe).

Weighting model:
    weight(gene) = 1.0 baseline
                 + 0.5  if lof.pLI >= 0.9  (strongly haploinsufficient)
                 + 0.25 if lof.pLI >= 0.5
                 - 0.25 if mis.oe  >  1.2  (missense-tolerant)
                 ...

Clamped to [0.5, 2.0].

The goal is modest multiplicative amplification/attenuation, not a
dominating factor — we're saying "a variant in a strongly-constrained
gene contributes more to pathway burden than an equivalent variant in a
tolerant gene," which is biologically well-established (Karczewski et al.
2020, Samocha et al. 2014).
"""
from __future__ import annotations

import csv
import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)

DEFAULT_CONSTRAINT_TSV = Path(__file__).parent / "data" / "gnomad_constraint.v4.1.tsv"
CACHE_JSON = Path(__file__).parent / "data" / "gnomad_constraint_cache.json"


def _to_float(x: str) -> Optional[float]:
    if x is None or x == "" or x == "NA":
        return None
    try:
        return float(x)
    except ValueError:
        return None


def build_constraint_cache(
    tsv_path: Path = DEFAULT_CONSTRAINT_TSV,
    cache_path: Path = CACHE_JSON,
    force: bool = False,
) -> Dict[str, Dict[str, float]]:
    """Parse the gnomAD constraint TSV once and cache as compact JSON.

    Picks the canonical/MANE-select transcript per gene. If multiple
    flagged canonical, the first wins (arbitrary but stable).
    """
    if cache_path.exists() and not force:
        with cache_path.open("r") as f:
            return json.load(f)

    if not tsv_path.exists():
        raise FileNotFoundError(
            f"gnomAD constraint TSV not found: {tsv_path}. "
            "Download from https://storage.googleapis.com/gcp-public-data--gnomad/"
            "release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
        )

    per_gene: Dict[str, Dict[str, float]] = {}
    with tsv_path.open("r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = (row.get("gene") or "").strip()
            if not gene:
                continue
            canonical = (row.get("canonical") or "").lower() == "true"
            mane = (row.get("mane_select") or "").lower() == "true"
            if not (canonical or mane):
                # Prefer canonical/MANE. For genes with none flagged, keep
                # first-seen; canonical/MANE wins if seen later.
                if gene in per_gene:
                    continue
            lof_pli = _to_float(row.get("lof.pLI"))
            lof_oe = _to_float(row.get("lof.oe"))
            mis_oe = _to_float(row.get("mis.oe"))
            if lof_pli is None and lof_oe is None and mis_oe is None:
                continue
            per_gene[gene] = {
                "lof_pLI": lof_pli if lof_pli is not None else 0.0,
                "lof_oe": lof_oe if lof_oe is not None else 1.0,
                "mis_oe": mis_oe if mis_oe is not None else 1.0,
            }

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w") as f:
        json.dump(per_gene, f, separators=(",", ":"))
    logger.info("Built constraint cache: %d genes -> %s", len(per_gene), cache_path)
    return per_gene


def weight_for_gene(
    gene: str,
    constraint_data: Dict[str, Dict[str, float]],
) -> float:
    """Compute multiplicative burden weight for a gene.

    Defaults to 1.0 if the gene is missing from constraint data (e.g.,
    aliases, novel). Clamped to [0.5, 2.0] to prevent outlier domination.
    """
    rec = constraint_data.get(gene)
    if not rec:
        return 1.0
    w = 1.0
    pli = rec.get("lof_pLI", 0.0)
    mis_oe = rec.get("mis_oe", 1.0)
    if pli >= 0.9:
        w += 0.5
    elif pli >= 0.5:
        w += 0.25
    if mis_oe > 1.2:
        w -= 0.25
    if mis_oe < 0.6:
        w += 0.25
    return max(0.5, min(2.0, w))


class ConstraintWeights:
    """Lazy-loading wrapper."""

    def __init__(
        self,
        tsv_path: Optional[Path] = None,
        cache_path: Optional[Path] = None,
    ):
        self.tsv_path = Path(tsv_path) if tsv_path else DEFAULT_CONSTRAINT_TSV
        self.cache_path = Path(cache_path) if cache_path else CACHE_JSON
        self._data: Optional[Dict[str, Dict[str, float]]] = None

    def _ensure(self) -> Dict[str, Dict[str, float]]:
        if self._data is None:
            self._data = build_constraint_cache(self.tsv_path, self.cache_path)
        return self._data

    def weight(self, gene: str) -> float:
        return weight_for_gene(gene, self._ensure())

    def info(self, gene: str) -> Optional[Dict[str, float]]:
        return self._ensure().get(gene)


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    force = "--rebuild" in sys.argv
    data = build_constraint_cache(force=force)
    print(f"{len(data):,} genes in constraint cache")
    probes = ["ATP5F1A","DLD","SLC25A5","GFM1","AMPD1","G6PD","ALDH2","KCNMA1",
              "COL5A2","NPHS2","MYO7A","CACNA1I","DYSF","BMPR2","TFG","SEC22A"]
    print(f"\n{'gene':10s} {'pLI':>5s} {'mis_oe':>7s} {'weight':>7s}")
    cw = ConstraintWeights()
    for g in probes:
        rec = data.get(g, {})
        w = cw.weight(g)
        print(f"  {g:10s} {rec.get('lof_pLI',0):5.2f} {rec.get('mis_oe',1.0):7.2f} {w:7.2f}")
