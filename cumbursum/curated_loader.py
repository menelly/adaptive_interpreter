#!/usr/bin/env python3
"""Curated-variant loader for CumBurSum.

Reads a simple 4-column TSV of already-curated variants:
    GENE  VARIANT  SYSTEMS  PRIORITY

and produces Variant objects with consequence-based burden scores so that
frameshifts, stop-gains, and canonical splice disruptions get LP-level
weights (not the low scores the raw CASCADE pipeline assigns because it
only reasons about missense mechanisms).

Consequence-to-burden mapping (defaults):
    frameshift (fs)       -> adj_lof = 1.00, adj_score = 1.00
    stop_gain (nonsense)  -> adj_lof = 1.00, adj_score = 1.00
    start_loss            -> adj_lof = 0.95, adj_score = 0.95
    splice, canonical ±1/±2 -> adj_lof = 0.90, adj_score = 0.90
    splice, deeper intron -> adj_lof = 0.60, adj_score = 0.60
    inframe indel         -> adj_lof = 0.50, adj_score = 0.50 (LOF-ish default)
    missense              -> scored from PRIORITY column:
                              High   -> adj_score = 0.70  (mix of DN + LOF)
                              Medium -> adj_score = 0.50
                              Low    -> adj_score = 0.30
                              (unknown) -> 0.30

For the missense rows, contribution is split across DN/LOF:
    adj_dn  = adj_score * 0.6
    adj_lof = adj_score * 0.4
because curated missense of "this matters" origin is more often DN-ish
than pure LOF. This is a blunt prior — fine for v0 sanity checking, but
CASCADE-scored missense values should replace these when available.
"""
from __future__ import annotations

import csv
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from cumbursum.burden_scorer import Variant

logger = logging.getLogger(__name__)

CONSEQUENCE_WEIGHTS = {
    "frameshift":   {"score": 1.00, "lof_frac": 1.00, "dn_frac": 0.00},
    "stop_gain":    {"score": 1.00, "lof_frac": 1.00, "dn_frac": 0.00},
    "start_loss":   {"score": 0.95, "lof_frac": 1.00, "dn_frac": 0.00},
    "splice_canonical": {"score": 0.90, "lof_frac": 1.00, "dn_frac": 0.00},
    "splice_deep":  {"score": 0.60, "lof_frac": 1.00, "dn_frac": 0.00},
    "inframe":      {"score": 0.50, "lof_frac": 0.70, "dn_frac": 0.30},
    "missense_high":   {"score": 0.70, "lof_frac": 0.40, "dn_frac": 0.60},
    "missense_medium": {"score": 0.50, "lof_frac": 0.40, "dn_frac": 0.60},
    "missense_low":    {"score": 0.30, "lof_frac": 0.40, "dn_frac": 0.60},
    "unknown":      {"score": 0.30, "lof_frac": 0.40, "dn_frac": 0.60},
}

RE_STOP_GAIN = re.compile(r"(p\.[A-Za-z]{3}\d+(?:\*|Ter))", re.IGNORECASE)
RE_FRAMESHIFT = re.compile(r"(fs(?:\*\d+)?|fs\b|frameshift)", re.IGNORECASE)
RE_START_LOSS = re.compile(r"p\.Met1", re.IGNORECASE)
RE_SPLICE_CANONICAL = re.compile(r"c\.-?\d+[+\-][12]\b")
RE_SPLICE_DEEP = re.compile(r"c\.-?\d+[+\-][3-9]|c\.-?\d+[+\-]\d{2}")
RE_INFRAME_INDEL = re.compile(r"(del|dup|ins|_\d+)", re.IGNORECASE)
RE_MISSENSE_HGVS = re.compile(r"p\.[A-Za-z]{3}\d+[A-Za-z]{3}")


def infer_consequence(variant_str: str) -> str:
    """Infer a coarse consequence category from an HGVS-ish variant string."""
    v = variant_str or ""
    # Order matters: check specific patterns first.
    if RE_FRAMESHIFT.search(v):
        return "frameshift"
    if RE_STOP_GAIN.search(v) or re.search(r"\*\)", v):
        return "stop_gain"
    if RE_START_LOSS.search(v):
        return "start_loss"
    if RE_SPLICE_CANONICAL.search(v):
        return "splice_canonical"
    if RE_SPLICE_DEEP.search(v):
        return "splice_deep"
    if RE_MISSENSE_HGVS.search(v):
        return "missense"
    if RE_INFRAME_INDEL.search(v):
        return "inframe"
    return "unknown"


def score_curated_variant(
    consequence: str,
    priority: str,
) -> Tuple[float, float, float, float]:
    """Return (adj_score, adj_dn, adj_lof, adj_gof) for a curated variant.

    GOF is always 0 here — the curated list is a "this is pathogenic"
    list, and GOF would require specific evidence we don't have at this
    layer.
    """
    key = consequence
    if consequence == "missense":
        pr = (priority or "").strip().lower()
        if pr.startswith("high"):
            key = "missense_high"
        elif pr.startswith("med"):
            key = "missense_medium"
        elif pr.startswith("low"):
            key = "missense_low"
        else:
            key = "missense_medium"
    weights = CONSEQUENCE_WEIGHTS.get(key, CONSEQUENCE_WEIGHTS["unknown"])
    score = weights["score"]
    dn = score * weights["dn_frac"]
    lof = score * weights["lof_frac"]
    return score, dn, lof, 0.0


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


def load_curated_variants(
    tsv_path: Path,
    exclude_mt: bool = True,
    exclude_noise_families: bool = True,
) -> List[Variant]:
    """Load a curated TSV with optional CASCADE score columns.

    Supported schemas:
      - Minimal curated: GENE / VARIANT / SYSTEMS / PRIORITY
      - Merged: same + CASCADE_final_score / CASCADE_final_classification /
        CASCADE_filtered_{dn,lof,gof}_score

    Behavior:
      - Non-missense consequences (frameshift, stop_gain, splice, etc.)
        always use the consequence-based override (CASCADE doesn't score
        these correctly anyway).
      - Missense: if CASCADE scores are present, use filtered_dn/lof/gof
        as adj_dn/adj_lof/adj_gof, and final_score as adj_score. Otherwise,
        fall back to PRIORITY-based shim.
    """
    from cumbursum.burden_scorer import _compile_noise_patterns, _is_noise_gene
    noise_patterns = _compile_noise_patterns() if exclude_noise_families else []

    tsv_path = Path(tsv_path)
    variants: List[Variant] = []
    cons_counts: Dict[str, int] = {}
    score_sources: Dict[str, int] = {"cascade": 0, "consequence_override": 0, "priority_shim": 0}
    skipped_mt = skipped_noise = 0

    with tsv_path.open("r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldmap = {k.lower().strip(): k for k in (reader.fieldnames or [])}
        gene_col = fieldmap.get("gene")
        var_col = fieldmap.get("variant")
        pri_col = fieldmap.get("priority")
        # CASCADE overlay columns (optional)
        cas_final_col = fieldmap.get("cascade_final_score")
        cas_class_col = fieldmap.get("cascade_final_classification")
        cas_dn_col = fieldmap.get("cascade_filtered_dn_score")
        cas_lof_col = fieldmap.get("cascade_filtered_lof_score")
        cas_gof_col = fieldmap.get("cascade_filtered_gof_score")
        has_cascade_cols = bool(cas_dn_col and cas_lof_col and cas_gof_col)

        if not (gene_col and var_col):
            raise ValueError(
                f"Curated TSV missing required columns. Found: {reader.fieldnames}"
            )
        for row in reader:
            gene = (row.get(gene_col) or "").strip()
            if not gene:
                continue
            if exclude_mt and gene.startswith("MT-"):
                skipped_mt += 1
                continue
            if noise_patterns and _is_noise_gene(gene, noise_patterns):
                skipped_noise += 1
                continue
            variant_str = (row.get(var_col) or "").strip()
            priority = (row.get(pri_col) or "").strip() if pri_col else ""
            consequence = infer_consequence(variant_str)
            cons_counts[consequence] = cons_counts.get(consequence, 0) + 1

            # Decide score source.
            cascade_final = (row.get(cas_final_col) or "").strip() if has_cascade_cols else ""
            use_cascade = (
                consequence == "missense"
                and has_cascade_cols
                and cascade_final  # non-empty = CASCADE actually scored this row
            )

            classification = ""
            if use_cascade:
                adj_dn = _to_float(row.get(cas_dn_col))
                adj_lof = _to_float(row.get(cas_lof_col))
                adj_gof = _to_float(row.get(cas_gof_col))
                adj_score = _to_float(row.get(cas_final_col))
                classification = (row.get(cas_class_col) or "").strip() if cas_class_col else ""
                score_sources["cascade"] += 1
            elif consequence != "missense" and consequence != "unknown":
                # Frameshift / stop_gain / splice / inframe -> deterministic override
                adj_score, adj_dn, adj_lof, adj_gof = score_curated_variant(consequence, priority)
                # Non-missense with adj_score = 1.0 is LP-equivalent by our override
                if adj_score >= 1.0:
                    classification = "LP (consequence)"
                elif adj_score >= 0.7:
                    classification = "LP-equiv (consequence)"
                else:
                    classification = "VUS (consequence)"
                score_sources["consequence_override"] += 1
            else:
                # Missense without CASCADE, or unknown consequence -> priority shim
                adj_score, adj_dn, adj_lof, adj_gof = score_curated_variant(consequence, priority)
                classification = f"priority:{priority}" if priority else "unscored"
                score_sources["priority_shim"] += 1

            variants.append(Variant(
                gene=gene,
                adj_score=adj_score,
                adj_dn=adj_dn,
                adj_lof=adj_lof,
                adj_gof=adj_gof,
                hgvs=variant_str,
                consequence=consequence,
                classification=classification,
            ))

    logger.info("Loaded %d curated variants from %s (MT skipped=%d, noise skipped=%d)",
                len(variants), tsv_path.name, skipped_mt, skipped_noise)
    logger.info("  Consequence breakdown: %s", dict(sorted(cons_counts.items())))
    logger.info("  Score sources: %s", score_sources)
    return variants


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--curated", required=True)
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()
    vs = load_curated_variants(Path(args.curated))
    print(f"\n{len(vs)} variants loaded. Top {args.top} by adj_score:")
    for v in sorted(vs, key=lambda v: v.adj_score, reverse=True)[:args.top]:
        print(f"  {v.gene:12s}  cons={v.consequence:16s}  adj_score={v.adj_score:.2f}  dn={v.adj_dn:.2f}  lof={v.adj_lof:.2f}  hgvs={v.hgvs}")
