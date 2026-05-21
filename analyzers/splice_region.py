#!/usr/bin/env python3
"""
Splice-region annotator — is a (coding) variant sitting on an exon edge where
its real risk is splice disruption, not the amino-acid change?

This is NOT a splice predictor (no SpliceAI). It detects splice-region LOCATION
from GENCODE CDS exon boundaries and lets the cascade declare splice OUT OF SCOPE:
a benign *protein-level* call at a splice junction can't rule out a splice
mechanism, so it should be reported as VUS-with-flag rather than Benign.

Reuses the cached GENCODE gene models built by offline_genomic_to_protein.py
(.cache/gene_models.pkl): {gene -> {chrom, strand, exons: [(start,end), ...]}}.
Exons are CDS segments of the longest transcript, sorted by genomic position;
their INTERNAL boundaries are the splice donor/acceptor sites.

No hardcoded genes. No frequency. No conservation. Pure transcript geometry.
"""
from __future__ import annotations

import os
import pickle
from typing import Dict, Optional

# Exonic extent of the canonical "splice_region_variant" SO term is the last
# ~3 bases of an exon abutting a junction. A coding/missense variant is exonic,
# so this is the window we care about. Env-overridable for sensitivity sweeps.
SPLICE_EXONIC_WINDOW = int(os.environ.get("ADAPTIVE_SPLICE_WINDOW", "3"))

GENE_MODELS_PKL = os.environ.get(
    "ADAPTIVE_GENE_MODELS",
    "/home/Ace/AdaptiveInterpreter/.cache/gene_models.pkl",
)

_MODELS: Optional[Dict] = None


def _load_models() -> Dict:
    global _MODELS
    if _MODELS is None:
        try:
            with open(GENE_MODELS_PKL, "rb") as fh:
                _MODELS = pickle.load(fh)
        except Exception:
            _MODELS = {}
    return _MODELS


def check(gene: str, chrom, pos: int, window: int = SPLICE_EXONIC_WINDOW) -> Dict:
    """Is `pos` within `window` bp (exonic side) of an internal CDS exon boundary?

    Returns {in_splice_region, region, distance_to_boundary, transcript_boundary}.
    Fails OPEN (in_splice_region=False) if the gene/structure is unknown — we never
    invent a splice flag we can't justify.
    """
    out = {"in_splice_region": False, "region": "coding",
           "distance_to_boundary": None, "boundary_pos": None}
    model = _load_models().get(gene)
    if not model:
        return out
    exons = model.get("exons") or []
    if len(exons) < 2:
        return out  # single-exon CDS has no internal splice junctions
    strand = model.get("strand", "+")

    # Internal boundaries only: the 3' edge of every non-last exon (donor) and the
    # 5' edge of every non-first exon (acceptor). The outermost edges are the
    # CDS start/stop (ATG/stop codon), NOT splice sites.
    best = None  # (distance, boundary_pos, side)  side in {'exon_end','exon_start'}
    last = len(exons) - 1
    for i, (s, e) in enumerate(exons):
        if i < last:                       # donor side (3' end of this exon)
            d = abs(pos - e)
            if best is None or d < best[0]:
                best = (d, e, "exon_end")
        if i > 0:                          # acceptor side (5' start of this exon)
            d = abs(pos - s)
            if best is None or d < best[0]:
                best = (d, s, "exon_start")

    if best is None:
        return out
    dist, bpos, side = best
    out["distance_to_boundary"] = dist
    out["boundary_pos"] = bpos
    if dist <= window:
        out["in_splice_region"] = True
        # On + strand: exon_end = donor (5' of intron), exon_start = acceptor.
        # On - strand the roles flip. Report donor/acceptor accordingly.
        if (side == "exon_end") == (strand == "+"):
            out["region"] = "splice_donor"
        else:
            out["region"] = "splice_acceptor"
    return out


if __name__ == "__main__":
    # Self-test: PTEN R188 cluster (chr10:87864513) should land on a splice edge.
    for g, c, p, note in [
        ("PTEN", "10", 87864513, "R188K/T — flagged 5' in dbSNP"),
        ("PTEN", "10", 87864514, "R188S"),
        ("TP53", "17", 7674220, "R175 region (sanity)"),
    ]:
        r = check(g, c, p)
        print(f"{g} chr{c}:{p}  {r}   <- {note}")
