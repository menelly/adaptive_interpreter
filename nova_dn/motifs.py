"""Lightweight motif heuristics: sequence-only, regex/logic based.
These are intentionally simple and dependency-free.
"""
from __future__ import annotations
import re
from typing import List, Tuple

# Precompile simple patterns
CATALYTIC_PATTERNS = [
    re.compile(r"H..H"),          # metal binding HxH
    re.compile(r"H[ST]D"),        # HSD/H-TD catalytic motif sketch
    re.compile(r"P\w{3}G\w{2}KT"),  # P-loopish sketch (very loose)
]


def is_collagen_gly_site(seq: str, pos1: int) -> bool:
    """True if 1-based position is a Gly in a Gly-X-Y repeating context.
    We just check local G at pos and that (pos+3) or (pos-3) is also G.
    """
    i = pos1 - 1
    if i < 0 or i >= len(seq):
        return False
    if seq[i].upper() != "G":
        return False
    # Look for periodicity
    if i + 3 < len(seq) and seq[i + 3].upper() == "G":
        return True
    if i - 3 >= 0 and seq[i - 3].upper() == "G":
        return True
    return False


def nglyc_gain_loss(seq: str, pos1: int, ref: str, alt: str) -> Tuple[bool, bool]:
    """Check N-linked glycosylation motif N-X-[ST] around the mutation.
    Returns (gained, lost).
    We consider windows starting at i-2..i with X != P.
    """
    s = list(seq)
    i = pos1 - 1
    if not (0 <= i < len(s)):
        return (False, False)
    # original context
    def has_motif(chars: List[str]) -> bool:
        for j in range(0, len(chars) - 2):
            a, b, c = chars[j : j + 3]
            if a == "N" and b != "P" and c in ("S", "T"):
                return True
        return False

    orig = s.copy()
    orig[i] = ref.upper()
    new = s.copy()
    new[i] = alt.upper()

    # check in a small window around i-2..i
    left = max(0, i - 2)
    right = min(len(s), i + 3)
    return (has_motif(new[left:right]) and not has_motif(orig[left:right]),
            has_motif(orig[left:right]) and not has_motif(new[left:right]))


def catalytic_motif_near(seq: str, pos1: int, window: int = 3) -> bool:
    i0 = max(0, pos1 - 1 - window)
    i1 = min(len(seq), pos1 - 1 + window + 1)
    region = seq[i0:i1]
    return any(p.search(region) for p in CATALYTIC_PATTERNS)


def rough_coiled_coil_flag(seq: str, pos1: int) -> bool:
    """Very rough placeholder for coiled-coil propensity: check for heptad-ish hydrophobic cadence around site.
    We detect if positions i, i+3, i+4 within window are mostly hydrophobic letters.
    """
    hydros = set("AILMVFWY")
    i = pos1 - 1
    window = seq[max(0, i - 7) : min(len(seq), i + 8)]
    if len(window) < 8:
        return False
    # sample positions modulo 7
    hits = 0
    for k in (0, 3, 4):
        j = i - len(window)//2 + k
        if 0 <= j < len(seq) and seq[j] in hydros:
            hits += 1
    return hits >= 2

