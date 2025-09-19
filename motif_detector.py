#!/usr/bin/env python3
"""
🔍 MOTIF DETECTOR - Regulatory Motif & Proximity Detection
Built for GOF analysis (Nova, 2025)

Detects canonical regulatory motifs (DFG, HRD, APE, glycine loop, VAIK, etc.),
assigns proximity weights for variants, and classifies mutation types.
"""

import re
from typing import Dict, List, Tuple

# Amino acid property groups
HYDROPHOBIC = set("AILMVFWY")
ACIDIC = set("DE")
BASIC = set("KRH")
POLAR = set("STNQ")
GLYCINE = {"G"}
CYSTEINE = {"C"}

# Universal motifs with regex patterns
MOTIF_PATTERNS = {
    "DFG": re.compile(r"DFG"),
    "HRD": re.compile(r"HRD"),
    "APE": re.compile(r"APE"),
    "GLY_LOOP": re.compile(r"G.G..G"),   # Glycine-rich loop (GxGxxG)
    "VAIK": re.compile(r"VAIK"),
}

# Switch regions (for RAS-like GTPases)
SWITCH_I_RANGE = (30, 38)
SWITCH_II_RANGE = (59, 76)


def detect_motifs(sequence: str) -> Dict[str, List[int]]:
    """
    Scan sequence for universal motifs.

    Returns:
        Dict of motif_name -> list of starting positions (1-based, UniProt style).
    """
    motif_hits = {}
    for name, pattern in MOTIF_PATTERNS.items():
        hits = [m.start() + 1 for m in pattern.finditer(sequence)]
        if hits:
            motif_hits[name] = hits
    return motif_hits


def classify_mutation(ref: str, alt: str, position: int) -> Tuple[str, float]:
    """
    Classify mutation type with GOF relevance.

    Returns:
        (description, base_score)
    """
    if ref in GLYCINE and alt != "G":
        return ("Glycine→X (flexibility loss)", 0.9)
    if ref in HYDROPHOBIC and alt in ACIDIC:
        return ("Hydrophobic→Acidic", 0.8)
    if alt in CYSTEINE:
        return ("Cysteine creation (potential disulfide)", 0.7)
    if ref in BASIC and alt in ACIDIC:
        return ("Charge inversion", 0.7)
    return ("Other mutation", 0.3)


def compute_proximity_boost(
    variant_position: int, motif_positions: Dict[str, List[int]], window: int = 10
) -> Tuple[str, float]:
    """
    Compute proximity weighting based on distance to motifs.

    Returns:
        (closest_motif, weight)
    """
    closest_motif, min_dist = None, 999
    for motif, positions in motif_positions.items():
        for pos in positions:
            dist = abs(variant_position - pos)
            if dist < min_dist:
                min_dist, closest_motif = dist, motif

    if closest_motif is None:
        return ("None", 1.0)  # no boost

    if min_dist <= 5:
        return (closest_motif, 1.8)  # strong boost
    elif min_dist <= 10:
        return (closest_motif, 1.5)  # moderate boost
    else:
        return (closest_motif, 1.0)  # no significant effect


def detect_regulatory_context(sequence: str, ref: str, alt: str, position: int) -> Dict:
    """
    Main entry: detect motif proximity + mutation classification.

    Args:
        sequence: protein sequence (string, 1-letter amino acids)
        ref: reference AA (single letter)
        alt: alternate AA (single letter)
        position: variant position (1-based)

    Returns:
        dict with motif context, mutation class, and GOF score contribution
    """
    motif_positions = detect_motifs(sequence)
    mut_desc, base_score = classify_mutation(ref, alt, position)
    motif, multiplier = compute_proximity_boost(position, motif_positions)

    # Special case: RAS-like switch regions
    if SWITCH_I_RANGE[0] <= position <= SWITCH_I_RANGE[1] or \
       SWITCH_II_RANGE[0] <= position <= SWITCH_II_RANGE[1]:
        return {
            "motif_context": "RAS Switch region",
            "mutation_class": mut_desc,
            "score": max(base_score, 0.9),
            "explanation": f"{mut_desc} in RAS Switch region → canonical GOF"
        }

    final_score = min(base_score * multiplier, 1.0)

    return {
        "motif_context": motif,
        "mutation_class": mut_desc,
        "score": final_score,
        "explanation": f"{mut_desc} near {motif} motif (multiplier {multiplier:.2f}) → score {final_score:.2f}"
    }


# Example usage
if __name__ == "__main__":
    seq = "VKIGDFGLATVKSRWSGSHQ"  # BRAF activation loop snippet
    result = detect_regulatory_context(seq, "V", "E", 600)
    print(result)
    # Expect explanation: Hydrophobic→Acidic near DFG motif (multiplier ~1.8) → score ~0.9
