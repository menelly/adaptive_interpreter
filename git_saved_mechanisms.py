"""Mechanism scoring heuristics (v0):
Each returns (score: float in 0..1, features: List[dict], one_liner: str)
Inputs: seq (str), pos1 (1-based), ref (str), alt (str), context (dict)
"""
from __future__ import annotations
from typing import Dict, List, Tuple
from .amino_acid_props import delta
from . import motifs


def _w(ctx: Dict | None, key: str, default: float) -> float:
    try:
        return float(((ctx or {}).get("_weights", {})).get(key, default))
    except Exception:
        return float(default)


def _top_features(feats: List[Dict], k: int = 2) -> List[str]:
    # sort by |weight| desc
    fs = sorted(feats, key=lambda f: abs(f.get("weight", 0)), reverse=True)
    return [f["feature"] for f in fs[:k]]


def score_interface_poisoning(seq: str, pos1: int, ref: str, alt: str, context: Dict | None) -> Tuple[float, List[Dict], str]:
    """
    Score likelihood of interface poisoning mechanism.

    Interface poisoning occurs when amino acid changes disrupt protein-protein interactions
    by altering charge distribution, hydrophobicity, or introducing structural constraints
    at binding interfaces.

    Args:
        seq: Full protein sequence
        pos1: 1-based position of variant
        ref: Reference amino acid
        alt: Alternate amino acid
        context: Optional protein annotation context

    Returns:
        Tuple of (score, feature_list, explanation_string)
    """
    ctx = context or {}
    d = delta(ref, alt)
    feats: List[Dict] = []
    # Interface likelihood (hint from context if provided)
    interface_likelihood = float(ctx.get("interface_likelihood", 0.0))
    if interface_likelihood:
        feats.append({"feature": "interface_likelihood", "value": interface_likelihood, "weight": _w(ctx, "interface_poisoning.interface_likelihood", 0.4)})
    # Charge/hydropathy flips
    feats.append({"feature": "|d_charge|", "value": d["abs_chg"], "weight": _w(ctx, "interface_poisoning.|d_charge|", 0.25)})
    feats.append({"feature": "|d_hydropathy|", "value": d["norm_abs_hyd"], "weight": _w(ctx, "interface_poisoning.|d_hydropathy|", 0.2)})
    # Proline breaker
    if d["proline_introduced"]:
        feats.append({"feature": "proline_introduced", "value": 1, "weight": _w(ctx, "interface_poisoning.proline_introduced", 0.25)})
    # Cys at interface (disulfide/salt-bridge change)
    if d["cysteine_gain"] or d["cysteine_loss"]:
        feats.append({"feature": "cys_gain_or_loss", "value": 1, "weight": _w(ctx, "interface_poisoning.cys_gain_or_loss", 0.2)})
    # Flexible loop dampener for interface claims
    if ctx.get("flexible_loop"):
        feats.append({"feature": "flexible_loop", "value": 1, "weight": _w(ctx, "interface_poisoning.flexible_loop", -0.1)})

    score = 0.0
    for f in feats:
        score += f["weight"] * (f["value"] if isinstance(f["value"], (int, float)) else 1.0)
    score = max(0.0, min(1.0, score))
    one = " + ".join(_top_features(feats)) or "mild interface risk"
    return score, feats, f"interface poisoning due to {one}"


def score_active_site_jamming(seq: str, pos1: int, ref: str, alt: str, context: Dict | None) -> Tuple[float, List[Dict], str]:
    """
    Score likelihood of active site jamming mechanism.

    Active site jamming occurs when variants disrupt catalytic sites, substrate binding
    pockets, or cofactor binding regions, preventing normal enzymatic or binding function.

    Args:
        seq: Full protein sequence
        pos1: 1-based position of variant
        ref: Reference amino acid
        alt: Alternate amino acid
        context: Optional protein annotation context

    Returns:
        Tuple of (score, feature_list, explanation_string)
    """
    ctx = context or {}
    d = delta(ref, alt)
    feats: List[Dict] = []
    # Motif/active site hint
    if ctx.get("active_site_proximity", 0):
        feats.append({"feature": "active_site_proximity", "value": float(ctx["active_site_proximity"]), "weight": _w(ctx, "active_site_jamming.active_site_proximity", 0.4)})
    if motifs.catalytic_motif_near(seq, pos1):
        feats.append({"feature": "catalytic_motif_near", "value": 1, "weight": _w(ctx, "active_site_jamming.catalytic_motif_near", 0.3)})
    # Volume/aromatic/proline in core (approximate)
    feats.append({"feature": "|d_volume|", "value": d["norm_abs_vol"], "weight": _w(ctx, "active_site_jamming.|d_volume|", 0.2)})
    if d["aromatic_gain"] or d["aromatic_loss"]:
        feats.append({"feature": "aromatic_swap", "value": 1, "weight": _w(ctx, "active_site_jamming.aromatic_swap", 0.15)})
    if d["proline_introduced"]:
        feats.append({"feature": "proline_introduced", "value": 1, "weight": _w(ctx, "active_site_jamming.proline_introduced", 0.25)})
    # Flexible loops tend to be solvent-exposed; slightly damp active-site jamming confidence
    if ctx.get("flexible_loop"):
        feats.append({"feature": "flexible_loop", "value": 1, "weight": _w(ctx, "active_site_jamming.flexible_loop", -0.1)})
    # Disulfide-rich secretory context: prefer trafficking over catalytic jamming
    if ctx.get("in_disulfide_pair"):
        feats.append({"feature": "in_disulfide_pair", "value": 1, "weight": _w(ctx, "active_site_jamming.in_disulfide_pair", -0.25)})
    if ctx.get("secretory_disulfide_rich"):
        feats.append({"feature": "secretory_disulfide_rich", "value": 1, "weight": _w(ctx, "active_site_jamming.secretory_disulfide_rich", -0.05)})

    score = 0.0
    for f in feats:
        score += f["weight"] * (f["value"] if isinstance(f["value"], (int, float)) else 1.0)
    score = max(0.0, min(1.0, score))
    one = " + ".join(_top_features(feats)) or "mild active-site risk"
    return score, feats, f"active-site jamming via {one}"


def score_structural_lattice_disruption(seq: str, pos1: int, ref: str, alt: str, context: Dict | None) -> Tuple[float, List[Dict], str]:
    """
    Score likelihood of structural lattice disruption mechanism.

    Structural lattice disruption occurs when variants break critical structural motifs
    such as collagen Gly-X-Y repeats, beta-sheet registers, or alpha-helix geometry.

    Args:
        seq: Full protein sequence
        pos1: 1-based position of variant
        ref: Reference amino acid
        alt: Alternate amino acid
        context: Optional protein annotation context

    Returns:
        Tuple of (score, feature_list, explanation_string)
    """
    d = delta(ref, alt)
    feats: List[Dict] = []
    # Collagen Gly-X-Y killer
    if ref.upper() == "G" and motifs.is_collagen_gly_site(seq, pos1):
        feats.append({"feature": "collagen_Gly_site", "value": 1, "weight": _w(context, "lattice_disruption.collagen_Gly_site", 0.6)})
    # Annotation says this specific Gly is critical
    if (context or {}).get("critical_collagen_gly"):
        feats.append({"feature": "critical_collagen_gly", "value": 1, "weight": _w(context, "lattice_disruption.critical_collagen_gly", 0.1)})
    # Coiled-coil mismatch rough flag
    if motifs.rough_coiled_coil_flag(seq, pos1):
        feats.append({"feature": "coiled_coil_flag", "value": 1, "weight": _w(context, "lattice_disruption.coiled_coil_flag", 0.25)})
    # Secondary structure mismatch proxies
    if d["proline_introduced"]:
        feats.append({"feature": "proline_in_helix/strand", "value": 1, "weight": _w(context, "lattice_disruption.proline_in_helix/strand", 0.25)})

    score = 0.0
    for f in feats:
        score += f["weight"] * (f["value"] if isinstance(f["value"], (int, float)) else 1.0)
    score = max(0.0, min(1.0, score))
    one = " + ".join(_top_features(feats)) or "mild lattice risk"
    return score, feats, f"lattice disruption via {one}"


def score_trafficking_maturation(seq: str, pos1: int, ref: str, alt: str, context: Dict | None) -> Tuple[float, List[Dict], str]:
    """
    Score likelihood of trafficking/maturation disruption mechanism.

    Trafficking/maturation disruption occurs when variants interfere with protein folding,
    disulfide bond formation, glycosylation, or other post-translational modifications
    required for proper protein maturation and cellular localization.

    Args:
        seq: Full protein sequence
        pos1: 1-based position of variant
        ref: Reference amino acid
        alt: Alternate amino acid
        context: Optional protein annotation context

    Returns:
        Tuple of (score, feature_list, explanation_string)
    """
    d = delta(ref, alt)
    feats: List[Dict] = []
    # Disulfide perturbation
    if d["cysteine_gain"] or d["cysteine_loss"]:
        feats.append({"feature": "disulfide_network_change", "value": 1, "weight": _w(context, "trafficking_maturation.disulfide_network_change", 0.5)})
    # N-glycan gain/loss
    g_gain, g_lost = motifs.nglyc_gain_loss(seq, pos1, ref, alt)
    if g_gain:
        feats.append({"feature": "NXS/T_gained", "value": 1, "weight": _w(context, "trafficking_maturation.NXS/T_gained", 0.3)})
    if g_lost:
        feats.append({"feature": "NXS/T_lost", "value": 1, "weight": _w(context, "trafficking_maturation.NXS/T_lost", 0.25)})
    # Disulfide-region sensitivity
    if (context or {}).get("in_disulfide_pair"):
        feats.append({"feature": "in_disulfide_pair", "value": 1, "weight": _w(context, "trafficking_maturation.in_disulfide_pair", 0.2)})
        aromatic = set("FWY")
        polar = set("STNQDEKRH")
        if ref.upper() in aromatic and alt.upper() in polar:
            feats.append({"feature": "aromatic_to_polar_in_disulfide_region", "value": 1, "weight": _w(context, "trafficking_maturation.aromatic_to_polar_in_disulfide_region", 0.25)})

    score = 0.0
    for f in feats:
        score += f["weight"] * (f["value"] if isinstance(f["value"], (int, float)) else 1.0)
    score = max(0.0, min(1.0, score))
    one = " + ".join(_top_features(feats)) or "mild trafficking risk"
    return score, feats, f"trafficking/maturation mischief via {one}"
