# ion_channel_analyzer.py
"""
⚡ NOVA'S ION CHANNEL LATTICE DISRUPTION ANALYZER
Analysis of ion channel structural disruption mechanisms.

Handles pore geometry, selectivity filters, gating mechanisms,
and transmembrane domain disruption.
"""

from typing import Dict, List, Tuple
from .utils import aggregate_lattice_features, format_lattice_explanation
from ..amino_acid_props import get_amino_acid_volume_change, get_charge_change


def analyze_ion_channel_lattice_disruption(seq: str, pos1: int, ref: str, alt: str, context: Dict = None) -> Tuple[float, List[Dict], str]:
    """
    ⚡ NOVA'S ION CHANNEL LATTICE DISRUPTION ANALYSIS
    
    Analyzes ion channel disruption through multiple mechanisms:
    - Pore geometry disruption
    - Selectivity filter disruption
    - Gating mechanism disruption
    - Transmembrane domain disruption
    
    Args:
        seq: Full protein sequence
        pos1: 1-based position of variant
        ref: Reference amino acid
        alt: Alternate amino acid
        context: Optional protein annotation context
        
    Returns:
        Tuple of (score, feature_list, explanation_string)
    """
    if context is None:
        context = {}
    
    features = []
    
    # 1. Pore-lining region analysis
    if is_in_pore_lining(seq, pos1, context):
        features.extend(analyze_pore_disruption(ref, alt))
    
    # 2. Selectivity filter analysis
    if is_in_selectivity_filter(seq, pos1, context):
        features.extend(analyze_selectivity_filter(ref, alt))
    
    # 3. Gating mechanism analysis
    features.extend(analyze_gating_disruption(seq, pos1, ref, alt))
    
    # 4. Transmembrane domain analysis
    if is_in_transmembrane(seq, pos1, context):
        features.extend(analyze_transmembrane_disruption(ref, alt))
    
    # 5. General channel geometry
    features.extend(analyze_channel_geometry(seq, pos1, ref, alt))
    
    # Aggregate using Nova's method
    score = aggregate_lattice_features(features)
    explanation = format_lattice_explanation("ion_channel_lattice_disruption", features)
    
    return score, features, explanation


def analyze_pore_disruption(ref: str, alt: str) -> List[Dict]:
    """
    🕳️ PORE GEOMETRY DISRUPTION ANALYSIS
    
    Changes in pore-lining residues can alter channel conductance.
    """
    features = []
    
    # Size changes in pore
    size_change = get_amino_acid_volume_change(ref, alt)
    if abs(size_change) > 50:  # Significant size change
        weight = 0.8 if size_change > 0 else 0.7  # Expansion vs contraction
        features.append({
            "feature": "pore_radius_disruption",
            "weight": weight,
            "description": f"Size change alters pore geometry"
        })
    
    # Charge changes in pore
    charge_change = get_charge_change(ref, alt)
    if abs(charge_change) > 0:
        features.append({
            "feature": "pore_charge_disruption",
            "weight": 0.7,
            "description": "Charge change affects ion selectivity"
        })
    
    return features


def analyze_selectivity_filter(ref: str, alt: str) -> List[Dict]:
    """
    🎯 SELECTIVITY FILTER DISRUPTION ANALYSIS
    
    Selectivity filters are highly conserved and critical for function.
    """
    features = []
    
    # Any change in selectivity filter is potentially catastrophic
    features.append({
        "feature": "selectivity_filter_disruption",
        "weight": 0.9,
        "description": "Change in highly conserved selectivity filter"
    })
    
    # Specific patterns for different channel types
    if ref.upper() == "G" and is_potassium_channel_signature(ref, alt):
        features.append({
            "feature": "k_channel_gly_disruption",
            "weight": 0.95,
            "description": "Glycine disruption in K+ channel signature"
        })
    
    return features


def analyze_gating_disruption(seq: str, pos1: int, ref: str, alt: str) -> List[Dict]:
    """
    🚪 GATING MECHANISM DISRUPTION ANALYSIS
    
    Gating hinges and conformational switches are critical.
    """
    features = []
    
    # Gating hinge glycines
    if is_gating_glycine(seq, pos1) and ref.upper() == "G" and alt.upper() != "G":
        features.append({
            "feature": "gating_hinge_loss",
            "weight": 0.7,
            "description": "Loss of flexible glycine at gating hinge"
        })
    
    # Proline introduction can disrupt conformational changes
    if alt.upper() == "P" and ref.upper() != "P":
        features.append({
            "feature": "gating_rigidity",
            "weight": 0.5,
            "description": "Proline introduction reduces conformational flexibility"
        })
    
    return features


def analyze_transmembrane_disruption(ref: str, alt: str) -> List[Dict]:
    """
    🧱 TRANSMEMBRANE DOMAIN DISRUPTION ANALYSIS
    
    Transmembrane regions have specific constraints.
    """
    features = []
    
    # Charged residues in transmembrane domain
    charge_change = get_charge_change(ref, alt)
    if abs(charge_change) > 0:
        features.append({
            "feature": "tm_charge_insertion",
            "weight": 0.6,
            "description": "Charged residue in transmembrane domain"
        })
    
    # Proline in transmembrane helices
    if alt.upper() == "P" and ref.upper() != "P":
        features.append({
            "feature": "tm_helix_breaker",
            "weight": 0.7,
            "description": "Proline disrupts transmembrane helix"
        })
    
    return features


def analyze_channel_geometry(seq: str, pos1: int, ref: str, alt: str) -> List[Dict]:
    """
    📐 GENERAL CHANNEL GEOMETRY ANALYSIS
    """
    features = []
    
    # Bulky residues can disrupt channel assembly
    if alt.upper() in ["W", "F", "Y"]:
        features.append({
            "feature": "bulky_channel_disruption",
            "weight": 0.4,
            "description": f"Bulky residue {alt} may disrupt channel geometry"
        })
    
    return features


# Helper functions for channel-specific analysis

def is_in_pore_lining(seq: str, pos1: int, context: Dict) -> bool:
    """
    Determine if position is in pore-lining region.
    
    This is simplified - real analysis would use structural data
    or sophisticated prediction algorithms.
    """
    # Check context annotations
    domains = context.get("domains", [])
    for domain in domains:
        if "pore" in domain.get("description", "").lower():
            start = domain.get("start", 0)
            end = domain.get("end", len(seq))
            if start <= pos1 <= end:
                return True
    
    return False


def is_in_selectivity_filter(seq: str, pos1: int, context: Dict) -> bool:
    """
    Determine if position is in selectivity filter.
    
    Simplified detection based on known motifs.
    """
    # Look for common selectivity filter motifs
    window_start = max(0, pos1 - 5)
    window_end = min(len(seq), pos1 + 5)
    window = seq[window_start:window_end].upper()
    
    # K+ channel signature (GYG, GFG, etc.)
    if "GYG" in window or "GFG" in window:
        return True
    
    return False


def is_gating_glycine(seq: str, pos1: int) -> bool:
    """
    Determine if glycine is at a gating hinge position.
    
    Simplified heuristic - real analysis would use structural data.
    """
    # Look for glycines in potential hinge regions
    # This is a placeholder - real implementation would be more sophisticated
    return seq[pos1-1].upper() == "G"


def is_in_transmembrane(seq: str, pos1: int, context: Dict) -> bool:
    """
    Determine if position is in transmembrane domain.
    """
    # Check context annotations
    domains = context.get("domains", [])
    for domain in domains:
        if "transmembrane" in domain.get("description", "").lower():
            start = domain.get("start", 0)
            end = domain.get("end", len(seq))
            if start <= pos1 <= end:
                return True
    
    return False


def is_potassium_channel_signature(ref: str, alt: str) -> bool:
    """
    Check if this looks like a K+ channel selectivity filter disruption.
    """
    # K+ channels often have GYG or similar motifs
    return ref.upper() == "G"
