# collagen_analyzer.py
"""
🧬 NOVA'S COLLAGEN LATTICE DISRUPTION ANALYZER
Deep analysis of collagen triple helix disruption mechanisms.

Handles Gly-X-Y triplet analysis, hydroxyproline chemistry,
crosslink site disruption, and collagen-specific geometry.
"""

from typing import Dict, List, Tuple
import math
from .utils import aggregate_lattice_features, format_lattice_explanation
from ..amino_acid_props import get_amino_acid_volume_change, get_charge_change


def analyze_collagen_lattice_disruption(seq: str, pos1: int, ref: str, alt: str, context: Dict = None) -> Tuple[float, List[Dict], str]:
    """
    🧬 NOVA'S COLLAGEN-SPECIFIC LATTICE DISRUPTION ANALYSIS
    
    Analyzes collagen triple helix disruption through multiple mechanisms:
    - Gly-X-Y triplet disruption
    - Hydroxyproline site loss
    - Crosslink site disruption
    - Triple helix geometry bombs
    
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
    
    # 1. Gly-X-Y triplet analysis
    features.extend(analyze_gly_xy_triplet(seq, pos1, ref, alt))
    
    # 2. Hydroxyproline site analysis
    features.extend(analyze_hydroxyproline_sites(seq, pos1, ref, alt))
    
    # 3. Crosslink site analysis
    features.extend(analyze_crosslink_sites(seq, pos1, ref, alt))
    
    # 4. Triple helix geometry analysis
    features.extend(analyze_triple_helix_geometry(seq, pos1, ref, alt))
    
    # 5. Regional importance scaling
    importance = get_collagen_region_importance(seq, pos1)
    for feature in features:
        feature["weight"] *= importance
    
    # Aggregate using Nova's root-sum-square method
    score = aggregate_lattice_features(features)
    explanation = format_lattice_explanation("collagen_lattice_disruption", features)
    
    return score, features, explanation


def analyze_gly_xy_triplet(seq: str, pos1: int, ref: str, alt: str) -> List[Dict]:
    """
    🎯 GLY-X-Y TRIPLET DISRUPTION ANALYSIS
    
    Collagen's Gly-X-Y repeat is mandatory for triple helix formation.
    """
    features = []
    
    # Get triplet position (0=Gly, 1=X, 2=Y)
    triplet_pos = get_triplet_position(seq, pos1)
    
    if triplet_pos == 0 and ref.upper() == "G":
        # Glycine loss at mandatory position - CATASTROPHIC
        features.append({
            "feature": "glycine_loss_mandatory",
            "weight": 0.9,
            "description": "Glycine loss in Gly-X-Y repeat"
        })
    
    elif triplet_pos == 2 and ref.upper() == "P":
        # Proline loss in Y position (often hydroxyproline)
        features.append({
            "feature": "proline_loss_Y_position", 
            "weight": 0.8,
            "description": "Proline loss in Y position of Gly-X-Y"
        })
    
    elif triplet_pos == 1:
        # X position constraints
        features.extend(analyze_x_position_constraints(ref, alt))
    
    return features


def analyze_hydroxyproline_sites(seq: str, pos1: int, ref: str, alt: str) -> List[Dict]:
    """
    🧪 HYDROXYPROLINE SITE ANALYSIS
    
    Hydroxyproline (4-hydroxyproline) is critical for collagen stability.
    """
    features = []
    
    if ref.upper() == "P" and is_likely_hydroxyproline_site(seq, pos1):
        features.append({
            "feature": "hydroxyproline_loss",
            "weight": 0.7,
            "description": "Loss of likely hydroxyproline site"
        })
    
    return features


def analyze_crosslink_sites(seq: str, pos1: int, ref: str, alt: str) -> List[Dict]:
    """
    🔗 CROSSLINK SITE DISRUPTION ANALYSIS
    
    Lysine and hydroxylysine form aldol condensation crosslinks.
    """
    features = []
    
    if ref.upper() == "K" and is_crosslink_site(seq, pos1):
        features.append({
            "feature": "crosslink_site_loss",
            "weight": 0.4,
            "description": "Loss of lysine crosslink site"
        })
    
    return features


def analyze_triple_helix_geometry(seq: str, pos1: int, ref: str, alt: str) -> List[Dict]:
    """
    📐 TRIPLE HELIX GEOMETRY DISRUPTION
    
    Certain amino acids cannot fit in the tight collagen geometry.
    """
    features = []
    
    # Bulky residues that break triple helix
    bulky_residues = ["W", "F", "Y"]
    if alt.upper() in bulky_residues and is_in_triple_helix_region(seq, pos1):
        features.append({
            "feature": "bulky_residue_helix",
            "weight": 0.9,
            "description": f"Bulky residue {alt} disrupts triple helix"
        })
    
    # Inappropriate proline introduction
    if alt.upper() == "P" and ref.upper() != "P":
        triplet_pos = get_triplet_position(seq, pos1)
        if triplet_pos != 2:  # Not Y position
            features.append({
                "feature": "inappropriate_proline",
                "weight": 0.6,
                "description": "Proline introduced outside Y position"
            })
    
    # Size/charge mismatches
    if is_in_triple_helix_region(seq, pos1):
        size_change = get_amino_acid_volume_change(ref, alt)
        if abs(size_change) > 50:  # Significant size change
            features.append({
                "feature": "size_mismatch_collagen",
                "weight": 0.4,
                "description": f"Size change disrupts packing"
            })
        
        charge_change = get_charge_change(ref, alt)
        if abs(charge_change) > 0:
            features.append({
                "feature": "charge_disruption_collagen",
                "weight": 0.3,
                "description": "Charge change in structured region"
            })
    
    return features


def get_triplet_position(seq: str, pos1: int) -> int:
    """
    🎯 GET POSITION WITHIN GLY-X-Y TRIPLET
    
    Returns:
        0 = Gly position
        1 = X position  
        2 = Y position
        -1 = Not in Gly-X-Y repeat
    """
    # Find the best frame with Gly-X-Y repeats around this position
    seq_upper = seq.upper()
    
    # Check frames around the position
    for frame in range(3):
        # Look for Gly pattern in this frame
        gly_positions = []
        for i in range(frame, len(seq_upper), 3):
            if seq_upper[i] == 'G':
                gly_positions.append(i + 1)  # Convert to 1-based
        
        # If we have enough Gly positions and our position fits the pattern
        if len(gly_positions) >= 3:
            for gly_pos in gly_positions:
                if gly_pos == pos1:
                    return 0  # Gly position
                elif gly_pos + 1 == pos1:
                    return 1  # X position
                elif gly_pos + 2 == pos1:
                    return 2  # Y position
    
    return -1  # Not in clear Gly-X-Y pattern


def is_likely_hydroxyproline_site(seq: str, pos1: int) -> bool:
    """Check if proline is likely to be hydroxylated (simplified heuristic)."""
    triplet_pos = get_triplet_position(seq, pos1)
    return triplet_pos == 2  # Y position prolines are often hydroxylated


def is_crosslink_site(seq: str, pos1: int) -> bool:
    """Check if lysine is in a crosslink-forming context (simplified)."""
    # Simplified: assume lysines in collagen regions can form crosslinks
    return is_in_triple_helix_region(seq, pos1)


def is_in_triple_helix_region(seq: str, pos1: int) -> bool:
    """Check if position is in a triple helix region."""
    triplet_pos = get_triplet_position(seq, pos1)
    return triplet_pos != -1  # In Gly-X-Y repeat


def get_collagen_region_importance(seq: str, pos1: int) -> float:
    """
    🎯 REGIONAL IMPORTANCE SCALING
    
    Different collagen regions have different importance:
    - N-terminal propeptide: lower impact
    - Central triple helix: MAXIMUM impact
    - C-terminal propeptide: moderate impact
    """
    seq_len = len(seq)
    relative_pos = pos1 / seq_len
    
    # Simplified regional scaling
    if 0.1 <= relative_pos <= 0.9:
        return 1.0  # Central region - maximum impact
    elif relative_pos < 0.1:
        return 0.7  # N-terminal - moderate impact
    else:
        return 0.8  # C-terminal - moderate impact


def analyze_x_position_constraints(ref: str, alt: str) -> List[Dict]:
    """
    Analyze constraints on X position in Gly-X-Y triplet.
    X position has fewer constraints than Gly or Y positions.
    """
    features = []
    
    # X position is more tolerant, but still has some constraints
    if alt.upper() in ["W", "F"]:  # Very bulky residues
        features.append({
            "feature": "bulky_x_position",
            "weight": 0.3,
            "description": "Bulky residue in X position"
        })
    
    return features
