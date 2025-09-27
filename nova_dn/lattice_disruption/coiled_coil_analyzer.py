# coiled_coil_analyzer.py
"""
🌀 NOVA'S COILED-COIL LATTICE DISRUPTION ANALYZER
Analysis of coiled-coil heptad repeat disruption mechanisms.

Handles heptad repeat packing, hydrophobic core disruption,
salt bridge analysis, and coiled-coil specific geometry.
"""

from typing import Dict, List, Tuple
from .utils import aggregate_lattice_features, format_lattice_explanation
from ..amino_acid_props import get_charge_change


def analyze_coiled_coil_lattice_disruption(seq: str, pos1: int, ref: str, alt: str, context: Dict = None) -> Tuple[float, List[Dict], str]:
    """
    🌀 NOVA'S COILED-COIL LATTICE DISRUPTION ANALYSIS
    
    Analyzes coiled-coil disruption through heptad repeat mechanisms:
    - Hydrophobic core disruption (a/d positions)
    - Salt bridge disruption (e/g positions)
    - Bulky residue crowding
    - Heptad register shifts
    
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
    
    # 1. Heptad position analysis
    heptad_pos = get_heptad_position(seq, pos1)
    
    if heptad_pos in ["a", "d"]:
        # Core hydrophobic positions
        features.extend(analyze_hydrophobic_core(ref, alt, heptad_pos))
    elif heptad_pos in ["e", "g"]:
        # Salt bridge positions
        features.extend(analyze_salt_bridges(ref, alt, heptad_pos))
    else:
        # Other positions (b, c, f)
        features.extend(analyze_surface_positions(ref, alt, heptad_pos))
    
    # 2. General coiled-coil geometry
    features.extend(analyze_coiled_coil_geometry(seq, pos1, ref, alt))
    
    # Aggregate using Nova's method
    score = aggregate_lattice_features(features)
    explanation = format_lattice_explanation("coiled_coil_lattice_disruption", features)
    
    return score, features, explanation


def analyze_hydrophobic_core(ref: str, alt: str, position: str) -> List[Dict]:
    """
    🔥 HYDROPHOBIC CORE DISRUPTION ANALYSIS
    
    Positions 'a' and 'd' form the hydrophobic core of coiled-coils.
    """
    features = []
    
    hydrophobic = set("AILMFWYV")
    
    # Loss of hydrophobic residue at core position
    if ref.upper() in hydrophobic and alt.upper() not in hydrophobic:
        weight = 0.8 if position == "a" else 0.7  # 'a' slightly more critical
        features.append({
            "feature": f"hydrophobic_core_loss_{position}",
            "weight": weight,
            "description": f"Loss of hydrophobic residue at {position} position"
        })
    
    # Bulky residue crowding in core
    if alt.upper() in ["W", "F", "Y"]:
        features.append({
            "feature": f"bulky_core_disruption_{position}",
            "weight": 0.9,
            "description": f"Bulky residue {alt} crowds {position} position"
        })
    
    # Charged residue in hydrophobic core
    charged = set("DEKR")
    if alt.upper() in charged:
        features.append({
            "feature": f"charged_core_disruption_{position}",
            "weight": 0.8,
            "description": f"Charged residue {alt} in hydrophobic {position} position"
        })
    
    return features


def analyze_salt_bridges(ref: str, alt: str, position: str) -> List[Dict]:
    """
    ⚡ SALT BRIDGE DISRUPTION ANALYSIS
    
    Positions 'e' and 'g' often form interhelical salt bridges.
    """
    features = []
    
    charge_change = get_charge_change(ref, alt)
    
    if abs(charge_change) > 0:
        features.append({
            "feature": f"salt_bridge_disruption_{position}",
            "weight": 0.6,
            "description": f"Charge change at {position} position disrupts salt bridges"
        })
    
    return features


def analyze_surface_positions(ref: str, alt: str, position: str) -> List[Dict]:
    """
    🌊 SURFACE POSITION ANALYSIS
    
    Positions b, c, f are typically surface-exposed.
    """
    features = []
    
    # Surface positions are more tolerant, but extreme changes still matter
    if alt.upper() == "P" and ref.upper() != "P":
        features.append({
            "feature": f"proline_surface_{position}",
            "weight": 0.3,
            "description": f"Proline introduction at surface {position} position"
        })
    
    return features


def analyze_coiled_coil_geometry(seq: str, pos1: int, ref: str, alt: str) -> List[Dict]:
    """
    📐 GENERAL COILED-COIL GEOMETRY ANALYSIS
    """
    features = []
    
    # Proline is generally disruptive in coiled-coils
    if alt.upper() == "P" and ref.upper() != "P":
        features.append({
            "feature": "proline_coiled_coil_breaker",
            "weight": 0.5,
            "description": "Proline disrupts coiled-coil geometry"
        })
    
    return features


def get_heptad_position(seq: str, pos1: int) -> str:
    """
    🎯 GET HEPTAD POSITION (a, b, c, d, e, f, g)
    
    Simplified heptad assignment - in reality this requires
    sophisticated coiled-coil prediction algorithms.
    
    Returns:
        Position letter (a-g) or "unknown"
    """
    # Simplified: assume we can detect the best heptad frame
    # In practice, this would use tools like COILS or Paircoil2
    
    # Find the most likely heptad frame by looking for hydrophobic patterns
    best_frame = find_best_heptad_frame(seq)
    
    if best_frame is not None:
        # Calculate position within heptad
        relative_pos = (pos1 - 1 - best_frame) % 7
        positions = ["a", "b", "c", "d", "e", "f", "g"]
        return positions[relative_pos]
    
    return "unknown"


def find_best_heptad_frame(seq: str) -> int:
    """
    🔍 FIND BEST HEPTAD FRAME
    
    Look for the frame that maximizes hydrophobic residues at a/d positions.
    """
    seq_upper = seq.upper()
    hydrophobic = set("AILMFWYV")
    
    best_score = 0
    best_frame = None
    
    for frame in range(7):
        score = 0
        count = 0
        
        # Check a positions (frame + 0, 7, 14, ...)
        for i in range(frame, len(seq_upper), 7):
            if i < len(seq_upper):
                if seq_upper[i] in hydrophobic:
                    score += 1
                count += 1
        
        # Check d positions (frame + 3, 10, 17, ...)
        for i in range(frame + 3, len(seq_upper), 7):
            if i < len(seq_upper):
                if seq_upper[i] in hydrophobic:
                    score += 1
                count += 1
        
        if count > 0:
            normalized_score = score / count
            if normalized_score > best_score:
                best_score = normalized_score
                best_frame = frame
    
    # Only return frame if we have reasonable confidence
    return best_frame if best_score > 0.4 else None
