# generic_analyzer.py
"""
🌍 NOVA'S GENERIC LATTICE DISRUPTION ANALYZER
Universal structural disruption analysis for proteins without
specialized family-specific analyzers.

Handles general geometry bombs, secondary structure disruption,
and universal packing constraints.
"""

from typing import Dict, List, Tuple
from .utils import aggregate_lattice_features, format_lattice_explanation
from ..amino_acid_props import get_amino_acid_volume_change, get_charge_change


def analyze_generic_lattice_disruption(seq: str, pos1: int, ref: str, alt: str, context: Dict = None) -> Tuple[float, List[Dict], str]:
    """
    🌍 NOVA'S GENERIC LATTICE DISRUPTION ANALYSIS
    
    Universal structural disruption analysis for any protein:
    - Secondary structure disruption
    - Buried region disruption
    - General geometry bombs
    - Packing constraint violations
    
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
    
    # 1. Secondary structure disruption
    features.extend(analyze_secondary_structure_disruption(seq, pos1, ref, alt))
    
    # 2. Buried region disruption
    features.extend(analyze_buried_region_disruption(seq, pos1, ref, alt, context))
    
    # 3. General geometry bombs
    features.extend(analyze_geometry_bombs(seq, pos1, ref, alt))
    
    # 4. Packing constraint violations
    features.extend(analyze_packing_constraints(seq, pos1, ref, alt))
    
    # Aggregate using Nova's method
    score = aggregate_lattice_features(features)
    explanation = format_lattice_explanation("generic_lattice_disruption", features)
    
    return score, features, explanation


def analyze_secondary_structure_disruption(seq: str, pos1: int, ref: str, alt: str) -> List[Dict]:
    """
    🌀 SECONDARY STRUCTURE DISRUPTION ANALYSIS
    
    Detect disruption of alpha helices and beta sheets.
    """
    features = []
    
    # Proline is a helix breaker
    if alt.upper() == "P" and ref.upper() != "P":
        if is_in_alpha_helix(seq, pos1):
            features.append({
                "feature": "helix_breaker",
                "weight": 0.5,
                "description": "Proline disrupts alpha helix"
            })
        elif is_in_beta_sheet(seq, pos1):
            features.append({
                "feature": "sheet_breaker", 
                "weight": 0.4,
                "description": "Proline disrupts beta sheet"
            })
    
    # Glycine loss can rigidify flexible regions
    if ref.upper() == "G" and alt.upper() != "G":
        if is_in_loop_region(seq, pos1):
            features.append({
                "feature": "flexibility_loss",
                "weight": 0.3,
                "description": "Glycine loss reduces flexibility"
            })
    
    return features


def analyze_buried_region_disruption(seq: str, pos1: int, ref: str, alt: str, context: Dict) -> List[Dict]:
    """
    🏔️ BURIED REGION DISRUPTION ANALYSIS
    
    Buried regions have strict constraints on size and charge.
    """
    features = []
    
    if is_buried_position(seq, pos1, context):
        # Size changes in buried regions
        size_change = get_amino_acid_volume_change(ref, alt)
        if abs(size_change) > 30:  # Moderate threshold for buried regions
            features.append({
                "feature": "buried_size_change",
                "weight": 0.6,
                "description": f"Size change in buried region"
            })
        
        # Charge changes in buried regions
        charge_change = get_charge_change(ref, alt)
        if abs(charge_change) > 0:
            features.append({
                "feature": "buried_charge_flip",
                "weight": 0.7,
                "description": "Charge change in buried region"
            })
        
        # Hydrophobic to polar changes
        if is_hydrophobic(ref) and is_polar(alt):
            features.append({
                "feature": "buried_polar_introduction",
                "weight": 0.6,
                "description": "Polar residue in hydrophobic core"
            })
    
    return features


def analyze_geometry_bombs(seq: str, pos1: int, ref: str, alt: str) -> List[Dict]:
    """
    💣 GEOMETRY BOMB ANALYSIS
    
    Detect substitutions that are generally disruptive.
    """
    features = []
    
    # Bulky residues in structured regions
    if alt.upper() in ["W", "F", "Y"] and is_in_structured_region(seq, pos1):
        features.append({
            "feature": "bulky_structured_disruption",
            "weight": 0.4,
            "description": f"Bulky residue {alt} in structured region"
        })
    
    # Cysteine changes (disulfide bond disruption)
    if ref.upper() == "C" or alt.upper() == "C":
        features.append({
            "feature": "cysteine_change",
            "weight": 0.5,
            "description": "Cysteine change may disrupt disulfide bonds"
        })
    
    return features


def analyze_packing_constraints(seq: str, pos1: int, ref: str, alt: str) -> List[Dict]:
    """
    📦 PACKING CONSTRAINT ANALYSIS
    
    General protein packing constraints.
    """
    features = []
    
    # Very large size changes
    size_change = get_amino_acid_volume_change(ref, alt)
    if abs(size_change) > 80:  # Very large change
        features.append({
            "feature": "extreme_size_change",
            "weight": 0.5,
            "description": f"Extreme size change may disrupt packing"
        })
    
    # Charge reversals
    ref_charge = get_residue_charge(ref)
    alt_charge = get_residue_charge(alt)
    if ref_charge * alt_charge < 0:  # Opposite charges
        features.append({
            "feature": "charge_reversal",
            "weight": 0.4,
            "description": "Charge reversal may disrupt interactions"
        })
    
    return features


# Helper functions for generic analysis

def is_in_alpha_helix(seq: str, pos1: int) -> bool:
    """
    Simplified alpha helix prediction.
    
    Real implementation would use sophisticated secondary structure prediction.
    """
    # Placeholder: assume helical regions based on amino acid propensities
    helix_favoring = set("AELKR")
    
    # Check local sequence context
    window_start = max(0, pos1 - 3)
    window_end = min(len(seq), pos1 + 3)
    window = seq[window_start:window_end].upper()
    
    helix_score = sum(1 for aa in window if aa in helix_favoring) / len(window)
    return helix_score > 0.6


def is_in_beta_sheet(seq: str, pos1: int) -> bool:
    """
    Simplified beta sheet prediction.
    """
    # Placeholder: assume sheet regions based on amino acid propensities
    sheet_favoring = set("VIFYWTC")
    
    window_start = max(0, pos1 - 3)
    window_end = min(len(seq), pos1 + 3)
    window = seq[window_start:window_end].upper()
    
    sheet_score = sum(1 for aa in window if aa in sheet_favoring) / len(window)
    return sheet_score > 0.6


def is_in_loop_region(seq: str, pos1: int) -> bool:
    """
    Simplified loop region detection.
    """
    # If not clearly helix or sheet, assume loop
    return not (is_in_alpha_helix(seq, pos1) or is_in_beta_sheet(seq, pos1))


def is_in_structured_region(seq: str, pos1: int) -> bool:
    """
    Check if position is in a structured region (helix or sheet).
    """
    return is_in_alpha_helix(seq, pos1) or is_in_beta_sheet(seq, pos1)


def is_buried_position(seq: str, pos1: int, context: Dict) -> bool:
    """
    Determine if position is buried in protein core.
    
    Simplified heuristic - real analysis would use structural data.
    """
    # Check for domain annotations suggesting buried regions
    domains = context.get("domains", [])
    for domain in domains:
        desc = domain.get("description", "").lower()
        if any(term in desc for term in ["core", "domain", "fold"]):
            start = domain.get("start", 0)
            end = domain.get("end", len(seq))
            if start <= pos1 <= end:
                return True
    
    # Fallback: assume central regions are more likely buried
    seq_len = len(seq)
    relative_pos = pos1 / seq_len
    return 0.2 <= relative_pos <= 0.8


def is_hydrophobic(aa: str) -> bool:
    """Check if amino acid is hydrophobic."""
    return aa.upper() in "AILMFWYV"


def is_polar(aa: str) -> bool:
    """Check if amino acid is polar."""
    return aa.upper() in "STNQCYH"


def get_residue_charge(aa: str) -> int:
    """Get formal charge of amino acid at physiological pH."""
    positive = "KR"
    negative = "DE"
    
    aa_upper = aa.upper()
    if aa_upper in positive:
        return 1
    elif aa_upper in negative:
        return -1
    else:
        return 0
