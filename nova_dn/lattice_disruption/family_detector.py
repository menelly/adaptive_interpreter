# family_detector.py
"""
🔍 PROTEIN FAMILY DETECTION FOR LATTICE DISRUPTION ROUTING
Determines which specialized analyzer to use based on sequence and context.
"""

from typing import Dict, Optional
import re


def detect_protein_family(seq: str, context: Dict = None) -> str:
    """
    🎯 DETECT PROTEIN FAMILY FOR LATTICE ANALYSIS ROUTING
    
    Args:
        seq: Protein sequence
        context: Optional annotation context
        
    Returns:
        Family string: "COLLAGEN", "COILED_COIL", "ION_CHANNEL", or "GENERIC"
    """
    if context is None:
        context = {}
    
    # Check for collagen first (most specific)
    if is_collagen_family(seq, context):
        return "COLLAGEN"
    
    # Check for coiled-coil proteins
    if is_coiled_coil_family(seq, context):
        return "COILED_COIL"
    
    # Check for ion channels
    if is_ion_channel_family(seq, context):
        return "ION_CHANNEL"
    
    # Default to generic analyzer
    return "GENERIC"


def is_collagen_family(seq: str, context: Dict) -> bool:
    """
    Detect collagen family proteins.
    
    Criteria:
    - Gly-X-Y repeat patterns
    - Function mentions collagen
    - Gene family annotations
    """
    # Check for Gly-X-Y repeats (collagen signature)
    if has_gly_xy_repeats(seq):
        return True
    
    # Check function/annotation context
    function_text = context.get("function", "").lower()
    if any(term in function_text for term in ["collagen", "triple helix", "extracellular matrix"]):
        return True
    
    # Check gene family
    gene_family = context.get("gene_family", "").upper()
    if "COLLAGEN" in gene_family:
        return True
    
    return False


def is_coiled_coil_family(seq: str, context: Dict) -> bool:
    """
    Detect coiled-coil proteins.
    
    Criteria:
    - Heptad repeat patterns
    - Function mentions coiled-coil
    - Known coiled-coil domains
    """
    # Check for heptad repeats (simplified detection)
    if has_heptad_repeats(seq):
        return True
    
    # Check function context
    function_text = context.get("function", "").lower()
    if any(term in function_text for term in ["coiled coil", "coiled-coil", "leucine zipper"]):
        return True
    
    return False


def is_ion_channel_family(seq: str, context: Dict) -> bool:
    """
    Detect ion channel proteins.
    
    Criteria:
    - Function mentions channel/transport
    - Transmembrane domains
    - Known channel families
    """
    function_text = context.get("function", "").lower()
    
    # Ion channel keywords
    channel_terms = [
        "ion channel", "sodium channel", "potassium channel", "calcium channel",
        "chloride channel", "voltage-gated", "ligand-gated", "transporter"
    ]
    
    if any(term in function_text for term in channel_terms):
        return True
    
    return False


def has_gly_xy_repeats(seq: str, min_repeats: int = 6) -> bool:
    """
    🧬 DETECT GLY-X-Y COLLAGEN REPEATS
    
    Look for Gly every 3rd position in stretches of at least min_repeats.
    """
    seq_upper = seq.upper()
    
    # Scan each frame (0, 1, 2)
    for frame in range(3):
        gly_count = 0
        max_gly_run = 0
        
        for i in range(frame, len(seq_upper), 3):
            if seq_upper[i] == 'G':
                gly_count += 1
                max_gly_run = max(max_gly_run, gly_count)
            else:
                gly_count = 0
        
        if max_gly_run >= min_repeats:
            return True
    
    return False


def has_heptad_repeats(seq: str, min_repeats: int = 3) -> bool:
    """
    🌀 DETECT HEPTAD REPEATS (COILED-COIL SIGNATURE)
    
    Look for hydrophobic residues at positions a and d in heptad repeats.
    Simplified detection - real coiled-coil prediction is more complex.
    """
    seq_upper = seq.upper()
    hydrophobic = set("AILMFWYV")
    
    # Scan for heptad patterns
    for frame in range(7):
        heptad_count = 0
        
        for i in range(frame, len(seq_upper) - 7, 7):
            # Check positions a (i) and d (i+3) for hydrophobic residues
            if seq_upper[i] in hydrophobic and seq_upper[i + 3] in hydrophobic:
                heptad_count += 1
            else:
                if heptad_count >= min_repeats:
                    return True
                heptad_count = 0
        
        if heptad_count >= min_repeats:
            return True
    
    return False
