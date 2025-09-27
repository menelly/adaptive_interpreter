# universal_router.py
"""
🌟 NOVA'S UNIVERSAL LATTICE DISRUPTION ROUTER
The main entry point for all lattice disruption analysis.

Routes variants to appropriate protein family analyzers based on
structural context and sequence features.
"""

from typing import Dict, Any, Tuple, List
from .utils import aggregate_lattice_features, format_lattice_explanation
from .family_detector import detect_protein_family
from .collagen_analyzer import analyze_collagen_lattice_disruption
from .coiled_coil_analyzer import analyze_coiled_coil_lattice_disruption
from .ion_channel_analyzer import analyze_ion_channel_lattice_disruption
from .generic_analyzer import analyze_generic_lattice_disruption


def analyze_lattice_disruption(seq: str, pos1: int, ref: str, alt: str, context: Dict = None) -> Tuple[float, List[Dict], str]:
    """
    🧬 NOVA'S UNIVERSAL LATTICE DISRUPTION ANALYZER
    
    Routes to appropriate protein family analyzer based on structural context.
    
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
    
    # Detect protein family for routing
    family = detect_protein_family(seq, context)
    
    # Route to appropriate analyzer
    if family == "COLLAGEN":
        return analyze_collagen_lattice_disruption(seq, pos1, ref, alt, context)
    elif family == "COILED_COIL":
        return analyze_coiled_coil_lattice_disruption(seq, pos1, ref, alt, context)
    elif family == "ION_CHANNEL":
        return analyze_ion_channel_lattice_disruption(seq, pos1, ref, alt, context)
    else:
        return analyze_generic_lattice_disruption(seq, pos1, ref, alt, context)


# Functions moved to utils.py to prevent circular imports
