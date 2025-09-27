# lattice_disruption/__init__.py
"""
🧬 NOVA'S UNIVERSAL LATTICE DISRUPTION FRAMEWORK
Built by Nova (OpenAI) & Ace - 2025

Universal structural sabotage analysis for dominant negative mechanisms.
Modular architecture supports protein family-specific analyzers.

Part of the Revolutionary Genomics Platform
"""

from .universal_router import analyze_lattice_disruption
from .family_detector import detect_protein_family
from .collagen_analyzer import analyze_collagen_lattice_disruption
from .coiled_coil_analyzer import analyze_coiled_coil_lattice_disruption
from .ion_channel_analyzer import analyze_ion_channel_lattice_disruption
from .generic_analyzer import analyze_generic_lattice_disruption

__all__ = [
    'analyze_lattice_disruption',
    'detect_protein_family',
    'analyze_collagen_lattice_disruption',
    'analyze_coiled_coil_lattice_disruption', 
    'analyze_ion_channel_lattice_disruption',
    'analyze_generic_lattice_disruption'
]

__version__ = "1.0.0"
__authors__ = ["Nova (OpenAI)", "Ace"]
