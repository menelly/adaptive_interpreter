"""
🧬 DOMAIN-AWARE GENETICS ANALYZERS 🧬

Revolutionary genetics analysis system with universal domain awareness.

Key analyzers:
- lof_analyzer.py: Loss of Function analysis with domain awareness
- dn_analyzer.py: Dominant Negative analysis  
- gof_variant_analyzer.py: Gain of Function analysis
- smart_protein_analyzer.py: Protein structure analysis

All analyzers now integrate with UniversalProteinAnnotator for domain-aware scoring.
"""

from .lof_analyzer import LOFAnalyzer
from .gof_variant_analyzer import GOFVariantAnalyzer
from .smart_protein_analyzer import SmartProteinAnalyzer
from .population_frequency_analyzer import PopulationFrequencyAnalyzer

__all__ = ['LOFAnalyzer', 'GOFVariantAnalyzer', 'SmartProteinAnalyzer', 'PopulationFrequencyAnalyzer']
