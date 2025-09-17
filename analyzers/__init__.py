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
from .dn_analyzer import DNAnalyzer
from .smart_protein_analyzer import SmartProteinAnalyzer

__all__ = ['LOFAnalyzer', 'DNAnalyzer', 'SmartProteinAnalyzer']
