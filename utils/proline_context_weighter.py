#!/usr/bin/env python3
"""
🧬 NOVA'S PROLINE CONTEXT WEIGHTER
Context-aware proline scoring system that stops the "Proline Panic Machine™"

This system provides biologically intelligent proline substitution scoring:
- HIGH IMPACT: Proline near cleavage sites, binding sites, structural bends
- LOW IMPACT: Proline in bulk triple helix regions
- CONTEXT-DEPENDENT: Based on functional importance, not blanket penalties

Created by Nova with love! 💜
"""

from typing import List, Dict, Any, Optional
import math

# 🎯 NOVA'S TUNABLE PARAMETERS - Adjust empirically
BASELINE_PROLINE_MULT = 1.10  # Small baseline nudge for proline substitutions
NEAR_WINDOW = 5               # Residues around critical features
INTERMEDIATE_WINDOW = 10      # Residues around interfaces
TRIPLE_HELIX_DOWN = 0.70     # Downweight for bulk triple helix
CLEAVAGE_NEAR_MULT = 1.50    # Upweight near cleavage sites
BINDING_NEAR_MULT = 1.50     # Upweight near binding sites
INTERFACE_NEAR_MULT = 1.40   # Upweight near interfaces
BURIAL_SCALE = 0.5           # How much burial affects penalty
CONSERVATION_SCALE = 0.6     # How much conservation affects penalty
MIN_MULT = 0.5               # Minimum multiplier cap
MAX_MULT = 2.0               # Maximum multiplier cap

class ProlineContextWeighter:
    """🚀 Nova's context-aware proline importance scoring system"""
    
    def __init__(self):
        """Initialize the Proline Context Weighter"""
        print("🧬 Initializing Nova's Proline Context Weighter - No more Proline Panic Machine™!")
    
    def _within_window(self, pos: int, feat_start: int, feat_end: int, window: int) -> bool:
        """Check if position is within window of feature"""
        return (pos >= feat_start - window) and (pos <= feat_end + window)
    
    def get_proline_multiplier(
        self,
        residue_pos: int,
        reference_aa: str,
        alt_aa: str,
        uniprot_features: List[Dict[str, Any]],
        conservation: Optional[float] = None,
        rsa: Optional[float] = None,
        secondary_structure: Optional[str] = None
    ) -> float:
        """
        🎯 Compute context-aware proline importance multiplier
        
        Args:
            residue_pos: Protein residue position (1-based)
            reference_aa: Reference amino acid
            alt_aa: Alternative amino acid
            uniprot_features: UniProt features with start, end, category
            conservation: Conservation score 0-1 (optional)
            rsa: Relative solvent accessibility 0-1 (optional)
            secondary_structure: Secondary structure prediction (optional)
            
        Returns:
            float: Multiplicative factor for variant scoring
        """
        
        # Only apply to proline substitutions
        if reference_aa.upper() != 'P' and alt_aa.upper() != 'P':
            return 1.0
            
        print(f"🧬 Nova's proline context analysis for position {residue_pos} ({reference_aa}→{alt_aa})")
        
        # Start with baseline proline multiplier
        mult = BASELINE_PROLINE_MULT
        context_factors = []
        
        # 🎯 1) Check for tolerant contexts first (downweight)
        for feature in uniprot_features:
            category = (feature.get("category") or feature.get("type") or "").lower()
            start_pos = int(feature.get("start") or feature.get("begin", {}).get("position", 0))
            end_pos = int(feature.get("end") or feature.get("end", {}).get("position", start_pos))
            
            if category in {"triple_helix", "triple-helix", "collagen_triple_helix", "region"}:
                description = (feature.get("description") or "").lower()
                if "triple" in description and "helix" in description:
                    if start_pos <= residue_pos <= end_pos:
                        mult *= TRIPLE_HELIX_DOWN
                        context_factors.append(f"triple_helix_bulk (×{TRIPLE_HELIX_DOWN})")
                        break
        
        # 🎯 2) Check proximity to high-impact features
        for feature in uniprot_features:
            category = (feature.get("category") or feature.get("type") or "").lower()
            start_pos = int(feature.get("start") or feature.get("begin", {}).get("position", 0))
            end_pos = int(feature.get("end") or feature.get("end", {}).get("position", start_pos))
            description = (feature.get("description") or "").lower()
            
            # Cleavage sites and processing domains
            if (category in {"cleavage_site", "processing", "processing_domain", "propep"} or 
                "cleavage" in description):
                if self._within_window(residue_pos, start_pos, end_pos, NEAR_WINDOW):
                    mult *= CLEAVAGE_NEAR_MULT
                    context_factors.append(f"near_cleavage (×{CLEAVAGE_NEAR_MULT})")
            
            # Active sites
            elif category in {"active_site", "act_site", "site"}:
                if self._within_window(residue_pos, start_pos, end_pos, NEAR_WINDOW):
                    mult *= CLEAVAGE_NEAR_MULT
                    context_factors.append(f"near_active_site (×{CLEAVAGE_NEAR_MULT})")
            
            # Binding sites
            elif category in {"binding_site", "binding"}:
                if self._within_window(residue_pos, start_pos, end_pos, NEAR_WINDOW):
                    mult *= BINDING_NEAR_MULT
                    context_factors.append(f"near_binding_site (×{BINDING_NEAR_MULT})")
            
            # Interfaces
            elif category in {"interface", "protein_binding_interface", "domain_interface"}:
                if self._within_window(residue_pos, start_pos, end_pos, INTERMEDIATE_WINDOW):
                    mult *= INTERFACE_NEAR_MULT
                    context_factors.append(f"near_interface (×{INTERFACE_NEAR_MULT})")
            
            # PTM sites
            elif category in {"modified_residue", "lipidation", "glycosylation", "disulfide_bond"}:
                if self._within_window(residue_pos, start_pos, end_pos, NEAR_WINDOW):
                    mult *= 1.2
                    context_factors.append(f"near_PTM (×1.2)")
        
        # 🎯 3) Structural context: burial (RSA)
        if rsa is not None:
            rsa = max(0.0, min(1.0, float(rsa)))
            burial_effect = (1.0 - rsa) * BURIAL_SCALE
            burial_mult = 1.0 + burial_effect
            mult *= burial_mult
            context_factors.append(f"burial_rsa_{rsa:.2f} (×{burial_mult:.2f})")
        
        # 🎯 4) Conservation
        if conservation is not None:
            conservation = max(0.0, min(1.0, float(conservation)))
            if conservation > 0.8:
                cons_effect = ((conservation - 0.8) / 0.2) * CONSERVATION_SCALE
                cons_mult = 1.0 + cons_effect
                mult *= cons_mult
                context_factors.append(f"high_conservation_{conservation:.2f} (×{cons_mult:.2f})")
        
        # 🎯 5) Secondary structure heuristics
        if secondary_structure:
            ss = secondary_structure.upper()[0]
            if ss == 'H':  # Helix
                mult *= 1.15
                context_factors.append("helix_context (×1.15)")
            elif ss == 'E':  # Beta strand
                mult *= 1.05
                context_factors.append("strand_context (×1.05)")
            else:  # Coil/loop
                mult *= 0.95
                context_factors.append("coil_context (×0.95)")
        
        # 🎯 6) Apply caps
        final_mult = max(MIN_MULT, min(MAX_MULT, mult))
        
        print(f"🎯 Proline context factors: {', '.join(context_factors) if context_factors else 'baseline_only'}")
        print(f"🧬 Final proline multiplier: {final_mult:.3f}")
        
        return round(final_mult, 3)

# 🧬 Quick test function
def test_proline_context():
    """Test Nova's proline context system"""
    weighter = ProlineContextWeighter()
    
    demo_features = [
        {"type": "REGION", "description": "triple helix region", "begin": {"position": "1000"}, "end": {"position": "1300"}},
        {"type": "PROPEP", "description": "C-terminal propeptide; cleavage site", "begin": {"position": "1240"}, "end": {"position": "1250"}},
        {"type": "BINDING", "description": "Ca-binding", "begin": {"position": "1400"}, "end": {"position": "1410"}},
    ]
    
    print("\n🧬 Testing Nova's Proline Context System:")
    print("=" * 50)
    
    # Test bulk triple helix (should be downweighted)
    mult1 = weighter.get_proline_multiplier(1150, "P", "A", demo_features, conservation=0.2, rsa=0.6)
    print(f"Bulk triple helix P→A: {mult1}")
    
    # Test near cleavage (should be upweighted)
    mult2 = weighter.get_proline_multiplier(1245, "P", "A", demo_features, conservation=0.9, rsa=0.2)
    print(f"Near cleavage P→A: {mult2}")
    
    # Test normal region (baseline)
    mult3 = weighter.get_proline_multiplier(1500, "P", "A", demo_features, conservation=0.1, rsa=0.7)
    print(f"Normal region P→A: {mult3}")

if __name__ == "__main__":
    test_proline_context()
