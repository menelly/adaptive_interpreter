#!/usr/bin/env python3
"""
🧬 Conservation ML Loader
Runtime system for ML-learned conservation multipliers

Built by Nova & Ace (2025) for DNModeling system
Loads JSON curves and interpolates conservation multipliers
"""

import json
import numpy as np
from typing import Dict, List, Tuple, Optional

class ConservationMLLoader:
    """Load and use ML-learned conservation curves"""
    
    def __init__(self, config_path: str = "conservation_multipliers.json"):
        self.config_path = config_path
        self.config = None
        self.load_config()
    
    def load_config(self):
        """Load conservation multiplier configuration"""
        try:
            with open(self.config_path, "r") as f:
                self.config = json.load(f)
            print(f"✅ Loaded conservation ML config: {self.config['meta']['model_version']}")
            print(f"🎯 Trained on: {self.config['meta']['trained_on']}")
        except FileNotFoundError:
            print(f"⚠️ Conservation ML config not found: {self.config_path}")
            print("🔄 Using fallback hardcoded multipliers")
            self.config = None
        except Exception as e:
            print(f"❌ Error loading conservation config: {e}")
            self.config = None
    
    def get_conservation_multiplier(self, gene_family: str, phylop: float, 
                                  phastcons: float = 0.0, gerp: float = 0.0) -> float:
        """
        Get ML-learned conservation multiplier for gene family + conservation scores
        
        Args:
            gene_family: Gene family classification
            phylop: PhyloP conservation score
            phastcons: PhastCons conservation score  
            gerp: GERP conservation score
            
        Returns:
            Conservation multiplier (1.0 = neutral, >1.0 = boost, <1.0 = penalty)
        """
        
        if self.config is None:
            return self._fallback_multiplier(phylop)
        
        family_config = self.config["families"].get(gene_family)
        if not family_config:
            print(f"⚠️ Unknown gene family: {gene_family}, using GENERAL")
            family_config = self.config["families"].get("GENERAL", {})
        
        if not family_config:
            return self._fallback_multiplier(phylop)
        
        # Get multipliers from each conservation metric
        multipliers = []
        
        # PhyloP curve
        phylop_curve = family_config.get("phylop_curve", [])
        if phylop_curve:
            phylop_mult = self._interpolate_curve(phylop_curve, phylop)
            multipliers.append(phylop_mult)
        
        # PhastCons curve  
        phastcons_curve = family_config.get("phastcons_curve", [])
        if phastcons_curve and phastcons != 0.0:
            phastcons_mult = self._interpolate_curve(phastcons_curve, phastcons)
            multipliers.append(phastcons_mult)
        
        # GERP curve
        gerp_curve = family_config.get("gerp_curve", [])
        if gerp_curve and gerp != 0.0:
            gerp_mult = self._interpolate_curve(gerp_curve, gerp)
            multipliers.append(gerp_mult)
        
        if not multipliers:
            return self._fallback_multiplier(phylop)
        
        # Combine multipliers (geometric mean to avoid extreme values)
        combined = float(np.prod(multipliers) ** (1.0 / len(multipliers)))
        
        # Clamp to reasonable range
        combined = max(0.5, min(3.0, combined))
        
        return combined
    
    def _interpolate_curve(self, curve: List[List[float]], score: float) -> float:
        """Linear interpolate a curve of [x,y] pairs"""
        if not curve:
            return 1.0
        
        # Extract x and y values
        xs, ys = zip(*curve)
        
        # Handle edge cases
        if score <= xs[0]:
            return ys[0]
        if score >= xs[-1]:
            return ys[-1]
        
        # Linear interpolation
        return float(np.interp(score, xs, ys))
    
    def _fallback_multiplier(self, phylop: float) -> float:
        """Fallback to hardcoded multipliers if ML config unavailable"""
        
        # Conservative fallback based on phyloP only
        if phylop >= 7.0:
            return 2.0  # Highly conserved
        elif phylop >= 3.0:
            return 1.5  # Moderately conserved
        elif phylop >= 1.0:
            return 1.2  # Somewhat conserved
        elif phylop >= -1.0:
            return 0.9  # Slightly non-conserved
        else:
            return 0.8  # Not conserved
    
    def get_family_info(self, gene_family: str) -> Dict:
        """Get information about a gene family's conservation curves"""
        if self.config is None:
            return {"error": "No ML config loaded"}
        
        family_config = self.config["families"].get(gene_family, {})
        
        info = {
            "family": gene_family,
            "has_phylop": len(family_config.get("phylop_curve", [])) > 0,
            "has_phastcons": len(family_config.get("phastcons_curve", [])) > 0,
            "has_gerp": len(family_config.get("gerp_curve", [])) > 0,
            "phylop_range": self._get_curve_range(family_config.get("phylop_curve", [])),
            "phastcons_range": self._get_curve_range(family_config.get("phastcons_curve", [])),
            "gerp_range": self._get_curve_range(family_config.get("gerp_curve", []))
        }
        
        return info
    
    def _get_curve_range(self, curve: List[List[float]]) -> Optional[Tuple[float, float]]:
        """Get min/max range of a curve"""
        if not curve:
            return None
        
        xs, ys = zip(*curve)
        return (min(xs), max(xs))
    
    def test_family_multipliers(self, gene_family: str, test_scores: List[float] = None):
        """Test conservation multipliers for a gene family"""
        
        if test_scores is None:
            test_scores = [-3, -1, 0, 1, 3, 5, 7, 9]
        
        print(f"\n🧬 Testing conservation multipliers for {gene_family}")
        print("=" * 60)
        
        family_info = self.get_family_info(gene_family)
        print(f"📊 Family info: {family_info}")
        
        print(f"\n{'PhyloP':<8} {'Multiplier':<12} {'Effect'}")
        print("-" * 30)
        
        for phylop in test_scores:
            mult = self.get_conservation_multiplier(gene_family, phylop)
            
            if mult > 1.2:
                effect = "BOOST"
            elif mult < 0.9:
                effect = "PENALTY"
            else:
                effect = "NEUTRAL"
            
            print(f"{phylop:<8.1f} {mult:<12.3f} {effect}")

def main():
    """Test the conservation ML loader"""
    loader = ConservationMLLoader()
    
    # Test different gene families
    test_families = ["ION_CHANNEL", "COLLAGEN_FIBRILLAR", "TUMOR_SUPPRESSOR", "GENERAL"]
    
    for family in test_families:
        loader.test_family_multipliers(family)
        print()

if __name__ == "__main__":
    main()
