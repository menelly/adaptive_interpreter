#!/usr/bin/env python3
"""
Functional Domain Weighting System
Based on Nova's biological intelligence approach

Replaces location-based domain weighting with function-based weighting
using UniProt annotations to determine biological importance.
"""

import json
import re
from typing import Dict, List, Optional, Tuple

# ---------------------------
# DOMAIN IMPORTANCE HIERARCHY
# ---------------------------
DOMAIN_WEIGHTS = {
    "cleavage_site": 1.5,
    "active_site": 1.5,
    "binding_site": 1.5,
    "processing_domain": 1.3,
    "catalytic_domain": 1.3,
    "structural_domain": 1.2,
    "triple_helix": 1.2,
    "mature_chain": 1.0,
    "functional_region": 1.0,
    "linker_region": 0.8,
    "disordered_region": 0.8,
    "signal_peptide": 0.5,
    "non_functional": 0.5,
}

class FunctionalDomainWeighter:
    """
    Biological intelligence-based domain weighting system.
    
    Uses UniProt functional annotations to determine the biological
    importance of protein regions, replacing simplistic location-based
    weighting with function-based intelligence.
    """
    
    def __init__(self):
        self.weights = DOMAIN_WEIGHTS.copy()
        
    def classify_domain(self, uniprot_feature: Dict) -> str:
        """
        Classify a UniProt feature into functional categories.
        
        Args:
            uniprot_feature: UniProt feature dict with 'type' and 'description'
            
        Returns:
            str: Functional category for weighting
        """
        desc = uniprot_feature.get("description", "").lower()
        ftype = uniprot_feature.get("type", "").upper()
        
        # Critical functional sites (1.5x weight)
        if any(keyword in desc for keyword in ["cleavage", "processing", "propeptide"]):
            return "cleavage_site"
        if ftype in ["ACT_SITE", "BINDING", "SITE"]:
            return "active_site"
        if any(keyword in desc for keyword in ["binding", "active site"]):
            return "binding_site"
            
        # High importance domains (1.3x weight)
        if any(keyword in desc for keyword in ["catalytic", "enzyme", "processing"]):
            return "catalytic_domain"
        if ftype == "PROPEP" and "processing" in desc:
            return "processing_domain"
            
        # Moderate importance (1.2x weight)
        if any(keyword in desc for keyword in ["helix", "triple", "structural"]):
            return "triple_helix"
        if ftype == "DOMAIN" or "structural" in desc:
            return "structural_domain"
            
        # Normal importance (1.0x weight)
        if ftype == "CHAIN":
            return "mature_chain"
        if any(keyword in desc for keyword in ["functional", "region"]):
            return "functional_region"
            
        # Low importance (0.8x weight)
        if any(keyword in desc for keyword in ["disordered", "low complexity", "linker"]):
            return "disordered_region"
        if "linker" in desc:
            return "linker_region"
            
        # Minimal importance (0.5x weight)
        if ftype == "SIGNAL" or "signal" in desc:
            return "signal_peptide"
        if any(keyword in desc for keyword in ["non-functional", "spacer"]):
            return "non_functional"
            
        # Default fallback
        return "functional_region"
    
    def get_domain_weight(self, uniprot_feature: Dict) -> float:
        """
        Return the domain weight for a given UniProt feature.
        
        Args:
            uniprot_feature: UniProt feature dict
            
        Returns:
            float: Weight multiplier for this domain
        """
        category = self.classify_domain(uniprot_feature)
        return self.weights.get(category, 1.0)
    
    def weight_variant_position(self, variant_pos: int, features: List[Dict]) -> float:
        """
        Given a variant position and UniProt features, return the highest weight.
        
        Args:
            variant_pos: 1-based amino acid position
            features: List of UniProt feature dicts with position ranges
            
        Returns:
            float: Weight multiplier for this position
        """
        applicable_weights = []
        
        for feat in features:
            try:
                # Handle different UniProt JSON formats
                if "begin" in feat and "end" in feat:
                    start = int(feat["begin"].get("position", 0))
                    end = int(feat["end"].get("position", start))
                elif "start" in feat and "end" in feat:
                    start = int(feat["start"])
                    end = int(feat["end"])
                else:
                    continue
                    
            except (ValueError, TypeError, KeyError):
                continue
            
            # Check if variant falls within this feature
            if start <= variant_pos <= end:
                weight = self.get_domain_weight(feat)
                applicable_weights.append(weight)
                
                # Debug info
                category = self.classify_domain(feat)
                print(f"🎯 Functional domain at pos {variant_pos}: {category} (weight: {weight})")
        
        if not applicable_weights:
            print(f"🔍 No functional domains found at pos {variant_pos}, using default weight: 1.0")
            return 1.0  # default if no domain found
            
        # Use maximum weight (strongest biological importance wins)
        max_weight = max(applicable_weights)
        print(f"🧬 Final functional weight for pos {variant_pos}: {max_weight}")
        return max_weight
    
    def get_feature_summary(self, features: List[Dict]) -> Dict[str, int]:
        """
        Get a summary of functional categories in the protein.
        
        Args:
            features: List of UniProt features
            
        Returns:
            Dict mapping categories to counts
        """
        summary = {}
        for feat in features:
            category = self.classify_domain(feat)
            summary[category] = summary.get(category, 0) + 1
        return summary

# ---------------------------
# CONVENIENCE FUNCTIONS
# ---------------------------
def create_functional_weighter() -> FunctionalDomainWeighter:
    """Create a new functional domain weighter instance."""
    return FunctionalDomainWeighter()

def weight_variant_functionally(variant_pos: int, uniprot_features: List[Dict]) -> float:
    """
    Convenience function to weight a variant based on functional domains.
    
    Args:
        variant_pos: 1-based amino acid position
        uniprot_features: List of UniProt feature annotations
        
    Returns:
        float: Functional weight multiplier
    """
    weighter = create_functional_weighter()
    return weighter.weight_variant_position(variant_pos, uniprot_features)

if __name__ == "__main__":
    # Demo usage with COL1A1-like features
    demo_features = [
        {
            "type": "PROPEP", 
            "description": "C-terminal propeptide; cleavage site",
            "begin": {"position": "1219"}, 
            "end": {"position": "1464"}
        },
        {
            "type": "CHAIN", 
            "description": "Mature collagen chain",
            "begin": {"position": "162"}, 
            "end": {"position": "1218"}
        },
        {
            "type": "REGION", 
            "description": "Triple-helical region",
            "begin": {"position": "162"}, 
            "end": {"position": "1218"}
        }
    ]
    
    weighter = create_functional_weighter()
    
    print("🧬 FUNCTIONAL DOMAIN WEIGHTING DEMO")
    print("=" * 50)
    
    # Test positions that were problematic in COL1A1
    test_positions = [1256, 1441, 1289, 500, 1000]
    
    for pos in test_positions:
        print(f"\n🔍 Testing position {pos}:")
        weight = weighter.weight_variant_position(pos, demo_features)
        print(f"   Final weight: {weight}")
    
    print(f"\n📊 Feature summary:")
    summary = weighter.get_feature_summary(demo_features)
    for category, count in summary.items():
        print(f"   {category}: {count} features")
