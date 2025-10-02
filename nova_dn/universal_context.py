#!/usr/bin/env python3
"""
Universal Context Integration for Nova DN Analyzer
Replaces hardcoded annotations with automated UniProt extraction
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from data_processing.universal_protein_annotator import UniversalProteinAnnotator
import json
from typing import Dict, Optional

class UniversalContext:
    def __init__(self):
        self.annotator = UniversalProteinAnnotator()
        self.cache = {}
    
    def get_context_for_protein(self, gene_name: str, uniprot_id: Optional[str] = None) -> Dict:
        """Get context for any protein - NO HARDCODING!"""
        cache_key = uniprot_id or gene_name
        
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        print(f"🔍 Auto-annotating {gene_name}...")
        annotations = self.annotator.annotate_protein(gene_name, uniprot_id)
        
        if "error" in annotations:
            print(f"⚠️ Could not annotate {gene_name}: {annotations['error']}")
            # Return minimal context
            annotations = {"uniprot_id": uniprot_id or "unknown"}
        
        self.cache[cache_key] = annotations
        return annotations
    
    def build_position_context(self, gene_name: str, position: int, uniprot_id: Optional[str] = None) -> Dict:
        """Build context flags for a specific position"""
        annotations = self.get_context_for_protein(gene_name, uniprot_id)
        
        context = {
            "active_site_proximity": 0.0,
            "interface_likelihood": 0.0,
            "flexible_loop": 0.0,
            "transmembrane_region": 0.0,
            "collagen_gly_site": 0.0,
            "disulfide_proximity": 0.0
        }
        
        # Active site proximity
        active_sites = annotations.get("known_active_or_binding_sites", [])
        if active_sites:
            min_dist = min(abs(position - site) for site in active_sites)
            if min_dist == 0:
                context["active_site_proximity"] = 1.0
            elif min_dist <= 5:
                context["active_site_proximity"] = 0.8
            elif min_dist <= 10:
                context["active_site_proximity"] = 0.4
        
        # DNA contact sites (special case of active sites)
        dna_sites = annotations.get("dna_contact_sites", [])
        if dna_sites:
            min_dist = min(abs(position - site) for site in dna_sites)
            if min_dist == 0:
                context["active_site_proximity"] = 1.0
            elif min_dist <= 3:
                context["active_site_proximity"] = 0.8
        
        # Interface likelihood
        interface_regions = annotations.get("interface_regions", [])
        if len(interface_regions) >= 2:
            start, end = interface_regions[0], interface_regions[1]
            if start <= position <= end:
                context["interface_likelihood"] = 0.6
        
        # Flexible loops
        flexible_loops = annotations.get("flexible_loops", [])
        if position in flexible_loops:
            context["flexible_loop"] = 1.0
        
        # Transmembrane regions
        tm_domain = annotations.get("transmembrane_domain", [])
        if len(tm_domain) >= 2:
            start, end = tm_domain[0], tm_domain[1]
            if start <= position <= end:
                context["transmembrane_region"] = 1.0
        
        # Collagen Gly sites
        collagen_repeats = annotations.get("collagen_repeats", {})
        if collagen_repeats:
            start = collagen_repeats.get("start", 0)
            end = collagen_repeats.get("end", 0)
            if start <= position <= end:
                # Check if it's actually a Gly position in the repeat
                relative_pos = (position - start) % 3
                if relative_pos == 0:  # First position of triplet should be Gly
                    context["collagen_gly_site"] = 1.0
        
        # Disulfide proximity
        disulfide_pairs = annotations.get("disulfide_pairs", [])
        for pair in disulfide_pairs:
            if len(pair) >= 2:
                if position in pair:
                    context["disulfide_proximity"] = 1.0
                    break
                else:
                    min_dist = min(abs(position - p) for p in pair)
                    if min_dist <= 5:
                        context["disulfide_proximity"] = 0.5
                        break
        
        return context


def test_universal_context():
    """Test the universal context system"""
    context = UniversalContext()
    
    test_cases = [
        ("TP53", 273, "P04637"),  # Known DNA contact
        ("TP53", 72, "P04637"),   # Flexible loop
        ("COL1A1", 1076, "P02452"), # Collagen Gly
        ("BRCA1", 1847, "P38398"),  # BRCT domain
        ("CFTR", 508, "P13569")     # NBD1 domain
    ]
    
    print("=== UNIVERSAL CONTEXT TESTING ===")
    for gene, pos, uniprot_id in test_cases:
        print(f"\n--- {gene} position {pos} ---")
        ctx = context.build_position_context(gene, pos, uniprot_id)
        
        for key, value in ctx.items():
            if value > 0:
                print(f"  {key}: {value:.2f}")
        
        if all(v == 0 for v in ctx.values()):
            print("  No special context detected")


if __name__ == "__main__":
    test_universal_context()
