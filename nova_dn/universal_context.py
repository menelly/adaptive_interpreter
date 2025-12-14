#!/usr/bin/env python3
"""
Universal Context Integration for Nova DN Analyzer
Replaces hardcoded annotations with automated UniProt extraction

🐙 2025-12-10: Added InterPro domain integration for REAL domain-based interface scoring
   - InterPro domains are the authoritative source (SPRY, BRCT, catalytic sites, etc.)
   - Predicted hydrophobic patches are fallback only
   - Domain TYPE determines interface_likelihood score
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from data_processing.universal_protein_annotator import UniversalProteinAnnotator
import json
from typing import Dict, Optional, List

# 🧬 Domain type weights for interface_likelihood scoring
# These are UNIVERSAL InterPro entry types - not gene-specific!
# More specific annotations → higher confidence in functional importance
DOMAIN_TYPE_WEIGHTS = {
    # HIGH - Specific functional domains (InterPro "domain" entries)
    "domain": 0.7,           # Specific domains (zinc finger, kinase, etc.)
    "repeat": 0.6,           # Repeat domains (often structural but important)

    # MEDIUM - Broader classifications (less specific about function)
    "family": 0.4,           # Gene family membership
    "homologous_superfamily": 0.3,  # Broad structural classification

    # DEFAULT for unknown types
    "default": 0.5
}

# 🎯 UNIVERSAL functional keywords that indicate HIGH importance across ALL genes
# These are biological function terms, NOT gene-specific domain names!
HIGH_IMPORTANCE_KEYWORDS = {
    # Enzymatic/catalytic functions
    "kinase", "phosphatase", "protease", "catalytic", "active site",
    "transferase", "hydrolase", "oxidoreductase", "ligase",
    # Binding functions
    "binding", "zinc finger", "ring-type",
    # DNA/RNA interaction
    "dna-binding", "rna-binding", "helix-turn-helix", "homeobox",
    # Channel/transport
    "pore", "channel", "transporter", "transmembrane",
    # Structural importance
    "substrate", "cofactor", "atp", "gtp", "nad",
}

# 🔻 UNIVERSAL features that indicate LOWER constraint (more tolerant of variation)
# These are structural/sequence features, NOT gene-specific domain names!
LOW_IMPORTANCE_KEYWORDS = {
    # Disordered/flexible regions - universally more tolerant
    "disordered", "low complexity", "compositionally biased",
    "proline-rich", "serine-rich", "glutamine-rich", "alanine-rich",
    # Generic/broad annotations
    "uncharacterized", "unknown function",
}

# NOTE: We intentionally do NOT hardcode gene-specific domains like "dapin", "spry", "brct"
# Instead, we rely on:
# 1. InterPro domain TYPE hierarchy (domain > family > superfamily)
# 2. Universal functional keywords
# 3. Conservation data (when available) to discriminate tolerant vs. constrained positions
# 4. Future: gnomAD MPC/MCR regional missense constraint scores


class UniversalContext:
    def __init__(self):
        self.annotator = UniversalProteinAnnotator()
        self.cache = {}
        self.interpro_cache = {}  # Separate cache for InterPro data

    def _load_interpro_domains(self, uniprot_id: str) -> Optional[Dict]:
        """
        Load REAL domain boundaries from InterPro cache.

        Returns dict with 'domains' list, or None if no cache exists.
        """
        if uniprot_id in self.interpro_cache:
            return self.interpro_cache[uniprot_id]

        # Try multiple possible cache locations
        possible_paths = [
            f"protein_annotations_cache/{uniprot_id}_interpro_domains.json",
            f"AdaptiveInterpreter/protein_annotations_cache/{uniprot_id}_interpro_domains.json",
            f"/home/Ace/AdaptiveInterpreter/protein_annotations_cache/{uniprot_id}_interpro_domains.json"
        ]

        for path in possible_paths:
            if os.path.exists(path):
                try:
                    with open(path, 'r') as f:
                        data = json.load(f)
                    self.interpro_cache[uniprot_id] = data
                    return data
                except Exception as e:
                    print(f"   ⚠️ Failed to load InterPro cache: {e}")

        # No cache found
        self.interpro_cache[uniprot_id] = None
        return None

    def _get_interface_likelihood_from_interpro(self, uniprot_id: str, position: int, protein_length: int = None) -> tuple[float, str]:
        """
        Get interface_likelihood score based on InterPro domain at position.

        Uses UNIVERSAL criteria only - no gene-specific hardcoding!
        - Domain TYPE hierarchy (domain > family > superfamily)
        - Universal functional keywords (kinase, binding, etc.)
        - Universal low-constraint keywords (disordered, low complexity, etc.)
        - Position-based heuristics (N-terminal regions often less constrained)

        Returns:
            (score, domain_description) - score is 0.0 if position not in any domain
        """
        interpro_data = self._load_interpro_domains(uniprot_id)

        if not interpro_data or not interpro_data.get('domains'):
            return (0.0, "")

        # Get protein length from InterPro data if not provided
        if protein_length is None:
            protein_length = interpro_data.get('length', 0)

        # Find all domains covering this position
        covering_domains = []
        for domain in interpro_data['domains']:
            start = domain.get('start', 0)
            end = domain.get('end', 0)
            if start <= position <= end:
                covering_domains.append(domain)

        if not covering_domains:
            return (0.0, "no_domain_coverage")

        # Score based on the MOST functionally important domain at this position
        best_score = 0.0
        best_domain = ""
        best_domain_size = 0
        has_high_importance = False
        all_low_importance = True  # Assume all are low until we find one that isn't

        for domain in covering_domains:
            domain_type = domain.get('type', 'default').lower()
            description = domain.get('description', '').lower()
            domain_size = domain.get('end', 0) - domain.get('start', 0) + 1

            # Start with base score from domain TYPE (universal InterPro hierarchy)
            base_score = DOMAIN_TYPE_WEIGHTS.get(domain_type, DOMAIN_TYPE_WEIGHTS['default'])

            # Check for UNIVERSAL high-importance functional keywords
            for keyword in HIGH_IMPORTANCE_KEYWORDS:
                if keyword in description:
                    base_score = max(base_score, 0.70)
                    has_high_importance = True
                    break

            # Check if this domain has UNIVERSAL low-importance features
            is_this_low = False
            for keyword in LOW_IMPORTANCE_KEYWORDS:
                if keyword in description:
                    is_this_low = True
                    break

            if not is_this_low:
                all_low_importance = False

            if base_score > best_score:
                best_score = base_score
                best_domain = domain.get('description', domain_type)
                best_domain_size = domain_size

        # KEY LOGIC: If ALL covering domains are low-importance (and none are high),
        # cap the score to indicate this position likely tolerates variation
        if all_low_importance and not has_high_importance:
            best_score = min(best_score, 0.40)

        # 🧬 UNIVERSAL DOMAIN SIZE HEURISTIC: Larger domains tend to be more constrained
        # This is based on the principle that larger functional domains have more
        # critical residues and are under stronger purifying selection.
        # Small domains (<100 aa) get a reduction, large domains (>150 aa) get a boost.
        if best_domain_size > 0 and not has_high_importance:
            if best_domain_size < 100:
                # Small domain - apply reduction (less likely to be critical)
                size_factor = 0.85
                best_domain += f" (small domain: {best_domain_size}aa)"
            elif best_domain_size > 150:
                # Large domain - apply boost (more likely to be critical)
                size_factor = 1.10
                best_domain += f" (large domain: {best_domain_size}aa)"
            else:
                size_factor = 1.0
            best_score *= size_factor

        # 🧬 UNIVERSAL POSITION HEURISTIC: N-terminal regions are often less constrained
        # This is a well-established biological principle - signal peptides, pro-domains,
        # and regulatory regions at the N-terminus often tolerate more variation.
        # Apply a modest reduction for positions in the first 15% of the protein.
        if protein_length > 0 and not has_high_importance:
            relative_position = position / protein_length
            if relative_position < 0.15:
                # N-terminal region - apply modest reduction
                n_term_factor = 0.90  # 10% reduction (reduced from 15% to avoid double-penalty)
                best_score *= n_term_factor
                if "(N-terminal region)" not in best_domain:
                    best_domain += " (N-terminal region)"

        return (best_score, best_domain)
    
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
        # 🐙 2025-12-10: Use InterPro domains as AUTHORITATIVE source
        # Fall back to predicted hydrophobic patches only if no InterPro data
        # Try both "length" and "sequence_length" keys (different sources use different names)
        protein_length = annotations.get("length") or annotations.get("sequence_length") or 0
        # Also try to get from sequence if available
        if protein_length == 0 and annotations.get("sequence"):
            protein_length = len(annotations.get("sequence", ""))
        interpro_score, interpro_domain = self._get_interface_likelihood_from_interpro(
            uniprot_id or annotations.get("uniprot_id", ""),
            position,
            protein_length
        )

        if interpro_score > 0:
            # InterPro domain found at this position - use domain-weighted score
            context["interface_likelihood"] = interpro_score
            context["_interpro_domain"] = interpro_domain  # For debugging
        elif interpro_domain == "no_domain_coverage":
            # InterPro data exists but position is NOT in any domain
            # This is GOOD - means position is likely in a linker/unstructured region
            context["interface_likelihood"] = 0.0
        else:
            # No InterPro data - fall back to predicted interface regions
            interface_regions = annotations.get("interface_regions", [])
            if interface_regions:
                # Check if it's old format (list of 2 ints) or new format (list of lists)
                if isinstance(interface_regions[0], list):
                    # New format: list of [start, end] ranges
                    for region in interface_regions:
                        if len(region) >= 2 and region[0] <= position <= region[1]:
                            context["interface_likelihood"] = 0.6
                            break
                elif len(interface_regions) >= 2:
                    # Old format: [start, end] as two separate elements
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
