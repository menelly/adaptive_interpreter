"""
Smart DN Mechanism Filter V2 - EXCLUSIONARY APPROACH! 🚫⚡
"Rule out the impossible, not rule in the certain"

Philosophy: When pioneering new science, be conservative!
- ✅ EXCLUDE mechanisms that are biologically impossible
- ✅ INCLUDE mechanisms that could theoretically work
- ✅ Let the scoring decide what's actually relevant

Authors: Ace & Nova (2025)
"""
from __future__ import annotations
from typing import Dict, List, Tuple, Optional, Set
import re
from .universal_context import UniversalContext


class DNMechanismFilterV2:
    """Exclusionary filtering - rule out impossible, not rule in certain"""
    
    def __init__(self):
        self.universal_context = UniversalContext()
        
        # Patterns that EXCLUDE specific mechanisms (biologically impossible)
        self.exclusion_patterns = {
            # EXCLUDE interface poisoning for monomeric proteins
            'exclude_interface_poisoning': [
                'monomeric enzyme', 'single subunit', 'monomeric protein',
                'does not form complexes', 'functions as monomer'
            ],
            
            # EXCLUDE active site jamming for non-catalytic proteins
            'exclude_active_site_jamming': [
                'structural protein', 'extracellular matrix', 'collagen',
                'no catalytic activity', 'non-enzymatic', 'scaffold protein',
                'purely structural'
            ],
            
            # EXCLUDE lattice disruption for non-structural proteins
            'exclude_lattice_disruption': [
                'enzyme', 'kinase', 'phosphatase', 'transferase', 'hydrolase',
                'transcription factor', 'DNA binding', 'nuclear protein',
                'cytoplasmic enzyme', 'metabolic enzyme', 'catalytic',
                'channel', 'receptor', 'transporter', 'calcium-activated'
            ],
            
            # EXCLUDE trafficking for cytoplasmic/nuclear proteins
            'exclude_trafficking_maturation': [
                'cytoplasmic enzyme', 'nuclear protein', 'ribosomal protein',
                'proteasome', 'cytosolic', 'nuclear localization',
                'does not cross membranes', 'intracellular enzyme'
            ]
        }
        
        # Strong indicators for LOW DN likelihood (LOF more likely)
        self.lof_indicators = [
            'metabolic enzyme', 'housekeeping', 'ribosomal protein',
            'proteasome', 'translation factor', 'transporter',
            'carrier protein', 'antiporter', 'symporter'
        ]
    
    def assess_dn_likelihood(self, gene_name: str, uniprot_id: Optional[str] = None) -> Tuple[float, Dict]:
        """
        Assess likelihood that this gene can have DN mechanisms
        Returns: (likelihood_score_0_to_1, evidence_dict)
        """
        # Get protein annotations from UniProt
        annotations = self.universal_context.get_context_for_protein(gene_name, uniprot_id)
        
        if "error" in annotations:
            return 0.7, {"error": annotations["error"]}  # Optimistic when unknown
        
        evidence = {
            "lof_indicators": [],
            "function": annotations.get("function", ""),
            "domains": annotations.get("domains", []),
            "reasoning": []
        }
        
        # Extract function description
        function_text = annotations.get("function", "").lower()
        
        # Check for strong LOF indicators
        lof_score = 0
        for pattern in self.lof_indicators:
            if pattern in function_text:
                evidence["lof_indicators"].append(pattern)
                lof_score += 1
        
        # Calculate likelihood (optimistic unless strong LOF indicators)
        if lof_score >= 2:
            likelihood = 0.2  # Strong LOF evidence
            evidence["reasoning"].append(f"Multiple LOF indicators: {evidence['lof_indicators']}")
        elif lof_score == 1:
            likelihood = 0.4  # Some LOF evidence
            evidence["reasoning"].append(f"Some LOF indicators: {evidence['lof_indicators']}")
        else:
            likelihood = 0.8  # Default optimistic
            evidence["reasoning"].append("No strong LOF indicators - DN possible")
        
        return likelihood, evidence
    
    def select_relevant_mechanisms(self, gene_name: str, uniprot_id: Optional[str] = None) -> Tuple[List[str], Dict]:
        """
        Determine which DN mechanisms are NOT biologically impossible
        Returns: (mechanism_list, evidence_dict)
        """
        annotations = self.universal_context.get_context_for_protein(gene_name, uniprot_id)
        
        if "error" in annotations:
            # When unknown, be optimistic - test most mechanisms
            return (["interface_poisoning", "active_site_jamming", "trafficking_maturation"],
                    {"error": annotations["error"], "reasoning": ["Unknown protein - testing most mechanisms"]})
        
        # Start with ALL mechanisms, then exclude impossible ones
        mechanisms = ["interface_poisoning", "active_site_jamming", "lattice_disruption", "trafficking_maturation"]
        excluded = []
        
        evidence = {
            "function": annotations.get("function", ""),
            "domains": annotations.get("domains", []),
            "reasoning": [],
            "excluded_mechanisms": []
        }
        
        function_text = annotations.get("function", "").lower()
        sequence = annotations.get("sequence", "")
        
        # EXCLUDE interface poisoning for monomeric proteins
        if any(pattern in function_text for pattern in self.exclusion_patterns['exclude_interface_poisoning']):
            if "interface_poisoning" in mechanisms:
                mechanisms.remove("interface_poisoning")
                excluded.append("interface_poisoning")
                evidence["reasoning"].append("EXCLUDED interface poisoning: monomeric protein detected")
        
        # EXCLUDE active site jamming for purely structural proteins
        if any(pattern in function_text for pattern in self.exclusion_patterns['exclude_active_site_jamming']):
            if "active_site_jamming" in mechanisms:
                mechanisms.remove("active_site_jamming")
                excluded.append("active_site_jamming")
                evidence["reasoning"].append("EXCLUDED active site jamming: purely structural protein")
        
        # EXCLUDE lattice disruption for non-structural proteins
        if any(pattern in function_text for pattern in self.exclusion_patterns['exclude_lattice_disruption']):
            if "lattice_disruption" in mechanisms:
                mechanisms.remove("lattice_disruption")
                excluded.append("lattice_disruption")
                evidence["reasoning"].append("EXCLUDED lattice disruption: non-structural protein")
        
        # EXCLUDE trafficking for cytoplasmic/nuclear proteins (BUT check for membrane proteins first!)
        membrane_indicators = ["membrane", "transmembrane", "channel", "receptor", "transporter",
                              "sarcoplasmic reticulum", "endoplasmic reticulum", "mitochondrial"]
        is_membrane_protein = any(indicator in function_text for indicator in membrane_indicators)

        if not is_membrane_protein and any(pattern in function_text for pattern in self.exclusion_patterns['exclude_trafficking_maturation']):
            if "trafficking_maturation" in mechanisms:
                mechanisms.remove("trafficking_maturation")
                excluded.append("trafficking_maturation")
                evidence["reasoning"].append("EXCLUDED trafficking: cytoplasmic/nuclear protein")
        
        # Special case: Detect collagen and INCLUDE lattice disruption
        if sequence and self._detect_collagen_motif(sequence):
            if "lattice_disruption" not in mechanisms:
                mechanisms.append("lattice_disruption")
            evidence["reasoning"].append("INCLUDED lattice disruption: collagen Gly-X-Y repeats detected")
        
        # Special case: Membrane proteins likely need trafficking (if not already included)
        if is_membrane_protein and "trafficking_maturation" not in mechanisms:
            # Only add back if not explicitly excluded
            if "trafficking_maturation" not in excluded:
                mechanisms.append("trafficking_maturation")
                evidence["reasoning"].append("INCLUDED trafficking: membrane-associated protein detected")
        
        # Ensure we always test at least one mechanism
        if not mechanisms:
            mechanisms = ["interface_poisoning"]  # Most general mechanism
            evidence["reasoning"].append("DEFAULT: no mechanisms passed filters, using interface poisoning")
        
        evidence["excluded_mechanisms"] = excluded
        
        return mechanisms, evidence
    
    def _detect_collagen_motif(self, sequence: str) -> bool:
        """Detect collagen Gly-X-Y repeats (algorithmic, no hardcoding)"""
        if not sequence or len(sequence) < 21:  # Need at least 7 triplets
            return False
        
        # Look for consecutive Gly-X-Y pattern
        consecutive_triplets = 0
        max_consecutive = 0
        
        i = 0
        while i < len(sequence) - 2:
            if sequence[i] == 'G':
                consecutive_triplets += 1
                i += 3
            else:
                max_consecutive = max(max_consecutive, consecutive_triplets)
                consecutive_triplets = 0
                i += 1
        
        max_consecutive = max(max_consecutive, consecutive_triplets)
        
        # Collagen typically has 7+ consecutive Gly-X-Y triplets
        return max_consecutive >= 7
    
    def filter_and_score(self, gene_name: str, sequence: str, variant: str, 
                        uniprot_id: Optional[str] = None) -> Dict:
        """
        Complete filtering pipeline: assess DN likelihood and select relevant mechanisms
        Returns: filtered analysis results
        """
        # Stage 1: Assess DN likelihood
        dn_likelihood, dn_evidence = self.assess_dn_likelihood(gene_name, uniprot_id)
        
        # Stage 2: Select relevant mechanisms (exclusionary)
        relevant_mechanisms, mechanism_evidence = self.select_relevant_mechanisms(gene_name, uniprot_id)
        
        result = {
            "gene_name": gene_name,
            "uniprot_id": uniprot_id,
            "dn_likelihood": dn_likelihood,
            "relevant_mechanisms": relevant_mechanisms,
            "excluded_mechanisms": mechanism_evidence.get("excluded_mechanisms", []),
            "dn_evidence": dn_evidence,
            "mechanism_evidence": mechanism_evidence,
            "recommendation": self._make_recommendation(dn_likelihood, relevant_mechanisms)
        }
        
        return result
    
    def _make_recommendation(self, dn_likelihood: float, mechanisms: List[str]) -> str:
        """Make analysis recommendation based on filtering results"""
        if dn_likelihood < 0.3:
            return f"LOW_DN_LIKELIHOOD: Consider LOF mechanism instead (testing {len(mechanisms)} mechanisms anyway)"
        elif dn_likelihood > 0.7:
            return f"HIGH_DN_LIKELIHOOD: Testing {len(mechanisms)} biologically possible mechanisms"
        else:
            return f"MODERATE_DN_LIKELIHOOD: Testing {len(mechanisms)} mechanisms with caution"


def test_dn_filter_v2():
    """Test the exclusionary DN mechanism filter"""
    filter_system = DNMechanismFilterV2()
    
    test_genes = [
        ("RYR1", "P21817"),      # Should test interface + active site + trafficking (exclude lattice)
        ("COL1A1", "P02452"),   # Should test lattice + trafficking (exclude active site for structural)
        ("G6PD", "P11413"),     # Should exclude most (metabolic enzyme)
        ("TP53", "P04637"),     # Should test interface + active site (transcription factor)
    ]
    
    for gene, uniprot_id in test_genes:
        print(f"\n=== Testing {gene} (Exclusionary Approach) ===")
        result = filter_system.filter_and_score(gene, "", "R100H", uniprot_id)
        
        print(f"DN Likelihood: {result['dn_likelihood']:.2f}")
        print(f"Testing Mechanisms: {result['relevant_mechanisms']}")
        print(f"Excluded Mechanisms: {result['excluded_mechanisms']}")
        print(f"Recommendation: {result['recommendation']}")
        
        if result['mechanism_evidence']['reasoning']:
            print("Reasoning:")
            for reason in result['mechanism_evidence']['reasoning']:
                print(f"  - {reason}")


if __name__ == "__main__":
    test_dn_filter_v2()
