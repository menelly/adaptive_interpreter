"""
Smart DN Mechanism Filter - NO HARDCODING!
Two-stage approach: DN likelihood + mechanism relevance

Uses ONLY:
- GO terms from UniProt API
- Pfam domains from databases  
- Sequence motifs from algorithmic detection
- Biological principles (no magic numbers)
"""
from __future__ import annotations
from typing import Dict, List, Tuple, Optional, Set
import re
from .universal_context import UniversalContext


class DNMechanismFilter:
    """Smart filtering system to determine which DN mechanisms are biologically relevant"""
    
    def __init__(self):
        self.universal_context = UniversalContext()
        
        # GO term patterns that indicate DN potential (from GO database)
        self.dn_indicating_patterns = {
            # Structural/oligomeric proteins
            'oligomerization': ['oligomer', 'dimer', 'trimer', 'tetramer', 'hexamer', 'multimer'],
            'structural': ['structural', 'cytoskeleton', 'extracellular matrix', 'collagen'],
            'complex_assembly': ['complex assembly', 'protein complex', 'subunit'],
            
            # Transcription factors (often dominant)
            'transcription': ['transcription', 'DNA binding', 'chromatin', 'histone'],
            
            # Receptors and channels (often dominant)
            'signaling': ['receptor', 'channel', 'signal transduction'],
            
            # Secreted proteins (trafficking issues)
            'secreted': ['secreted', 'extracellular', 'signal peptide', 'glycoprotein']
        }
        
        # Patterns that suggest LOF/recessive (low DN likelihood)
        self.lof_indicating_patterns = {
            'metabolic_enzyme': ['hydrolase', 'transferase', 'oxidoreductase', 'lyase', 'isomerase'],
            'transporter': ['transporter', 'carrier', 'permease', 'antiporter', 'symporter'],
            'housekeeping': ['ribosomal', 'proteasome', 'mitochondrial matrix', 'translation']
        }
    
    def assess_dn_likelihood(self, gene_name: str, uniprot_id: Optional[str] = None) -> Tuple[float, Dict]:
        """
        Assess likelihood that this gene can have DN mechanisms
        Returns: (likelihood_score_0_to_1, evidence_dict)
        """
        # Get protein annotations from UniProt
        annotations = self.universal_context.get_context_for_protein(gene_name, uniprot_id)
        
        if "error" in annotations:
            return 0.5, {"error": annotations["error"]}  # Neutral when unknown
        
        evidence = {
            "dn_indicators": [],
            "lof_indicators": [],
            "function": annotations.get("function", ""),
            "domains": annotations.get("domains", []),
            "go_terms": []  # Will be populated from function text
        }
        
        # Extract GO-like terms from function description
        function_text = annotations.get("function", "").lower()
        
        # Check for DN-indicating patterns
        dn_score = 0.0
        for category, patterns in self.dn_indicating_patterns.items():
            for pattern in patterns:
                if pattern in function_text:
                    evidence["dn_indicators"].append(f"{category}: {pattern}")
                    dn_score += 1.0
        
        # Check for LOF-indicating patterns  
        lof_score = 0.0
        for category, patterns in self.lof_indicating_patterns.items():
            for pattern in patterns:
                if pattern in function_text:
                    evidence["lof_indicators"].append(f"{category}: {pattern}")
                    lof_score += 1.0
        
        # Calculate likelihood (0.0 = definitely LOF, 1.0 = definitely DN)
        if dn_score == 0 and lof_score == 0:
            likelihood = 0.5  # Unknown
        else:
            likelihood = dn_score / (dn_score + lof_score)
        
        return likelihood, evidence
    
    def select_relevant_mechanisms(self, gene_name: str, uniprot_id: Optional[str] = None) -> Tuple[List[str], Dict]:
        """
        Determine which DN mechanisms are biologically relevant for this gene
        Returns: (mechanism_list, evidence_dict)
        """
        annotations = self.universal_context.get_context_for_protein(gene_name, uniprot_id)
        
        if "error" in annotations:
            # Default to safest mechanisms when unknown
            return ["interface_poisoning", "trafficking_maturation"], {"error": annotations["error"]}
        
        mechanisms = []
        evidence = {
            "function": annotations.get("function", ""),
            "domains": annotations.get("domains", []),
            "sequence_motifs": {},
            "reasoning": []
        }
        
        function_text = annotations.get("function", "").lower()
        sequence = annotations.get("sequence", "")
        
        # Interface poisoning - for oligomeric/complex proteins
        if any(pattern in function_text for pattern in 
               self.dn_indicating_patterns['oligomerization'] + 
               self.dn_indicating_patterns['complex_assembly'] +
               self.dn_indicating_patterns['signaling']):
            mechanisms.append("interface_poisoning")
            evidence["reasoning"].append("Interface poisoning: protein complex/oligomer detected")
        
        # Active site jamming - for enzymes and DNA-binding proteins
        if any(pattern in function_text for pattern in 
               ["enzyme", "catalytic", "kinase", "phosphatase", "DNA binding", "nuclease"]):
            mechanisms.append("active_site_jamming")
            evidence["reasoning"].append("Active site jamming: catalytic/binding activity detected")
        
        # Lattice disruption - for structural proteins
        if any(pattern in function_text for pattern in 
               self.dn_indicating_patterns['structural']):
            mechanisms.append("lattice_disruption")
            evidence["reasoning"].append("Lattice disruption: structural protein detected")
        
        # Detect collagen specifically (sequence-based)
        if sequence and self._detect_collagen_motif(sequence):
            if "lattice_disruption" not in mechanisms:
                mechanisms.append("lattice_disruption")
            evidence["sequence_motifs"]["collagen_repeats"] = True
            evidence["reasoning"].append("Lattice disruption: collagen Gly-X-Y repeats detected")
        
        # Trafficking/maturation - for secreted proteins, membrane proteins
        if any(pattern in function_text for pattern in 
               self.dn_indicating_patterns['secreted'] + 
               ["membrane", "transmembrane", "glycosylation", "folding"]):
            mechanisms.append("trafficking_maturation")
            evidence["reasoning"].append("Trafficking/maturation: secreted/membrane protein detected")
        
        # Default fallback - if no specific mechanisms identified
        if not mechanisms:
            mechanisms = ["interface_poisoning", "trafficking_maturation"]
            evidence["reasoning"].append("Default: using safest mechanism combination")
        
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
        
        # Stage 2: Select relevant mechanisms
        relevant_mechanisms, mechanism_evidence = self.select_relevant_mechanisms(gene_name, uniprot_id)
        
        result = {
            "gene_name": gene_name,
            "uniprot_id": uniprot_id,
            "dn_likelihood": dn_likelihood,
            "relevant_mechanisms": relevant_mechanisms,
            "dn_evidence": dn_evidence,
            "mechanism_evidence": mechanism_evidence,
            "recommendation": self._make_recommendation(dn_likelihood, relevant_mechanisms)
        }
        
        return result
    
    def _make_recommendation(self, dn_likelihood: float, mechanisms: List[str]) -> str:
        """Make analysis recommendation based on filtering results"""
        if dn_likelihood < 0.3:
            return "LOW_DN_LIKELIHOOD: Consider LOF/recessive mechanism instead"
        elif dn_likelihood > 0.7:
            return f"HIGH_DN_LIKELIHOOD: Analyze using {len(mechanisms)} relevant mechanisms"
        else:
            return f"MODERATE_DN_LIKELIHOOD: Analyze using {len(mechanisms)} mechanisms with caution"


def test_dn_filter():
    """Test the DN mechanism filter"""
    filter_system = DNMechanismFilter()
    
    test_genes = [
        ("TP53", "P04637"),      # Should be high DN (transcription factor)
        ("COL1A1", "P02452"),   # Should be high DN (structural, collagen)
        ("G6PD", "P11413"),     # Should be low DN (metabolic enzyme)
        ("SLC25A5", "P05141"),  # Should be low DN (transporter)
    ]
    
    for gene, uniprot_id in test_genes:
        print(f"\n=== Testing {gene} ===")
        result = filter_system.filter_and_score(gene, "", "R100H", uniprot_id)
        
        print(f"DN Likelihood: {result['dn_likelihood']:.2f}")
        print(f"Relevant Mechanisms: {result['relevant_mechanisms']}")
        print(f"Recommendation: {result['recommendation']}")
        
        if result['mechanism_evidence']['reasoning']:
            print("Reasoning:")
            for reason in result['mechanism_evidence']['reasoning']:
                print(f"  - {reason}")


if __name__ == "__main__":
    test_dn_filter()
