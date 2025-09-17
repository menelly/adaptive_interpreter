#!/usr/bin/env python3
"""
🧬 BIOLOGICAL ROUTER - Smart Mechanism-Based Analyzer Selection
Routes variants to appropriate analyzers based on biological principles!

Built by Ace & Nova (2025) for intelligent genetics analysis

Logic:
- Easy classifications → Direct routing (CFTR → LOF only)
- Uncertain cases → Run all analyzers (conservative approach)
- High confidence → Skip unnecessary analyses
- Low confidence → Full analysis for safety

This makes the system faster AND more biologically principled!
"""

from typing import Dict, List, Tuple, Optional, Set
from pathlib import Path
import sys

# Import our existing mechanism filter
sys.path.append(str(Path(__file__).parent / "nova_dn"))
from nova_dn.dn_mechanism_filter import DNMechanismFilter
from nova_dn.universal_context import UniversalContext


class BiologicalRouter:
    """Smart router that determines which analyzers to run based on gene biology"""
    
    def __init__(self):
        self.dn_filter = DNMechanismFilter()
        self.universal_context = UniversalContext()
        
        # 🚀 NO MORE HARDCODED GENES! 🚀
        # Trust the universal GO term analysis system!
        # Every gene gets analyzed the same way - no special cases!
        self.high_confidence_genes = {}
        
        # GO term patterns for classification (when gene not in high-confidence list)
        self.go_classification_patterns = {
            'LOF_INDICATORS': [
                # Enzymatic activities (loss of function)
                'hydrolase activity', 'transferase activity', 'oxidoreductase activity',
                'lyase activity', 'isomerase activity', 'ligase activity',
                'transporter activity', 'carrier activity', 'permease activity',
                'metabolic process', 'catabolic process', 'biosynthetic process',
                'ribosomal protein', 'proteasome', 'translation factor',
                # Tumor suppressors (loss of function)
                'tumor suppressor', 'tumour suppressor', 'negative regulator',
                'cell cycle arrest', 'apoptosis', 'dna repair',
                # Housekeeping functions
                'housekeeping', 'maintenance', 'homeostasis'
            ],
            'GOF_INDICATORS': [
                # Oncogenes and growth promotion
                'oncogene', 'proto-oncogene', 'growth factor receptor',
                'tyrosine kinase activity', 'serine/threonine kinase activity',
                'signal transduction', 'cell proliferation', 'cell cycle progression',
                # Activating functions
                'positive regulator', 'activator', 'stimulator'
            ],
            'DN_INDICATORS': [
                # Structural proteins
                'structural constituent', 'extracellular matrix',
                'collagen', 'fibrillin', 'elastin', 'cytoskeleton',
                # Complex-forming proteins
                'transcription factor', 'DNA binding', 'chromatin binding',
                'protein complex assembly', 'oligomerization',
                'tetramer', 'dimer', 'multimer', 'quaternary structure'
            ]
        }
    
    def route_variant(self, gene: str, variant: str, variant_type: str = 'missense', 
                     uniprot_id: Optional[str] = None) -> Dict:
        """
        Determine which analyzers to run for this variant
        
        Args:
            gene: Gene symbol (e.g., 'TP53')
            variant: Variant notation (e.g., 'p.R273H')
            variant_type: Type of variant ('missense', 'frameshift', 'nonsense', etc.)
            uniprot_id: Optional UniProt ID
            
        Returns:
            {
                'analyzers_to_run': ['DN', 'LOF', 'GOF'],
                'routing_strategy': 'HIGH_CONFIDENCE' | 'GO_BASED' | 'CONSERVATIVE',
                'confidence': 0.95,
                'rationale': 'CFTR is ion channel with recessive inheritance',
                'gene_classification': 'LOF_ONLY'
            }
        """
        
        # Handle non-missense variants first
        if variant_type in ['frameshift', 'nonsense', 'splice']:
            return {
                'analyzers_to_run': ['LOF'],
                'routing_strategy': 'VARIANT_TYPE_BASED',
                'confidence': 0.99,
                'rationale': f'{variant_type} variants are almost always loss-of-function',
                'gene_classification': 'LOF_BY_VARIANT_TYPE'
            }
        
        # Check high-confidence gene classifications first
        for classification, genes in self.high_confidence_genes.items():
            if gene in genes:
                analyzers = self._get_analyzers_for_classification(classification)
                return {
                    'analyzers_to_run': analyzers,
                    'routing_strategy': 'HIGH_CONFIDENCE',
                    'confidence': 0.95,
                    'rationale': genes[gene],
                    'gene_classification': classification
                }
        
        # Use GO terms for classification
        go_result = self._classify_by_go_terms(gene, uniprot_id)
        if go_result['confidence'] >= 0.65:  # Lowered threshold to catch more tumor suppressors
            # 🎯 REN'S BRILLIANT IDEA: For medium confidence, run primary + backup analyzers!
            if go_result['confidence'] < 0.85:
                go_result['backup_analyzers'] = ['DN', 'LOF', 'GOF']  # Run all for backup
                go_result['primary_analyzer'] = go_result['analyzers_to_run'][0]  # Primary recommendation
                go_result['analyzers_to_run'] = ['DN', 'LOF', 'GOF']  # Actually run all
                go_result['routing_strategy'] = 'GO_BASED_WITH_BACKUP'
                go_result['rationale'] = f"{go_result['rationale']} (running backup analyzers due to medium confidence)"
            return go_result
        
        # Conservative fallback: run all analyzers
        return {
            'analyzers_to_run': ['DN', 'LOF', 'GOF'],
            'routing_strategy': 'CONSERVATIVE',
            'confidence': 0.5,
            'rationale': f'Insufficient evidence to classify {gene}, running all analyzers for safety',
            'gene_classification': 'UNKNOWN'
        }
    
    def _get_analyzers_for_classification(self, classification: str) -> List[str]:
        """Get list of analyzers to run for a given classification"""
        mapping = {
            'LOF_ONLY': ['LOF'],
            'GOF_ONLY': ['GOF'], 
            'DN_ONLY': ['DN'],
            'MULTI_MECHANISM': ['DN', 'LOF', 'GOF']
        }
        return mapping.get(classification, ['DN', 'LOF', 'GOF'])
    
    def _classify_by_go_terms(self, gene: str, uniprot_id: Optional[str] = None) -> Dict:
        """Classify gene based on GO terms and functional annotations"""
        
        # Get protein annotations
        annotations = self.universal_context.get_context_for_protein(gene, uniprot_id)
        
        if "error" in annotations:
            return {
                'analyzers_to_run': ['DN', 'LOF', 'GOF'],
                'routing_strategy': 'CONSERVATIVE',
                'confidence': 0.3,
                'rationale': f'Could not retrieve annotations for {gene}',
                'gene_classification': 'ANNOTATION_ERROR'
            }
        
        function_text = annotations.get("function", "").lower()
        go_terms = annotations.get("go_terms", [])
        
        # Score each mechanism type
        scores = {'LOF': 0, 'GOF': 0, 'DN': 0}
        
        for mechanism, patterns in self.go_classification_patterns.items():
            mechanism_type = mechanism.split('_')[0]  # LOF, GOF, or DN
            
            for pattern in patterns:
                if pattern in function_text:
                    scores[mechanism_type] += 1
                    
                # Also check GO terms
                for go_term in go_terms:
                    if pattern in go_term.lower():
                        scores[mechanism_type] += 1
        
        # Determine classification based on scores
        max_score = max(scores.values())
        if max_score == 0:
            # No clear indicators
            return {
                'analyzers_to_run': ['DN', 'LOF', 'GOF'],
                'routing_strategy': 'CONSERVATIVE',
                'confidence': 0.4,
                'rationale': f'No clear mechanism indicators for {gene}',
                'gene_classification': 'NO_CLEAR_INDICATORS'
            }
        
        # Find dominant mechanism(s)
        dominant_mechanisms = [mech for mech, score in scores.items() if score == max_score]
        
        if len(dominant_mechanisms) == 1:
            # Clear single mechanism
            mechanism = dominant_mechanisms[0]
            confidence = min(0.85, 0.6 + (max_score * 0.05))  # Cap at 0.85 for GO-based
            
            return {
                'analyzers_to_run': [mechanism],
                'routing_strategy': 'GO_BASED',
                'confidence': confidence,
                'rationale': f'GO terms suggest {mechanism} mechanism for {gene}',
                'gene_classification': f'{mechanism}_BY_GO'
            }
        
        else:
            # Multiple mechanisms or tie
            return {
                'analyzers_to_run': ['DN', 'LOF', 'GOF'],
                'routing_strategy': 'CONSERVATIVE',
                'confidence': 0.6,
                'rationale': f'Multiple mechanism indicators for {gene}',
                'gene_classification': 'MULTI_INDICATORS'
            }


def test_biological_router():
    """Test the biological router with known genes"""
    router = BiologicalRouter()
    
    test_cases = [
        ('CFTR', 'p.F508del', 'missense'),      # Should route to LOF only
        ('FGFR3', 'p.G380R', 'missense'),      # Should route to GOF only  
        ('COL1A1', 'p.G1076S', 'missense'),    # Should route to DN only
        ('TP53', 'p.R273H', 'missense'),       # Should route to multiple
        ('UNKNOWN_GENE', 'p.A123B', 'missense'), # Should be conservative
        ('CFTR', 'p.W1282X', 'nonsense'),      # Should route to LOF (variant type)
    ]
    
    print("🧬 BIOLOGICAL ROUTER TEST RESULTS")
    print("=" * 60)
    
    for gene, variant, var_type in test_cases:
        result = router.route_variant(gene, variant, var_type)
        print(f"\n{gene} {variant} ({var_type}):")
        print(f"  Analyzers: {', '.join(result['analyzers_to_run'])}")
        print(f"  Strategy: {result['routing_strategy']}")
        print(f"  Confidence: {result['confidence']:.2f}")
        print(f"  Rationale: {result['rationale']}")


if __name__ == "__main__":
    test_biological_router()
