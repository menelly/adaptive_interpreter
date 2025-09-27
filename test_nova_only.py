#!/usr/bin/env python3
"""
🧬 NOVA-ONLY TESTING SCRIPT
Test specific disagreement variants through NOVA system only
to see if it performs better than CASCADE system

Built by Ace for Ren (2025)
"""

import sys
from nova_dn.mixed_mechanism_resolver import UnifiedMechanismResolver
from cascade_analyzer import CascadeAnalyzer

def test_nova_variant(gene: str, variant: str, expected_clinvar: str):
    """Test a single variant through NOVA system only"""
    
    print(f"\n🧬 TESTING NOVA ONLY: {gene} {variant}")
    print(f"   Expected ClinVar: {expected_clinvar}")
    print("=" * 50)
    
    # Initialize NOVA resolver
    resolver = UnifiedMechanismResolver()
    
    # We need to get sequence data first - use cascade analyzer for that
    cascade = CascadeAnalyzer()
    
    try:
        # Get basic sequence info using UniversalProteinAnnotator
        from universal_protein_annotator import UniversalProteinAnnotator
        annotator = UniversalProteinAnnotator()

        # Get UniProt ID
        uniprot_id = annotator._find_uniprot_id(gene)
        if not uniprot_id:
            print(f"❌ Could not get UniProt ID for {gene}")
            return None

        # Get sequence from UniProt features
        features = annotator.get_uniprot_features(uniprot_id)
        if "error" in features:
            print(f"❌ Could not get features for {gene} ({uniprot_id}): {features['error']}")
            return None

        sequence = features.get('sequence', '')
        if not sequence:
            print(f"❌ Could not get sequence for {gene} ({uniprot_id})")
            return None

        print(f"✅ Got sequence for {gene} ({uniprot_id}): {len(sequence)} residues")
        
        # Parse variant (e.g., p.R811L -> R, 811, L)
        import re
        match = re.match(r'p\.([A-Z])(\d+)([A-Z*])', variant)
        if not match:
            print(f"❌ Could not parse variant: {variant}")
            return None
            
        ref, pos_str, alt = match.groups()
        pos1 = int(pos_str)
        
        print(f"🎯 Parsed variant: {ref}{pos1}{alt}")
        
        # Build context
        context = {
            'gene_name': gene,
            'gene_family': 'GENERAL',  # Default for now
            'domains': [],
            'uniprot_id': uniprot_id,
        }
        
        # Determine variant type
        variant_type = "nonsense" if alt == "*" else "missense"
        
        # 🚀 RUN NOVA ANALYSIS ONLY!
        print(f"🧬 Running NOVA analysis...")
        nova_result = resolver.resolve_mechanisms(
            sequence, pos1, ref, alt, context, variant_type
        )
        
        print(f"🎯 NOVA RESULT:")
        print(f"   Final Score: {nova_result.get('final_score', 0.0):.3f}")
        print(f"   Classification: {nova_result.get('classification', 'Unknown')}")
        print(f"   Rationale: {nova_result.get('rationale', 'No rationale')}")
        
        # Compare to expected
        nova_class = nova_result.get('classification', 'Unknown')
        if expected_clinvar.upper() in ['PATHOGENIC', 'LIKELY_PATHOGENIC']:
            expected_class = 'Pathogenic-range'
        elif expected_clinvar.upper() in ['BENIGN', 'LIKELY_BENIGN']:
            expected_class = 'Benign-range'
        else:
            expected_class = 'VUS-range'
            
        print(f"🎯 COMPARISON:")
        print(f"   NOVA: {nova_class}")
        print(f"   Expected: {expected_class}")
        
        return nova_result
        
    except Exception as e:
        print(f"❌ NOVA analysis failed: {e}")
        return None

def main():
    """Test the disagreement variants"""
    
    # Test cases from Ren's selection
    test_variants = [
        ("MYO7A", "p.R811L", "LIKELY_BENIGN"),
        ("GJB2", "p.V198M", "PATHOGENIC"), 
        ("HNF1A", "p.L107R", "PATHOGENIC"),
        ("HNF1A", "p.R272H", "PATHOGENIC"),
        ("SMAD4", "p.Q249H", "LIKELY_BENIGN"),
        ("NOTCH1", "p.L2314P", "LIKELY_BENIGN"),
    ]
    
    print("🧬 NOVA-ONLY TESTING ON DISAGREEMENT VARIANTS")
    print("Testing variants that CASCADE system disagreed with ClinVar")
    print("=" * 60)
    
    results = []
    for gene, variant, expected in test_variants:
        result = test_nova_variant(gene, variant, expected)
        results.append({
            'gene': gene,
            'variant': variant,
            'expected': expected,
            'nova_result': result
        })
    
    print(f"\n🎯 SUMMARY OF NOVA-ONLY RESULTS:")
    print("=" * 60)
    for r in results:
        if r['nova_result']:
            nova_class = r['nova_result'].get('classification', 'Unknown')
            nova_score = r['nova_result'].get('final_score', 0.0)
            print(f"{r['gene']} {r['variant']}: NOVA={nova_class} ({nova_score:.3f}) vs Expected={r['expected']}")
        else:
            print(f"{r['gene']} {r['variant']}: NOVA=FAILED vs Expected={r['expected']}")

if __name__ == "__main__":
    main()
