#!/usr/bin/env python3
"""
Test Nova's weighted classification system using REAL UniProt/GO data
NO MORE FAKE DATA!!!
"""

import sys
sys.path.append('.')

from plausibility_filter import classify_gene_family
from universal_protein_annotator import UniversalProteinAnnotator
# No need for UniProtMapper - using annotator directly

def test_real_gene_classification():
    """Test gene classification using REAL UniProt and GO data"""
    
    # Test genes from Ren's list
    test_genes = [
        "GJB2",      # Should be ION_CHANNEL (connexin)
        "TGFBR2",    # Should be TUMOR_SUPPRESSOR (TGF-beta)
        "HNF1A",     # Should be TRANSCRIPTION_FACTOR
        "SCN2A",     # Should be ION_CHANNEL (sodium channel)
        "TP53",      # Should be TUMOR_SUPPRESSOR
        "ATP5F1A",   # Should be METABOLIC_ENZYME (DN effects proven!)
        "TFG",       # Should be SCAFFOLD_ADAPTOR (DN effects proven!)
        "G6PD",      # Should be METABOLIC_ENZYME
        "SERPINA1",  # Should be SIGNALING_REGULATOR
        "MYO7A",     # Should be MOTOR_PROTEIN
    ]
    
    print("🧬 TESTING NOVA'S SYSTEM WITH REAL UNIPROT/GO DATA")
    print("=" * 60)
    print("🚫 NO MORE FAKE DATA - USING REAL ANNOTATIONS!")
    print()
    
    # Initialize the real annotation system
    annotator = UniversalProteinAnnotator()

    results = {}

    for gene in test_genes:
        print(f"🔍 Processing {gene}...")

        try:
            # Get real UniProt ID using the same method as cascade analyzer
            uniprot_id = annotator._find_uniprot_id(gene)
            if not uniprot_id:
                print(f"   ❌ Could not find UniProt ID for {gene}")
                continue

            print(f"   📋 UniProt ID: {uniprot_id}")

            # Get REAL UniProt features and GO terms
            features = annotator.get_uniprot_features(uniprot_id)
            
            # Extract function description and GO terms
            function_description = features.get('function', f'{gene} protein')
            go_terms = features.get('go_terms', [])
            
            print(f"   📝 Function: {function_description[:80]}...")
            print(f"   🏷️  GO terms: {len(go_terms)} terms")
            
            # Use Nova's weighted classification system
            classification = classify_gene_family(gene, function_description, go_terms)
            results[gene] = classification
            
            print(f"   🎯 Classification: {classification}")
            print()
            
        except Exception as e:
            print(f"   ❌ Error processing {gene}: {e}")
            print()
            continue
    
    print("=" * 60)
    print("📊 FINAL RESULTS:")
    
    # Count classifications
    classification_counts = {}
    for gene, classification in results.items():
        classification_counts[classification] = classification_counts.get(classification, 0) + 1
        print(f"   {gene:12} → {classification}")
    
    print()
    print("📈 CLASSIFICATION DISTRIBUTION:")
    for family, count in sorted(classification_counts.items()):
        print(f"   {family:25} → {count:2} genes")
    
    print()
    general_count = classification_counts.get("GENERAL", 0)
    total_genes = len(results)
    if total_genes > 0:
        non_general_pct = ((total_genes - general_count) / total_genes) * 100
        print(f"🎯 GENES THAT AVOIDED 'GENERAL': {total_genes - general_count}/{total_genes} ({non_general_pct:.1f}%)")
    
    return results

if __name__ == "__main__":
    test_real_gene_classification()
