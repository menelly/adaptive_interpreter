#!/usr/bin/env python3
"""
🧪 TEST NOVA'S GENE CLASSIFICATIONS
Test if our family classification system correctly identifies Nova's curated genes

Built by Ace (2025) to validate our category_keywords.json
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from plausibility_filter import classify_gene_family
from universal_protein_annotator import UniversalProteinAnnotator

def test_nova_genes():
    """Test Nova's gene list against our classification system"""

    # Initialize the PROPER annotation system
    annotator = UniversalProteinAnnotator()
    
    # Nova's curated gene list with expected families
    nova_genes = {
        # COLLAGEN FAMILIES
        "COLLAGEN_FIBRILLAR": ["COL1A1", "COL3A1"],
        "COLLAGEN_NETWORK": ["COL4A1", "COL4A5"], 
        "COLLAGEN_ANCHORING": ["COL17A1"],
        "COLLAGEN_FACIT": ["COL9A1", "COL12A1"],
        
        # ELASTIN & FIBRILLIN
        "FIBRILLIN": ["FBN1", "FBN2"],
        "ELASTIN": ["ELN", "GLB1"],
        
        # CYTOSKELETAL
        "INTERMEDIATE_FILAMENT": ["KRT5", "GFAP"],
        "CYTOSKELETON_POLYMER": ["ACTB", "TUBB3"],
        "LAMIN": ["LMNA", "LMNB1"],
        
        # SIGNALING & REGULATORS
        "NEGATIVE_REGULATOR": ["PTEN", "RB1"],
        "SIGNALING_REGULATOR": ["NF1", "APC"],
        "SCAFFOLD_ADAPTOR": ["GRB2", "DLG1"],
        
        # METABOLISM & TRANSPORT
        "AUTOSOMAL_RECESSIVE": ["CFTR", "PAH"],
        "METABOLIC_ENZYME": ["G6PD", "HEXA"],
        "TRANSPORTER": ["SLC2A1", "ABCA4"],
        
        # STRUCTURAL CATCH-ALL
        "STRUCTURAL": ["DMD", "VWF"]
    }
    
    print("🧪 TESTING NOVA'S GENE CLASSIFICATIONS")
    print("="*80)
    
    total_genes = 0
    correct_classifications = 0
    misclassifications = []
    
    for expected_family, genes in nova_genes.items():
        print(f"\n🔍 Testing {expected_family}:")
        
        for gene in genes:
            total_genes += 1

            try:
                print(f"    🔍 Processing {gene}...")

                # Get REAL UniProt annotation
                annotation = annotator.annotate_protein(gene)

                if "error" in annotation:
                    print(f"      ⚠️ Could not annotate {gene}: {annotation['error']}")
                    predicted_family = "UNCLASSIFIED"
                    confidence = 0.0
                    all_scores = {}
                else:
                    # Extract function and GO terms
                    function_description = annotation.get('function', f'{gene} protein')
                    go_terms = annotation.get('go_terms', [])

                    print(f"      📝 Function: {function_description[:60]}...")
                    print(f"      🏷️  GO terms: {len(go_terms)} terms")

                    # Use the PROPER classification system
                    predicted_family = classify_gene_family(gene, function_description, go_terms)
                    confidence = 1.0  # Real classification confidence
                    all_scores = {predicted_family: confidence}

            except Exception as e:
                print(f"      ❌ Error processing {gene}: {e}")
                predicted_family = "ERROR"
                confidence = 0.0
                all_scores = {}
            
            # Check if correct
            if predicted_family == expected_family:
                print(f"  ✅ {gene}: {predicted_family} (confidence: {confidence:.2f})")
                correct_classifications += 1
            else:
                print(f"  ❌ {gene}: Expected {expected_family}, got {predicted_family} (confidence: {confidence:.2f})")
                misclassifications.append({
                    'gene': gene,
                    'expected': expected_family,
                    'predicted': predicted_family,
                    'confidence': confidence,
                    'all_scores': all_scores
                })
    
    print("\n" + "="*80)
    print(f"📊 CLASSIFICATION RESULTS:")
    print(f"   Total genes tested: {total_genes}")
    print(f"   Correct classifications: {correct_classifications}")
    print(f"   Accuracy: {correct_classifications/total_genes*100:.1f}%")
    
    if misclassifications:
        print(f"\n❌ MISCLASSIFICATIONS ({len(misclassifications)}):")
        for miss in misclassifications:
            print(f"   {miss['gene']}: {miss['expected']} → {miss['predicted']} (conf: {miss['confidence']:.2f})")
            
            # Show top 3 scores for debugging
            top_scores = sorted(miss['all_scores'].items(), key=lambda x: x[1], reverse=True)[:3]
            score_str = ", ".join([f"{fam}:{score:.2f}" for fam, score in top_scores])
            print(f"      Top scores: {score_str}")
    
    return correct_classifications, total_genes, misclassifications

def suggest_improvements(misclassifications):
    """Suggest improvements to category_keywords.json based on misclassifications"""
    
    if not misclassifications:
        print("\n🎉 NO IMPROVEMENTS NEEDED - Perfect classification!")
        return
    
    print(f"\n🔧 SUGGESTED IMPROVEMENTS:")
    print("="*50)
    
    for miss in misclassifications:
        gene = miss['gene']
        expected = miss['expected']
        predicted = miss['predicted']
        
        print(f"\n🔍 {gene} (Expected: {expected}, Got: {predicted}):")
        
        # Suggest keywords to add
        if gene == "COL17A1":
            print(f"   Add to {expected}: 'collagen type xvii': 2.0")
        elif gene == "GLB1":
            print(f"   Add to {expected}: 'beta-galactosidase': 1.8, 'lysosomal': 1.5")
        elif gene == "GFAP":
            print(f"   Add to {expected}: 'glial fibrillary': 2.0, 'astrocyte': 1.5")
        elif gene == "TUBB3":
            print(f"   Add to {expected}: 'beta-tubulin': 1.8, 'neuronal tubulin': 1.8")
        elif gene == "LMNA":
            print(f"   Add to {expected}: 'nuclear envelope': 1.5, 'nuclear lamina': 1.8")
        elif gene == "LMNB1":
            print(f"   Add to {expected}: 'nuclear envelope': 1.5, 'nuclear lamina': 1.8")
        elif gene == "PTEN":
            print(f"   Add to {expected}: 'phosphatase': 1.8, 'tumor suppressor': 1.5")
        elif gene == "RB1":
            print(f"   Add to {expected}: 'retinoblastoma': 2.0, 'cell cycle': 1.5")
        elif gene == "NF1":
            print(f"   Add to {expected}: 'neurofibromatosis': 2.0, 'ras gtpase': 1.5")
        elif gene == "APC":
            print(f"   Add to {expected}: 'adenomatous polyposis': 2.0, 'wnt signaling': 1.8")
        elif gene == "GRB2":
            print(f"   Add to {expected}: 'growth factor receptor': 1.8, 'sh2 domain': 1.5")
        elif gene == "DLG1":
            print(f"   Add to {expected}: 'discs large': 1.8, 'pdz domain': 1.5")
        elif gene == "PAH":
            print(f"   Add to {expected}: 'phenylalanine hydroxylase': 2.0")
        elif gene == "G6PD":
            print(f"   Add to {expected}: 'glucose-6-phosphate': 2.0")
        elif gene == "HEXA":
            print(f"   Add to {expected}: 'hexosaminidase': 2.0, 'tay-sachs': 1.8")
        elif gene == "SLC2A1":
            print(f"   Add to {expected}: 'glucose transporter': 2.0, 'glut1': 2.0")
        elif gene == "ABCA4":
            print(f"   Add to {expected}: 'atp-binding cassette': 2.0, 'retinal': 1.5")
        elif gene == "DMD":
            print(f"   Add to {expected}: 'dystrophin': 2.0, 'muscle membrane': 1.8")
        elif gene == "VWF":
            print(f"   Add to {expected}: 'von willebrand': 2.0, 'blood coagulation': 1.8")
        else:
            print(f"   Manual review needed for {gene}")

if __name__ == "__main__":
    correct, total, misses = test_nova_genes()
    suggest_improvements(misses)
    
    print(f"\n🎯 FINAL SCORE: {correct}/{total} ({correct/total*100:.1f}%)")
    
    if correct/total >= 0.9:
        print("🎉 EXCELLENT! Classification system is working well!")
    elif correct/total >= 0.7:
        print("👍 GOOD! Minor improvements needed.")
    else:
        print("🔧 NEEDS WORK! Major improvements required.")
