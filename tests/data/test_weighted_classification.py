#!/usr/bin/env python3
"""
Test the new weighted gene family classification system
"""

import sys
sys.path.append('.')

from plausibility_filter import classify_gene_family

def test_problem_genes():
    """Test our problem genes that were getting misclassified"""
    
    test_cases = [
        {
            "gene": "GJB2",
            "function": "Gap junction beta-2 protein (Connexin-26) (Cx26). Forms connexin channels that allow intercellular communication.",
            "go_terms": ["gap junction", "cell communication", "ion transport"],
            "expected": "ION_CHANNEL"
        },
        {
            "gene": "TGFBR2", 
            "function": "Transforming growth factor beta receptor type 2. Serine/threonine kinase receptor for TGF-beta signaling pathway.",
            "go_terms": ["transforming growth factor", "receptor activity", "signal transduction"],
            "expected": "TUMOR_SUPPRESSOR"
        },
        {
            "gene": "HNF1A",
            "function": "Hepatocyte nuclear factor 1-alpha. Transcription factor with homeodomain that regulates gene expression.",
            "go_terms": ["transcription factor", "dna binding", "sequence-specific dna binding"],
            "expected": "TRANSCRIPTION_FACTOR"
        },
        {
            "gene": "SCN2A",
            "function": "Sodium channel protein type 2 subunit alpha. Voltage-gated sodium channel involved in neuronal excitability.",
            "go_terms": ["sodium channel", "voltage-gated channel", "ion transport"],
            "expected": "ION_CHANNEL"
        },
        {
            "gene": "TP53",
            "function": "Cellular tumor antigen p53. Tumor suppressor protein that acts as guardian of the genome.",
            "go_terms": ["tumor suppressor", "dna damage response", "transcription factor"],
            "expected": "TUMOR_SUPPRESSOR"  # Should be primary despite TF activity
        }
    ]
    
    print("🧬 TESTING NOVA'S WEIGHTED GENE FAMILY CLASSIFICATION SYSTEM")
    print("=" * 70)
    
    for case in test_cases:
        gene = case["gene"]
        function = case["function"]
        go_terms = case["go_terms"]
        expected = case["expected"]
        
        result = classify_gene_family(gene, function, go_terms)
        
        status = "✅ CORRECT" if result == expected else "❌ INCORRECT"
        
        print(f"\n🔍 {gene}:")
        print(f"   Function: {function[:80]}...")
        print(f"   Expected: {expected}")
        print(f"   Got:      {result}")
        print(f"   Status:   {status}")
        
        if result != expected:
            print(f"   🚨 MISCLASSIFICATION DETECTED!")

if __name__ == "__main__":
    test_problem_genes()
