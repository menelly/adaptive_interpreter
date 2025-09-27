#!/usr/bin/env python3
"""
🧬 TEST SMART AR DETECTION

Test that our smart AR detection works for CFTR and other genes.
"""

import sys
sys.path.append('.')

from plausibility_filter import classify_gene_family

def test_ar_detection():
    """Test smart AR detection with real function descriptions"""
    
    print("🧬 TESTING SMART AR DETECTION\n")
    
    # Test cases with realistic function descriptions
    test_cases = [
        {
            'gene': 'CFTR',
            'function': 'Epithelial ion channel that plays an important role in the regulation of epithelial ion and water transport. Cystic fibrosis transmembrane conductance regulator.',
            'go_terms': ['chloride channel activity', 'ion transport'],
            'expected': 'AUTOSOMAL_RECESSIVE'
        },
        {
            'gene': 'FBN1', 
            'function': 'Fibrillin-1 is a major component of microfibrils. Marfan syndrome and related connective tissue disorders.',
            'go_terms': ['extracellular matrix', 'structural constituent'],
            'expected': 'FIBRILLIN'
        },
        {
            'gene': 'COL1A1',
            'function': 'Type I collagen alpha-1 chain. Osteogenesis imperfecta and Ehlers-Danlos syndrome.',
            'go_terms': ['collagen', 'extracellular matrix'],
            'expected': 'COLLAGEN_FIBRILLAR'
        },
        {
            'gene': 'SCN1A',
            'function': 'Voltage-gated sodium channel subunit alpha Nav1.1. Epilepsy and seizure disorders.',
            'go_terms': ['sodium channel activity', 'voltage-gated'],
            'expected': 'ION_CHANNEL'
        },
        {
            'gene': 'ABCA4',
            'function': 'ATP-binding cassette transporter. Stargardt disease and retinal dystrophy.',
            'go_terms': ['atp binding', 'transport'],
            'expected': 'AUTOSOMAL_RECESSIVE'
        }
    ]
    
    for case in test_cases:
        result = classify_gene_family(case['gene'], case['function'], case['go_terms'])
        status = "✅" if result == case['expected'] else "❌"
        
        print(f"{status} {case['gene']}: {result} (expected: {case['expected']})")
        if result != case['expected']:
            print(f"   Function: {case['function'][:80]}...")
        print()

if __name__ == "__main__":
    test_ar_detection()
