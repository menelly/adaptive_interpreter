#!/usr/bin/env python3
"""
🧬 TEST CFTR AR CLASSIFICATION FIX

Test that CFTR now gets AUTOSOMAL_RECESSIVE classification in the full cascade system.
"""

import sys
sys.path.append('.')

from cascade_batch_processor import CascadeBatchProcessor

def test_cftr_ar_fix():
    """Test that CFTR gets AR classification in cascade system"""
    
    processor = CascadeBatchProcessor()
    
    print("🧬 TESTING CFTR AR CLASSIFICATION IN CASCADE SYSTEM\n")
    
    # Test a CFTR variant
    test_data = {
        'Gene': 'CFTR',
        'Variant': 'p.F311L',
        'HGVS': 'NM_000492.4(CFTR):c.931T>C (p.Phe311Leu)',
        'ClinVar': 'LIKELY_PATHOGENIC',
        'Frequency': 0.0
    }
    
    try:
        # Process the variant
        result = processor.process_single_variant(
            test_data['Gene'], 
            test_data['Variant'], 
            test_data['HGVS'],
            test_data['ClinVar'],
            test_data['Frequency']
        )
        
        print(f"Gene: {test_data['Gene']}")
        print(f"Variant: {test_data['Variant']}")
        print(f"Status: {result.get('status', 'UNKNOWN')}")
        
        if 'gene_family' in result:
            print(f"Gene Family: {result['gene_family']}")
            
        if 'plausibility_result' in result:
            print(f"Plausibility Result: {result['plausibility_result']}")
            
        # Check if we see AR classification
        output_text = str(result)
        if 'AUTOSOMAL_RECESSIVE' in output_text:
            print("✅ CFTR correctly classified as AUTOSOMAL_RECESSIVE!")
        elif 'ION_CHANNEL' in output_text:
            print("❌ CFTR still classified as ION_CHANNEL")
        else:
            print("❓ Classification unclear")
            
        print(f"\nFinal Score: {result.get('final_score', 'N/A')}")
        print(f"Classification: {result.get('classification', 'N/A')}")
        
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    test_cftr_ar_fix()
