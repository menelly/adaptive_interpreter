#!/usr/bin/env python3
"""
🧬 TEST NOVA VARIANT PARSING FIX

Test that Nova can now parse variants with "p." prefix correctly.
"""

import sys
sys.path.append('.')

from cascade_batch_processor import CascadeBatchProcessor

def test_nova_parsing():
    """Test Nova variant parsing with different formats"""
    
    processor = CascadeBatchProcessor()
    
    # Mock cascade result for testing
    mock_cascade_result = {
        'status': 'SUCCESS',
        'gene': 'FBN1',
        'variant': 'p.D2411N',
        'uniprot_id': 'P35555',
        'domains': [],
        'final_score': 0.5
    }
    
    test_variants = [
        'p.D2411N',  # With p. prefix (problematic)
        'D2411N',    # Without prefix (should work)
        'p.P840L',   # Another with prefix
        'P840L'      # Another without prefix
    ]
    
    print("🧬 TESTING NOVA VARIANT PARSING\n")
    
    for variant in test_variants:
        print(f"Testing variant: {variant}")
        try:
            result = processor.run_nova_analysis('FBN1', variant, mock_cascade_result)
            print(f"   ✅ SUCCESS: {result}")
        except Exception as e:
            print(f"   ❌ FAILED: {e}")
        print()

if __name__ == "__main__":
    test_nova_parsing()
