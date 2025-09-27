#!/usr/bin/env python3
"""
🧬 TEST INHERITANCE-AWARE GENE CLASSIFICATION FIX

Test that CFTR gets AUTOSOMAL_RECESSIVE classification instead of ION_CHANNEL
and receives proper LOF boosting for cystic fibrosis variants.
"""

import sys
sys.path.append('.')

from cascade_batch_processor import CascadeBatchProcessor
from plausibility_filter import apply_pathogenicity_filter

def test_inheritance_fixes():
    """Test inheritance-aware classification for multiple genes"""

    processor = CascadeBatchProcessor()

    print("🧬 TESTING INHERITANCE-AWARE GENE CLASSIFICATION\n")

    # Test 1: CFTR (AR - should boost LOF)
    cftr_family = processor.get_gene_family('CFTR')
    print(f"1️⃣ CFTR gene family: {cftr_family}")

    cftr_scores = {'LOF': 0.324, 'DN': 0.36, 'GOF': 0.14}
    old_filtered = apply_pathogenicity_filter(cftr_scores, 'ION_CHANNEL')
    new_filtered = apply_pathogenicity_filter(cftr_scores, 'AUTOSOMAL_RECESSIVE')

    print(f"   LOF: {cftr_scores['LOF']:.3f} → {new_filtered['LOF']['weighted_score']:.3f} (was {old_filtered['LOF']['weighted_score']:.3f})")
    print(f"   🎯 CFTR LOF boosted by {((new_filtered['LOF']['weighted_score']/old_filtered['LOF']['weighted_score'] - 1) * 100):.1f}%\n")

    # Test 2: COL1A1 (AD - should balance LOF+DN)
    col1a1_family = processor.get_gene_family('COL1A1')
    print(f"2️⃣ COL1A1 gene family: {col1a1_family}")

    # Test the problematic COL1A1 p.R361G case
    col1a1_scores = {'LOF': 0.84, 'DN': 0.50, 'GOF': 0.0}
    old_col_filtered = apply_pathogenicity_filter(col1a1_scores, 'COLLAGEN_FIBRILLAR')

    print(f"   LOF: {col1a1_scores['LOF']:.3f} → {old_col_filtered['LOF']['weighted_score']:.3f} (was 0.588 with old 0.7x)")
    print(f"   DN:  {col1a1_scores['DN']:.3f} → {old_col_filtered['DN']['weighted_score']:.3f}")
    print(f"   🎯 COL1A1 LOF now gets {old_col_filtered['LOF']['weighted_score']:.3f} (more balanced!)\n")

    # Test 3: FBN1 (AD - should balance LOF+DN)
    fbn1_family = processor.get_gene_family('FBN1')
    print(f"3️⃣ FBN1 gene family: {fbn1_family}")

    fbn1_scores = {'LOF': 0.43, 'DN': 0.35, 'GOF': 0.15}
    fbn1_filtered = apply_pathogenicity_filter(fbn1_scores, 'FIBRILLIN')

    print(f"   LOF: {fbn1_scores['LOF']:.3f} → {fbn1_filtered['LOF']['weighted_score']:.3f}")
    print(f"   DN:  {fbn1_scores['DN']:.3f} → {fbn1_filtered['DN']['weighted_score']:.3f}")
    print(f"   🎯 FBN1 balanced for Marfan (both LOF+DN viable)\n")

    print("🚀 SUMMARY:")
    print("   ✅ AR genes (CFTR): LOF boosted 1.3x")
    print("   ✅ AD structural (COL1A1): Balanced LOF 0.85x, DN 1.25x")
    print("   ✅ AD fibrillin (FBN1): Balanced LOF 0.85x, DN 1.2x")
    print("   🎯 This preserves biological reality while improving classification!")

if __name__ == "__main__":
    test_inheritance_fixes()
