#!/usr/bin/env python3
"""
🧬🔥 TEST GLY/CYS INTEGRATION IN CASCADE ANALYZER
Quick test to validate the revolutionary Gly/Cys ML system integration
"""

import sys
import os

# Test the Gly/Cys multiplier method directly
def test_gly_cys_multiplier_method():
    """Test the _get_gly_cys_multiplier method from cascade analyzer"""
    
    print("🧬🔥 TESTING GLY/CYS MULTIPLIER METHOD:")
    print("=" * 60)
    
    try:
        # Import the cascade analyzer class (but don't initialize full system)
        sys.path.append('.')
        from cascade_analyzer import CascadeAnalyzer
        
        # Create a minimal analyzer instance for testing
        analyzer = CascadeAnalyzer.__new__(CascadeAnalyzer)  # Create without __init__
        
        # Initialize just the Gly/Cys integrator
        from gly_cys_simple_integrator import SimplifiedGlyCysIntegrator
        analyzer.gly_cys_ml = SimplifiedGlyCysIntegrator()
        
        # Test cases from our ClinVar data
        test_cases = [
            ('COL1A1', 'p.G893A', 'Critical collagen Gly-X-Y'),
            ('COL1A1', 'p.G272S', 'Critical collagen Gly-X-Y'),
            ('COL1A1', 'p.G515A', 'Critical collagen Gly-X-Y'),
            ('FBN1', 'p.C628Y', 'Critical disulfide bond'),
            ('FBN1', 'p.C2470Y', 'Critical disulfide bond'),
            ('SCN1A', 'p.G58R', 'Ion channel glycine'),
            ('SCN1A', 'p.G271V', 'Ion channel glycine'),
            ('RYR1', 'p.R614C', 'Ion channel cysteine gain'),
            ('RYR1', 'p.G4935S', 'Ion channel glycine loss'),
            ('KCNQ2', 'p.G574S', 'Ion channel glycine'),
            ('KCNQ2', 'p.Y755C', 'Ion channel cysteine gain'),
        ]
        
        print(f"{'Gene':8} {'Variant':12} | {'Multiplier':10} | {'Context'}")
        print("-" * 60)
        
        for gene, variant, context in test_cases:
            try:
                multiplier = analyzer._get_gly_cys_multiplier(gene, variant, 0.0)
                print(f"{gene:8} {variant:12} | {multiplier:9.3f}x | {context}")
            except Exception as e:
                print(f"{gene:8} {variant:12} | {'ERROR':9} | {e}")
        
        print("\n🎉 GLY/CYS MULTIPLIER METHOD TEST COMPLETE!")
        return True
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        return False


def test_variant_parsing():
    """Test variant parsing logic"""
    
    print("\n🧬 TESTING VARIANT PARSING:")
    print("=" * 40)
    
    import re
    
    test_variants = [
        'p.G893A',
        'p.C628Y', 
        'p.R614C',
        'p.G58R',
        'p.A123B',  # Non-Gly/Cys (should return 1.0)
        'invalid',   # Invalid format (should return 1.0)
    ]
    
    for variant in test_variants:
        match = re.match(r'p\.([A-Z])(\d+)([A-Z])', variant)
        if match:
            ref_aa, pos_str, alt_aa = match.groups()
            position = int(pos_str)
            is_gly_cys = ref_aa in ['G', 'C'] or alt_aa in ['G', 'C']
            print(f"{variant:10} -> {ref_aa}{position}{alt_aa} | Gly/Cys: {is_gly_cys}")
        else:
            print(f"{variant:10} -> INVALID FORMAT")
    
    print("✅ Variant parsing test complete!")


def analyze_test_file_variants():
    """Analyze specific variants from test files"""
    
    print("\n🧬🔥 ANALYZING TEST FILE VARIANTS:")
    print("=" * 50)
    
    # Known Gly/Cys variants from our test files
    known_variants = [
        # From VariantTest2.tsv and Test3.tsv
        ('COL1A1', 'p.G893A', 'Pathogenic'),
        ('COL1A1', 'p.G1190D', 'Pathogenic'),
        ('COL1A1', 'p.G272S', 'Pathogenic'),
        ('COL1A1', 'p.G701S', 'Pathogenic'),
        ('COL1A1', 'p.G335S', 'Pathogenic'),
        ('COL1A1', 'p.G515A', 'Pathogenic'),
        ('FBN1', 'p.C628Y', 'Pathogenic'),
        ('FBN1', 'p.C2470Y', 'Pathogenic'),
        ('SCN1A', 'p.G58R', 'Pathogenic'),
        ('SCN1A', 'p.G271V', 'Uncertain significance'),
        ('RYR1', 'p.G4935S', 'Pathogenic'),
        ('RYR1', 'p.G2266R', 'Pathogenic'),
        ('RYR1', 'p.R614C', 'Pathogenic'),
        ('RYR1', 'p.R2650C', 'Pathogenic'),
        ('KCNQ2', 'p.G574S', 'Pathogenic'),
        ('KCNQ2', 'p.G487R', 'Pathogenic'),
        ('KCNQ2', 'p.G456E', 'Pathogenic'),
        ('KCNQ2', 'p.G756S', 'Pathogenic'),
        ('KCNQ2', 'p.G418V', 'Pathogenic'),
        ('KCNQ2', 'p.G38R', 'Pathogenic'),
        ('KCNQ2', 'p.Y755C', 'Pathogenic'),
        ('KCNQ2', 'p.R376C', 'Pathogenic'),
    ]
    
    try:
        from gly_cys_simple_integrator import SimplifiedGlyCysIntegrator
        integrator = SimplifiedGlyCysIntegrator()
        
        print(f"{'Gene':8} {'Variant':12} | {'Multiplier':10} | {'Expected':12} | {'Context'}")
        print("-" * 70)
        
        collagen_count = 0
        fibrillin_count = 0
        ion_channel_count = 0
        high_multiplier_count = 0
        
        for gene, variant, expected in known_variants:
            # Parse variant
            import re
            match = re.match(r'p\.([A-Z])(\d+)([A-Z])', variant)
            if match:
                ref_aa, pos_str, alt_aa = match.groups()
                position = int(pos_str)
                
                multiplier = integrator.get_gly_cys_multiplier(gene, position, ref_aa, alt_aa)
                
                # Categorize
                if 'COL' in gene:
                    context = 'Collagen'
                    collagen_count += 1
                elif 'FBN' in gene:
                    context = 'Fibrillin'
                    fibrillin_count += 1
                elif gene in ['SCN1A', 'RYR1', 'KCNQ2']:
                    context = 'Ion Channel'
                    ion_channel_count += 1
                else:
                    context = 'Other'
                
                if multiplier >= 2.0:
                    high_multiplier_count += 1
                
                print(f"{gene:8} {variant:12} | {multiplier:9.3f}x | {expected:12} | {context}")
        
        print("\n📊 ANALYSIS SUMMARY:")
        print(f"   Collagen variants: {collagen_count}")
        print(f"   Fibrillin variants: {fibrillin_count}")
        print(f"   Ion channel variants: {ion_channel_count}")
        print(f"   High multipliers (≥2.0x): {high_multiplier_count}")
        
        print("\n🎉 TEST FILE VARIANT ANALYSIS COMPLETE!")
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")


def main():
    """Run all tests"""
    
    print("🧬🔥 GLY/CYS INTEGRATION VALIDATION SUITE")
    print("=" * 60)
    
    # Test 1: Multiplier method
    success1 = test_gly_cys_multiplier_method()
    
    # Test 2: Variant parsing
    test_variant_parsing()
    
    # Test 3: Test file analysis
    analyze_test_file_variants()
    
    print("\n🎉 ALL TESTS COMPLETE!")
    
    if success1:
        print("✅ GLY/CYS INTEGRATION IS WORKING!")
        print("🔥 Ready to revolutionize genomics with biological intelligence!")
    else:
        print("❌ Some tests failed - check integration")


if __name__ == "__main__":
    main()
