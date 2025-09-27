#!/usr/bin/env python3
"""
🧬 Human-Friendly Variant Analysis Test
Clean, readable output for clinical review
"""

from cascade_analyzer import CascadeAnalyzer
import sys

def test_variants_clean():
    """Test variants with clean human-readable output"""
    
    analyzer = CascadeAnalyzer()
    
    # Test variants with known ClinVar classifications
    test_cases = [
        # Ion channel variants (conservation issues)
        ("KCNQ2", "p.A265V", 0.0, "Pathogenic"),
        ("SCN1A", "p.L373S", 0.0, "Pathogenic"),
        
        # Collagen variants (hotspot regions)
        ("COL1A1", "p.R361G", 0.0, "Likely benign"),
        ("COL1A1", "p.G359E", 0.0, "VUS"),
        
        # Start codon variants (should be auto-pathogenic)
        ("FKRP", "p.M1L", 0.0, "Pathogenic"),
        ("FKRP", "p.M1V", 0.0, "Pathogenic"),
        
        # Tumor suppressor
        ("TP53", "p.R273H", 0.0, "Pathogenic"),
    ]
    
    print("🧬 HUMAN-FRIENDLY VARIANT ANALYSIS RESULTS")
    print("="*80)
    print(f"{'Gene':8} {'Variant':12} | {'Mechanism Scores':25} | {'Result':15} | ClinVar")
    print("-"*80)
    
    for gene, variant, freq, clinvar in test_cases:
        try:
            result = analyzer.analyze_cascade(gene, variant, freq, variant_type='missense')
            
            # Print one-line summary
            summary = analyzer.format_human_readable(gene, variant, result)
            print(f"{summary} | {clinvar}")
            
        except Exception as e:
            print(f"{gene:8} {variant:12} | ERROR: {str(e)[:40]}...")
    
    print("-"*80)
    print("🔥 = Hotspot detected | 🧬 = Conservation boost applied")
    print("B=Benign, LB=Likely Benign, VUS=Uncertain, VUS-P=VUS-Pathogenic, LP=Likely Pathogenic, P=Pathogenic")

def test_detailed_analysis():
    """Show detailed analysis for a few key variants"""
    
    analyzer = CascadeAnalyzer()
    
    key_variants = [
        ("KCNQ2", "p.A265V", 0.0, "Pathogenic"),
        ("COL1A1", "p.R361G", 0.0, "Likely benign"),
        ("FKRP", "p.M1L", 0.0, "Pathogenic"),
    ]
    
    for gene, variant, freq, clinvar in key_variants:
        try:
            result = analyzer.analyze_cascade(gene, variant, freq, variant_type='missense')
            analyzer.print_human_summary(gene, variant, result, clinvar)
        except Exception as e:
            print(f"\n❌ ERROR analyzing {gene} {variant}: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--detailed":
        test_detailed_analysis()
    else:
        test_variants_clean()
