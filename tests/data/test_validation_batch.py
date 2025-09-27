#!/usr/bin/env python3
"""
🧬 TEST VALIDATION BATCH - Quick test of biological routing batch processing

This creates a small test CSV and runs it through the biological cascade system
to validate our logging and output format before Ren runs the big validation!

Built by Ace (2025) for DNModeling validation
"""

import pandas as pd
from cascade_batch_processor import CascadeBatchProcessor
from pathlib import Path

def create_test_validation_csv():
    """Create a small test CSV with diverse gene types for validation"""
    
    test_variants = [
        # LOF genes (should route to LOF only)
        {'HGVS': 'NM_000492.4(CFTR):c.1521_1523delCTT (p.Phe508del)', 'gnomAD frequency': '0.0001', 'ClinVar_Classification': 'Pathogenic'},
        {'HGVS': 'NM_000492.4(CFTR):c.1652G>A (p.Gly551Asp)', 'gnomAD frequency': '0.00001', 'ClinVar_Classification': 'Pathogenic'},
        {'HGVS': 'NM_001042351.3(G6PD):c.202G>A (p.Val68Met)', 'gnomAD frequency': '0.001', 'ClinVar_Classification': 'Pathogenic'},

        # GOF genes (should route to GOF only)
        {'HGVS': 'NM_000142.5(FGFR3):c.1138G>A (p.Gly380Arg)', 'gnomAD frequency': '0.0', 'ClinVar_Classification': 'Pathogenic'},
        {'HGVS': 'NM_000142.5(FGFR3):c.742C>T (p.Arg248Cys)', 'gnomAD frequency': '0.0', 'ClinVar_Classification': 'Pathogenic'},
        {'HGVS': 'NM_020975.6(RET):c.1900T>C (p.Cys634Arg)', 'gnomAD frequency': '0.0', 'ClinVar_Classification': 'Pathogenic'},

        # DN genes (should route to DN only)
        {'HGVS': 'NM_000088.4(COL1A1):c.3226G>A (p.Gly1076Ser)', 'gnomAD frequency': '0.0', 'ClinVar_Classification': 'Pathogenic'},
        {'HGVS': 'NM_000138.5(FBN1):c.3725G>A (p.Cys1242Tyr)', 'gnomAD frequency': '0.0', 'ClinVar_Classification': 'Pathogenic'},

        # Multi-mechanism genes (should route to multiple analyzers)
        {'HGVS': 'NM_000546.6(TP53):c.817C>T (p.Arg273His)', 'gnomAD frequency': '0.0001', 'ClinVar_Classification': 'Pathogenic'},
        {'HGVS': 'NM_007294.4(BRCA1):c.181T>G (p.Cys61Gly)', 'gnomAD frequency': '0.0', 'ClinVar_Classification': 'Pathogenic'},
        {'HGVS': 'NM_000540.3(RYR1):c.7523G>A (p.Arg2508His)', 'gnomAD frequency': '0.0', 'ClinVar_Classification': 'VUS'},

        # Benign variants (should have lower scores)
        {'HGVS': 'NM_000546.6(TP53):c.215C>G (p.Pro72Arg)', 'gnomAD frequency': '0.25', 'ClinVar_Classification': 'Benign'},
        {'HGVS': 'NM_000492.4(CFTR):c.1408A>G (p.Met470Val)', 'gnomAD frequency': '0.5', 'ClinVar_Classification': 'Benign'},
    ]
    
    df = pd.DataFrame(test_variants)
    test_csv = 'test_validation_input.csv'
    df.to_csv(test_csv, index=False)
    
    print(f"📝 Created test validation CSV: {test_csv}")
    print(f"   Variants: {len(test_variants)}")

    # Count unique genes from HGVS column
    genes = [hgvs.split(':')[0] for hgvs in df['HGVS']]
    print(f"   Genes: {len(set(genes))}")
    
    return test_csv

def run_test_validation():
    """Run the test validation batch"""
    
    print("🧬 BIOLOGICAL ROUTING VALIDATION TEST")
    print("=" * 50)
    
    # Create test CSV
    test_csv = create_test_validation_csv()
    
    # Initialize batch processor
    processor = CascadeBatchProcessor()
    
    # Run batch processing
    output_file = 'test_validation_results.tsv'
    print(f"\n🚀 Running biological cascade batch processing...")
    print(f"   Input: {test_csv}")
    print(f"   Output: {output_file}")
    
    result = processor.process_csv(test_csv, output_file)
    
    if 'error' in result:
        print(f"❌ Batch processing failed: {result['error']}")
        return
    
    # Print results summary
    stats = result['stats']
    print(f"\n📊 VALIDATION TEST RESULTS:")
    print(f"   Total processed: {stats['processed']}")
    print(f"   Successful: {stats['success']}")
    print(f"   Failed: {stats['failed']}")
    
    # Show sample results
    if Path(output_file).exists():
        print(f"\n📋 Sample results from {output_file}:")
        df_results = pd.read_csv(output_file, sep='\t')
        
        # Show routing strategies
        routing_summary = df_results['routing_strategy'].value_counts()
        print(f"\n🧬 Routing Strategy Summary:")
        for strategy, count in routing_summary.items():
            print(f"   {strategy}: {count} variants")
        
        # Show analyzer usage
        print(f"\n⚡ Analyzer Usage:")
        for _, row in df_results.iterrows():
            gene = row['gene']
            analyzers = row['analyzers_run']
            strategy = row['routing_strategy']
            confidence = row['routing_confidence']
            print(f"   {gene}: {analyzers} ({strategy}, {confidence:.2f})")
        
        print(f"\n✅ Full results saved to: {output_file}")
        print(f"   Ready for your big validation run! 🎉")
    
    else:
        print(f"❌ Output file not created: {output_file}")

if __name__ == "__main__":
    run_test_validation()
