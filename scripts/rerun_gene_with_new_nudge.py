#!/usr/bin/env python3
"""
🔄 Rerun a gene through CASCADE with updated code (upward-only nudge)

Uses existing cascade output as input (has protein HGVS already parsed).

Usage: python rerun_gene_with_new_nudge.py PTEN
"""

import sys
import pandas as pd
from pathlib import Path

# Setup paths
sys.path.insert(0, str(Path(__file__).parent.parent))

from analyzers.cascade_analyzer import CascadeAnalyzer

def rerun_gene(gene: str, input_dir: str, output_dir: str):
    """Rerun a gene through CASCADE with new code"""
    
    input_path = Path(input_dir) / f"{gene}.cascade.tsv"
    output_path = Path(output_dir) / f"{gene}.cascade.tsv"
    
    if not input_path.exists():
        print(f"❌ Input file not found: {input_path}")
        return
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Load existing results
    df = pd.read_csv(input_path, sep='\t')
    print(f"🧬 Rerunning {gene}: {len(df)} variants")
    
    # Initialize analyzer
    analyzer = CascadeAnalyzer()
    
    results = []
    
    for idx, row in df.iterrows():
        variant = row['variant']
        gnomad_freq = float(row.get('gnomad_freq', 0) or 0)
        
        if idx % 50 == 0:
            print(f"   Progress: {idx}/{len(df)} ({100*idx/len(df):.0f}%)")
        
        try:
            # Run through cascade with new code
            result = analyzer.analyze_cascade_biological(
                gene=gene,
                variant=variant,
                gnomad_freq=gnomad_freq
            )
            
            # Preserve original ClinVar info
            result['clinical_sig'] = row.get('clinical_sig', '')
            result['review_status'] = row.get('review_status', '')
            result['hgvs'] = row.get('hgvs', '')
            result['variant_type'] = row.get('variant_type', '')
            result['molecular_consequence'] = row.get('molecular_consequence', '')
            
            results.append(result)
            
        except Exception as e:
            print(f"   ❌ Error on {variant}: {e}")
            # Keep original row on error
            results.append(row.to_dict())
    
    # Convert to DataFrame and save
    results_df = pd.DataFrame(results)
    
    # Add agreement columns
    from utils.classifier import VariantClassifier
    classifier = VariantClassifier()
    
    def map_clinvar_bucket(sig):
        sig = str(sig).lower()
        if 'pathogenic' in sig and 'conflicting' not in sig:
            return 'P'
        if 'benign' in sig and 'conflicting' not in sig:
            return 'B'
        if 'conflicting' in sig:
            return 'CONFLICTING'
        return 'VUS'
    
    def map_ai_bucket(classification):
        if classification in ['P', 'LP']:
            return 'P'
        if classification in ['B', 'LB']:
            return 'B'
        return 'VUS'
    
    def agreement_label(cv, ai):
        if cv in ['VUS', 'CONFLICTING', '']:
            return 'N/A'
        if cv == ai:
            return 'AGREE'
        if ai == 'VUS':
            return 'PARTIAL'
        return 'DISAGREE'
    
    results_df['clinvar_bucket'] = results_df['clinical_sig'].apply(map_clinvar_bucket)
    results_df['ai_bucket'] = results_df['final_classification'].apply(map_ai_bucket)
    results_df['agreement'] = results_df.apply(
        lambda r: agreement_label(r['clinvar_bucket'], r['ai_bucket']), axis=1
    )
    
    results_df.to_csv(output_path, sep='\t', index=False)
    print(f"\n✅ Saved to: {output_path}")
    
    # Print summary
    clear = results_df[results_df['clinvar_bucket'].isin(['P', 'B'])]
    print(f"\n📊 AGREEMENT SUMMARY (clear P/B only):")
    print(clear['agreement'].value_counts().to_string())

if __name__ == '__main__':
    gene = sys.argv[1] if len(sys.argv) > 1 else 'PTEN'
    input_dir = '/home/Ace/analysis/acmg_RERUN_POST_THRESHOLD_FIX'
    output_dir = '/home/Ace/analysis/nudge_test'
    
    rerun_gene(gene, input_dir, output_dir)

