#!/usr/bin/env python3
"""
Quick cascade analysis on GeneCards TSV export
"""
import pandas as pd
import sys
sys.path.insert(0, '/home/Ace/AdaptiveInterpreter')

from analyzers.cascade_analyzer import analyze_variant_cascade

# Load the TSV
df = pd.read_csv('/home/Ace/analysis/Semi-dominant hypothesis - Sheet1.tsv', sep='\t')

print(f"Loaded {len(df)} variants")
print(f"Columns: {list(df.columns)}")

# Run analysis on each AA change
results = []
for idx, row in df.iterrows():
    aa_chg = row['AA Chg']
    if pd.notna(aa_chg) and aa_chg.strip():
        try:
            result = analyze_variant_cascade(aa_chg)
            result['original_row'] = idx
            result['clinvar_sig'] = row.get('Clinical significance and condition', '')[:50]
            results.append(result)
        except Exception as e:
            print(f"Error on {aa_chg}: {e}")
    
    if idx % 100 == 0:
        print(f"Processed {idx}/{len(df)}...")

print(f"\nAnalyzed {len(results)} variants")

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv('/home/Ace/analysis/cascade_results.tsv', sep='\t', index=False)
print("Saved to /home/Ace/analysis/cascade_results.tsv")

