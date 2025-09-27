#!/usr/bin/env python3
from cascade_analyzer import CascadeAnalyzer
import pandas as pd

# Read the test file
df = pd.read_csv('tests/Multivariant - VariantTest2.tsv', sep='\t')
analyzer = CascadeAnalyzer()

print('🧬 HUMAN-FRIENDLY VARIANT ANALYSIS RESULTS')
print('='*80)
print('Gene     Variant      | Mechanism Scores          | Result          | ClinVar')
print('-'*80)

for _, row in df.iterrows():
    gene = row.iloc[0]
    variant = row.iloc[1] 
    clinvar = row.iloc[2] if len(row) > 2 else "Unknown"
    
    try:
        result = analyzer.analyze_cascade(gene, variant, 0.0, variant_type='missense')
        summary = analyzer.format_human_readable(gene, variant, result)
        print(f'{summary} | {clinvar}')
    except Exception as e:
        print(f'{gene:8} {variant:12} | ERROR: {str(e)[:40]}... | {clinvar}')

print('-'*80)
print('🔥 = Hotspot | 🧬 = Conservation | B/LB/VUS/VUS-P/LP/P')
