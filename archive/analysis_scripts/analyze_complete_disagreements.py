#!/usr/bin/env python3
"""
Analyze Complete Disagreements - Find Pathogenic→Benign or Benign→Pathogenic flips
"""

import pandas as pd
import sys
from pathlib import Path

def classify_clinvar(clinvar_text):
    """Classify ClinVar text into simplified categories"""
    if not clinvar_text or pd.isna(clinvar_text):
        return "UNKNOWN"

    text = str(clinvar_text).lower()

    # Pathogenic categories - look for actual pathogenic classifications
    if any(term in text for term in ['pathogenic:', 'likely pathogenic:']):
        return "PATHOGENIC"

    # Benign categories - look for actual benign classifications
    if any(term in text for term in ['benign:', 'likely benign:']):
        return "BENIGN"

    # VUS
    if 'uncertain significance:' in text:
        return "VUS"

    return "UNKNOWN"

def classify_our_result(classification):
    """Classify our results into simplified categories"""
    if not classification or pd.isna(classification):
        return "UNKNOWN"
    
    classification = str(classification).upper()
    
    if classification in ['P', 'LP']:
        return "PATHOGENIC"
    elif classification in ['B', 'LB']:
        return "BENIGN"
    elif classification in ['VUS', 'VUS-P']:
        return "VUS"
    
    return "UNKNOWN"

def analyze_complete_disagreements(results_file):
    """Find complete disagreements - Pathogenic→Benign or Benign→Pathogenic flips"""
    
    print("🔍 ANALYZING COMPLETE DISAGREEMENTS")
    print("=" * 60)
    
    # Read results
    df = pd.read_csv(results_file, sep='\t')
    
    # Classify both our results and ClinVar
    df['clinvar_category'] = df['expected_clinvar'].apply(classify_clinvar)
    df['our_category'] = df['final_classification'].apply(classify_our_result)
    
    # Find complete disagreements
    complete_disagreements = df[
        ((df['clinvar_category'] == 'PATHOGENIC') & (df['our_category'] == 'BENIGN')) |
        ((df['clinvar_category'] == 'BENIGN') & (df['our_category'] == 'PATHOGENIC'))
    ].copy()
    
    print(f"📊 TOTAL VARIANTS ANALYZED: {len(df)}")
    print(f"🔥 COMPLETE DISAGREEMENTS FOUND: {len(complete_disagreements)}")
    print()
    
    if len(complete_disagreements) == 0:
        print("✅ NO COMPLETE DISAGREEMENTS FOUND!")
        print("🎉 This means we have no Pathogenic→Benign or Benign→Pathogenic flips!")
        return
    
    # Analyze disagreement types
    patho_to_benign = complete_disagreements[
        (complete_disagreements['clinvar_category'] == 'PATHOGENIC') & 
        (complete_disagreements['our_category'] == 'BENIGN')
    ]
    
    benign_to_patho = complete_disagreements[
        (complete_disagreements['clinvar_category'] == 'BENIGN') & 
        (complete_disagreements['our_category'] == 'PATHOGENIC')
    ]
    
    print(f"🔴 PATHOGENIC → BENIGN: {len(patho_to_benign)} variants")
    print(f"🟢 BENIGN → PATHOGENIC: {len(benign_to_patho)} variants")
    print()
    
    # Show details for each type
    if len(patho_to_benign) > 0:
        print("🔴 PATHOGENIC → BENIGN DISAGREEMENTS:")
        print("-" * 50)
        for idx, row in patho_to_benign.iterrows():
            print(f"Gene: {row['gene']}")
            print(f"Variant: {row['variant']}")
            print(f"Our Classification: {row['final_classification']}")
            print(f"ClinVar: {row['expected_clinvar']}")
            print(f"Final Score: {row['final_score']:.3f}")
            print(f"Explanation: {row['explanation']}")
            print()
    
    if len(benign_to_patho) > 0:
        print("🟢 BENIGN → PATHOGENIC DISAGREEMENTS:")
        print("-" * 50)
        for idx, row in benign_to_patho.iterrows():
            print(f"Gene: {row['gene']}")
            print(f"Variant: {row['variant']}")
            print(f"Our Classification: {row['final_classification']}")
            print(f"ClinVar: {row['expected_clinvar']}")
            print(f"Final Score: {row['final_score']:.3f}")
            print(f"Explanation: {row['explanation']}")
            print()
    
    # Save complete disagreements to file
    output_file = results_file.replace('.tsv', '_complete_disagreements.tsv')
    complete_disagreements.to_csv(output_file, sep='\t', index=False)
    print(f"💾 Complete disagreements saved to: {output_file}")
    
    return complete_disagreements

def main():
    if len(sys.argv) != 2:
        print("Usage: python analyze_complete_disagreements.py <results_file.tsv>")
        sys.exit(1)
    
    results_file = sys.argv[1]
    
    if not Path(results_file).exists():
        print(f"Error: File {results_file} not found!")
        sys.exit(1)
    
    analyze_complete_disagreements(results_file)

if __name__ == "__main__":
    main()
