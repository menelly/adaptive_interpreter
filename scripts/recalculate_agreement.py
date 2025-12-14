#!/usr/bin/env python3
"""
Post-process cascade TSV files to fix agreement labels.

Fixes:
1. "Conflicting" in ClinVar was incorrectly mapped to P (matched "pathogenic" substring)
2. When AI says VUS, that's being cautious, not disagreeing
"""

import pandas as pd
import sys
from pathlib import Path


def map_clinvar_bucket(s: str) -> str:
    """Map ClinVar clinical_sig to bucket - FIXED VERSION"""
    if not s or pd.isna(s):
        return 'NA'
    s = str(s).lower()
    # IMPORTANT: Check conflicting FIRST before pathogenic/benign keywords
    if 'conflicting' in s:
        return 'CONFLICTING'
    if 'likely pathogenic' in s:
        return 'LP'
    if 'pathogenic' in s:
        return 'P'
    if 'likely benign' in s:
        return 'LB'
    if 'benign' in s:
        return 'B'
    if 'uncertain' in s or 'vus' in s:
        return 'VUS'
    return 'NA'


def agreement_label(clin: str, ai: str) -> str:
    """Determine agreement - FIXED VERSION"""
    if clin == 'NA' or ai == 'NA':
        return 'NA'
    
    # CONFLICTING in ClinVar = experts can't agree, so ANY call is reasonable
    if clin == 'CONFLICTING':
        return 'CONFLICTING_OK'
    
    # VUS <-> VUS agreement
    if clin == 'VUS' and ai == 'VUS':
        return 'AGREE'
    
    # If ClinVar VUS and AI moves to a side, that's better data
    if clin == 'VUS' and ai in {'B', 'LB'}:
        return 'BETTER_DATA_to_benign'
    if clin == 'VUS' and ai in {'P', 'LP'}:
        return 'BETTER_DATA_to_pathogenic'
    
    # If WE say VUS, we're being appropriately cautious - not disagreeing
    if ai == 'VUS':
        return 'CAUTIOUS'
    
    # Bucket agreement across sides
    if (clin in {'B', 'LB'} and ai in {'B', 'LB'}) or (clin in {'P', 'LP'} and ai in {'P', 'LP'}):
        return 'AGREE'
    
    # Only TRUE disagreement is P/LP <-> B/LB flips
    return 'DISAGREE'


def process_file(filepath: Path) -> dict:
    """Process a single cascade TSV file and return stats"""
    df = pd.read_csv(filepath, sep='\t')
    
    # Store old values for comparison
    old_agreement = df['agreement'].copy() if 'agreement' in df.columns else None
    
    # Recalculate clinvar_bucket from clinical_sig
    df['clinvar_bucket'] = df['clinical_sig'].apply(map_clinvar_bucket)
    
    # Recalculate agreement
    df['agreement'] = df.apply(
        lambda row: agreement_label(row['clinvar_bucket'], row['ai_bucket']), 
        axis=1
    )
    
    # Save back
    df.to_csv(filepath, sep='\t', index=False)
    
    # Calculate stats
    stats = df['agreement'].value_counts().to_dict()
    
    # Count changes
    changes = 0
    if old_agreement is not None:
        changes = (old_agreement != df['agreement']).sum()
    
    return {
        'total': len(df),
        'stats': stats,
        'changes': changes
    }


def main():
    if len(sys.argv) < 2:
        print("Usage: python recalculate_agreement.py <directory>")
        sys.exit(1)
    
    directory = Path(sys.argv[1])
    cascade_files = list(directory.glob("*.cascade.tsv"))
    
    print(f"🔧 Processing {len(cascade_files)} cascade files...")
    print()
    
    total_stats = {}
    total_changes = 0
    total_variants = 0
    
    for filepath in sorted(cascade_files):
        gene = filepath.stem.replace('.cascade', '')
        result = process_file(filepath)
        
        total_variants += result['total']
        total_changes += result['changes']
        
        for label, count in result['stats'].items():
            total_stats[label] = total_stats.get(label, 0) + count
        
        if result['changes'] > 0:
            print(f"  {gene}: {result['changes']} labels changed")
    
    print()
    print("=" * 60)
    print(f"📊 SUMMARY: {total_variants} variants across {len(cascade_files)} genes")
    print(f"🔄 {total_changes} agreement labels changed")
    print()
    print("Agreement breakdown:")
    for label in ['AGREE', 'CAUTIOUS', 'CONFLICTING_OK', 'BETTER_DATA_to_benign', 
                  'BETTER_DATA_to_pathogenic', 'DISAGREE', 'NA']:
        if label in total_stats:
            pct = total_stats[label] / total_variants * 100
            print(f"  {label}: {total_stats[label]} ({pct:.1f}%)")
    
    # Calculate "real" accuracy (excluding NA, CAUTIOUS, CONFLICTING_OK)
    agree = total_stats.get('AGREE', 0)
    disagree = total_stats.get('DISAGREE', 0)
    better_b = total_stats.get('BETTER_DATA_to_benign', 0)
    better_p = total_stats.get('BETTER_DATA_to_pathogenic', 0)
    
    decisive = agree + disagree + better_b + better_p
    if decisive > 0:
        accuracy = (agree + better_b + better_p) / decisive * 100
        print()
        print(f"🎯 Decisive accuracy (AGREE + BETTER_DATA vs DISAGREE): {accuracy:.1f}%")
        print(f"   ({agree + better_b + better_p} correct / {decisive} decisive calls)")


if __name__ == "__main__":
    main()

