#!/usr/bin/env python3
"""
🎯 Simulate New Threshold Impact on Agreement

Compares current vs proposed thresholds across all cascade files.
"""

import pandas as pd
from pathlib import Path
import sys

# Current thresholds
CURRENT = {'lb_vus': 0.35, 'vus_lp': 0.70}
# Proposed thresholds (data-driven)
PROPOSED = {'lb_vus': 0.32, 'vus_lp': 0.65}

def classify_with_thresholds(score: float, lb_vus: float, vus_lp: float) -> str:
    """Classify score using given thresholds"""
    if score >= vus_lp + 0.10:  # P threshold (LP + 0.10)
        return 'P'
    elif score >= vus_lp:
        return 'LP'
    elif score >= lb_vus + 0.15:  # VUS-P threshold
        return 'VUS-P'
    elif score >= lb_vus:
        return 'VUS'
    elif score >= lb_vus - 0.15:  # LB threshold
        return 'LB'
    else:
        return 'B'

def map_to_bucket(classification: str) -> str:
    """Map to P/VUS/B bucket"""
    if classification in ['P', 'LP']:
        return 'P'
    elif classification in ['B', 'LB']:
        return 'B'
    return 'VUS'

def agreement(clinvar_bucket: str, ai_bucket: str) -> str:
    """Check agreement"""
    if pd.isna(clinvar_bucket) or clinvar_bucket == '' or clinvar_bucket == 'nan':
        return 'NO_CLINVAR'
    if clinvar_bucket == 'VUS' or clinvar_bucket == 'CONFLICTING':
        return 'CLINVAR_VUS'
    if clinvar_bucket == ai_bucket:
        return 'AGREE'
    if ai_bucket == 'VUS':
        return 'PARTIAL'
    return 'DISAGREE'

def analyze_thresholds(cascade_dir: str):
    """Compare current vs proposed thresholds"""
    
    cascade_path = Path(cascade_dir)
    all_results = []
    
    print(f"🧬 Loading cascade files from: {cascade_dir}\n")
    
    for tsv in sorted(cascade_path.glob("*.cascade.tsv")):
        try:
            df = pd.read_csv(tsv, sep='\t')
            
            for _, row in df.iterrows():
                score = row['final_score']
                clinvar_bucket = str(row.get('clinvar_bucket', ''))
                
                # Current classification
                current_class = classify_with_thresholds(score, CURRENT['lb_vus'], CURRENT['vus_lp'])
                current_bucket = map_to_bucket(current_class)
                
                # Proposed classification
                proposed_class = classify_with_thresholds(score, PROPOSED['lb_vus'], PROPOSED['vus_lp'])
                proposed_bucket = map_to_bucket(proposed_class)
                
                all_results.append({
                    'gene': row['gene'],
                    'variant': row['variant'],
                    'score': score,
                    'clinvar_bucket': clinvar_bucket,
                    'current_class': current_class,
                    'current_bucket': current_bucket,
                    'proposed_class': proposed_class,
                    'proposed_bucket': proposed_bucket,
                    'current_agreement': agreement(clinvar_bucket, current_bucket),
                    'proposed_agreement': agreement(clinvar_bucket, proposed_bucket),
                })
        except Exception as e:
            print(f"  Failed {tsv.name}: {e}")
    
    results_df = pd.DataFrame(all_results)
    
    # Filter to variants with clear ClinVar (P or B buckets only)
    clear_df = results_df[results_df['clinvar_bucket'].isin(['P', 'B'])]
    
    print(f"📊 TOTAL VARIANTS: {len(results_df)}")
    print(f"📊 CLEAR CLINVAR (P/B only): {len(clear_df)}")
    
    print(f"\n{'='*60}")
    print(f"📈 AGREEMENT COMPARISON (Clear ClinVar P/B only)")
    print(f"{'='*60}")
    print(f"\n{'Metric':<20} {'Current':<15} {'Proposed':<15} {'Change':<10}")
    print(f"{'-'*60}")
    
    for metric in ['AGREE', 'PARTIAL', 'DISAGREE']:
        current_count = len(clear_df[clear_df['current_agreement'] == metric])
        proposed_count = len(clear_df[clear_df['proposed_agreement'] == metric])
        change = proposed_count - current_count
        change_str = f"+{change}" if change > 0 else str(change)
        emoji = "✅" if (metric == 'AGREE' and change > 0) or (metric == 'DISAGREE' and change < 0) else ""
        print(f"{metric:<20} {current_count:<15} {proposed_count:<15} {change_str:<10} {emoji}")
    
    # Show variants that flip from DISAGREE to better
    print(f"\n🔄 DISAGREE → BETTER (variants that improve):")
    improved = clear_df[
        (clear_df['current_agreement'] == 'DISAGREE') & 
        (clear_df['proposed_agreement'] != 'DISAGREE')
    ]
    print(f"   Count: {len(improved)}")
    
    if len(improved) > 0:
        print(f"\n   {'Gene':<10} {'Variant':<18} {'Score':<8} {'Current':<8} {'Proposed':<10} {'ClinVar':<8}")
        print(f"   {'-'*70}")
        for _, r in improved.head(20).iterrows():
            print(f"   {r['gene']:<10} {r['variant']:<18} {r['score']:<8.3f} {r['current_class']:<8} {r['proposed_class']:<10} {r['clinvar_bucket']:<8}")
    
    # Show variants that flip from better to DISAGREE (if any)
    print(f"\n⚠️  BETTER → DISAGREE (variants that get worse):")
    worsened = clear_df[
        (clear_df['current_agreement'] != 'DISAGREE') & 
        (clear_df['proposed_agreement'] == 'DISAGREE')
    ]
    print(f"   Count: {len(worsened)}")
    
    if len(worsened) > 0:
        print(f"\n   {'Gene':<10} {'Variant':<18} {'Score':<8} {'Current':<8} {'Proposed':<10} {'ClinVar':<8}")
        print(f"   {'-'*70}")
        for _, r in worsened.head(20).iterrows():
            print(f"   {r['gene']:<10} {r['variant']:<18} {r['score']:<8.3f} {r['current_class']:<8} {r['proposed_class']:<10} {r['clinvar_bucket']:<8}")

if __name__ == '__main__':
    cascade_dir = sys.argv[1] if len(sys.argv) > 1 else '/home/Ace/analysis/acmg_RERUN_POST_THRESHOLD_FIX'
    analyze_thresholds(cascade_dir)

