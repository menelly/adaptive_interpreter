#!/usr/bin/env python3
"""
🧬 Conservation Nudge Impact Analysis

Compares classification outcomes:
1. NO nudge (pure score-based)
2. Upward-only nudge (conserved → more pathogenic, but no downward push)
3. Current behavior (bidirectional nudge)

Usage: python analyze_conservation_nudge_impact.py <cascade_tsv_file>
"""

import pandas as pd
import sys
from pathlib import Path

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.classifier import VariantClassifier

def classify_score(score: float, classifier: VariantClassifier) -> str:
    """Classify based on score alone"""
    return classifier.interpret_score(score)

def apply_upward_only_nudge(classification: str, phylop: float) -> str:
    """Apply conservation nudge ONLY in upward direction"""
    levels = ['B', 'LB', 'VUS', 'VUS-P', 'LP', 'P']
    if classification not in levels:
        return classification
    
    idx = levels.index(classification)
    
    # Only nudge UP for highly conserved positions
    if phylop is not None and phylop >= 5.0:
        new_idx = min(idx + 1, len(levels) - 1)
        return levels[new_idx]
    
    # NO downward nudge - return unchanged
    return classification

def map_to_bucket(classification: str) -> str:
    """Map classification to P/VUS/B bucket"""
    if classification in ['P', 'LP']:
        return 'P'
    elif classification in ['B', 'LB']:
        return 'B'
    else:
        return 'VUS'

def agreement(clinvar_bucket: str, ai_bucket: str) -> str:
    """Check agreement between ClinVar and AI"""
    if pd.isna(clinvar_bucket) or clinvar_bucket == '':
        return 'NO_CLINVAR'
    if clinvar_bucket == ai_bucket:
        return 'AGREE'
    if clinvar_bucket == 'VUS' or ai_bucket == 'VUS':
        return 'PARTIAL'
    return 'DISAGREE'

def analyze_gene(tsv_path: str):
    """Analyze conservation nudge impact for one gene"""
    
    df = pd.read_csv(tsv_path, sep='\t')
    classifier = VariantClassifier()
    
    gene = df['gene'].iloc[0] if len(df) > 0 else "Unknown"
    print(f"\n{'='*80}")
    print(f"🧬 CONSERVATION NUDGE IMPACT ANALYSIS: {gene}")
    print(f"{'='*80}")
    print(f"Total variants: {len(df)}")
    
    # Get phyloP scores if available (may need to fetch or use conservation_impact as proxy)
    # For now, we'll work backwards from the summary to detect nudges
    
    results = []
    
    for _, row in df.iterrows():
        score = row['final_score']
        current_class = row['final_classification']
        clinvar_bucket = row.get('clinvar_bucket', '')
        summary = row.get('summary', '')
        
        # Pure score-based classification
        no_nudge_class = classify_score(score, classifier)
        
        # Upward-only nudge (we don't have phyloP in file, but we can detect if nudge happened)
        # If current != no_nudge, a nudge was applied
        # For upward-only: if current > no_nudge, keep it; if current < no_nudge, revert
        upward_only_class = no_nudge_class
        
        # Detect nudge direction from summary or by comparing
        levels = ['B', 'LB', 'VUS', 'VUS-P', 'LP', 'P']
        if current_class in levels and no_nudge_class in levels:
            current_idx = levels.index(current_class)
            no_nudge_idx = levels.index(no_nudge_class)
            
            if current_idx > no_nudge_idx:
                # Upward nudge was applied - KEEP IT
                upward_only_class = current_class
            elif current_idx < no_nudge_idx:
                # Downward nudge was applied - REVERT IT
                upward_only_class = no_nudge_class
        
        results.append({
            'variant': row['variant'],
            'score': score,
            'current': current_class,
            'no_nudge': no_nudge_class,
            'upward_only': upward_only_class,
            'clinvar_bucket': clinvar_bucket,
            'current_bucket': map_to_bucket(current_class),
            'no_nudge_bucket': map_to_bucket(no_nudge_class),
            'upward_only_bucket': map_to_bucket(upward_only_class),
        })
    
    results_df = pd.DataFrame(results)
    
    # Calculate agreement stats for each scenario
    for scenario in ['current', 'no_nudge', 'upward_only']:
        bucket_col = f'{scenario}_bucket'
        results_df[f'{scenario}_agreement'] = results_df.apply(
            lambda r: agreement(r['clinvar_bucket'], r[bucket_col]), axis=1
        )
    
    # Print summary stats
    print(f"\n📊 AGREEMENT WITH CLINVAR:")
    print(f"{'Scenario':<20} {'AGREE':<10} {'PARTIAL':<10} {'DISAGREE':<10}")
    print("-" * 50)
    
    for scenario in ['current', 'no_nudge', 'upward_only']:
        agree_col = f'{scenario}_agreement'
        counts = results_df[agree_col].value_counts()
        agree = counts.get('AGREE', 0)
        partial = counts.get('PARTIAL', 0)
        disagree = counts.get('DISAGREE', 0)
        print(f"{scenario:<20} {agree:<10} {partial:<10} {disagree:<10}")
    
    # Show variants that CHANGE between scenarios
    print(f"\n🔄 VARIANTS WHERE REMOVING DOWNWARD NUDGE CHANGES THINGS:")
    changed = results_df[results_df['current'] != results_df['upward_only']]
    print(f"Total changed: {len(changed)}")
    
    if len(changed) > 0:
        print(f"\n{'Variant':<20} {'Score':<8} {'Current':<8} {'UpwardOnly':<12} {'ClinVar':<8} {'Improves?':<10}")
        print("-" * 70)
        for _, r in changed.head(30).iterrows():
            current_agree = r['current_agreement']
            upward_agree = r['upward_only_agreement']
            improves = "✅ YES" if (current_agree == 'DISAGREE' and upward_agree != 'DISAGREE') else ""
            if current_agree == 'PARTIAL' and upward_agree == 'AGREE':
                improves = "✅ YES"
            print(f"{r['variant']:<20} {r['score']:<8.3f} {r['current']:<8} {r['upward_only']:<12} {r['clinvar_bucket']:<8} {improves}")
    
    return results_df

if __name__ == '__main__':
    if len(sys.argv) < 2:
        # Default to MYH7
        tsv_path = '/home/Ace/analysis/acmg_RERUN_POST_THRESHOLD_FIX/MYH7.cascade.tsv'
    else:
        tsv_path = sys.argv[1]
    
    analyze_gene(tsv_path)

