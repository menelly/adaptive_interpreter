#!/usr/bin/env python3
"""
🎯 Calculate Optimal Classification Thresholds from ClinVar Ground Truth

Uses ROC/AUC analysis on P/LP vs B/LB variants (no VUS/conflicting)
to find optimal decision boundaries.

Usage: python calculate_optimal_thresholds.py [cascade_dir]
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import sys

def load_all_cascades(cascade_dir: str) -> pd.DataFrame:
    """Load all cascade TSV files from directory"""
    cascade_path = Path(cascade_dir)
    all_dfs = []
    
    for tsv in cascade_path.glob("*.cascade.tsv"):
        try:
            df = pd.read_csv(tsv, sep='\t')
            all_dfs.append(df)
            print(f"  Loaded {tsv.name}: {len(df)} variants")
        except Exception as e:
            print(f"  Failed to load {tsv.name}: {e}")
    
    if all_dfs:
        return pd.concat(all_dfs, ignore_index=True)
    return pd.DataFrame()

def filter_clear_clinvar(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only clear P/LP/B/LB ClinVar calls - no VUS or conflicting"""
    # Map clinical_sig to binary
    pathogenic_terms = ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic']
    benign_terms = ['Benign', 'Likely_benign', 'Benign/Likely_benign']
    
    df = df.copy()
    df['clinvar_binary'] = None
    
    for idx, row in df.iterrows():
        clin = str(row.get('clinical_sig', '')).strip()
        if any(p in clin for p in pathogenic_terms) and 'Conflicting' not in clin:
            df.at[idx, 'clinvar_binary'] = 1  # Pathogenic
        elif any(b in clin for b in benign_terms) and 'Conflicting' not in clin:
            df.at[idx, 'clinvar_binary'] = 0  # Benign
    
    # Filter to only clear calls
    filtered = df[df['clinvar_binary'].notna()].copy()
    filtered['clinvar_binary'] = filtered['clinvar_binary'].astype(int)
    
    return filtered

def find_optimal_thresholds(df: pd.DataFrame):
    """Calculate ROC and find optimal thresholds"""
    
    scores = df['final_score'].values
    labels = df['clinvar_binary'].values
    
    print(f"\n📊 DATASET SUMMARY:")
    print(f"  Total clear variants: {len(df)}")
    print(f"  Pathogenic (P/LP): {sum(labels)}")
    print(f"  Benign (B/LB): {len(labels) - sum(labels)}")
    
    # Calculate ROC curve
    fpr, tpr, thresholds = roc_curve(labels, scores)
    roc_auc = auc(fpr, tpr)
    
    print(f"\n🎯 ROC-AUC: {roc_auc:.4f}")
    
    # Find optimal threshold using Youden's J statistic
    j_scores = tpr - fpr
    optimal_idx = np.argmax(j_scores)
    optimal_threshold = thresholds[optimal_idx]
    
    print(f"\n📍 OPTIMAL SINGLE THRESHOLD (Youden's J):")
    print(f"  Threshold: {optimal_threshold:.3f}")
    print(f"  Sensitivity (TPR): {tpr[optimal_idx]:.3f}")
    print(f"  Specificity (1-FPR): {1-fpr[optimal_idx]:.3f}")
    
    # Find thresholds for 3-class system (B/LB, VUS, LP/P)
    # We want: high confidence benign, uncertain middle, high confidence pathogenic
    
    print(f"\n🎚️  THRESHOLD ANALYSIS FOR 3-CLASS SYSTEM:")
    print(f"  (Looking for LB/VUS boundary and VUS/LP boundary)")
    print()
    
    # Find threshold where we catch 95% of pathogenic (high sensitivity for LP boundary)
    idx_95_sens = np.argmin(np.abs(tpr - 0.95))
    lp_threshold = thresholds[idx_95_sens]
    
    # Find threshold where we have 95% specificity (for LB boundary)  
    idx_95_spec = np.argmin(np.abs((1-fpr) - 0.95))
    lb_threshold = thresholds[idx_95_spec]
    
    print(f"  For 95% sensitivity (catch most P/LP):")
    print(f"    LP threshold ≈ {lp_threshold:.3f} (TPR={tpr[idx_95_sens]:.3f}, FPR={fpr[idx_95_sens]:.3f})")
    
    print(f"\n  For 95% specificity (avoid false P calls on B/LB):")
    print(f"    LB threshold ≈ {lb_threshold:.3f} (TPR={tpr[idx_95_spec]:.3f}, FPR={fpr[idx_95_spec]:.3f})")
    
    # Show current thresholds for comparison
    print(f"\n📋 CURRENT THRESHOLDS (for reference):")
    print(f"  LB/VUS boundary: 0.35")
    print(f"  VUS/LP boundary: 0.70")
    
    # Calculate performance at current thresholds
    print(f"\n📈 PERFORMANCE AT VARIOUS THRESHOLDS:")
    print(f"  {'Threshold':<12} {'Sensitivity':<12} {'Specificity':<12} {'Youden J':<12}")
    print(f"  {'-'*48}")
    
    for thresh in [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80]:
        idx = np.argmin(np.abs(thresholds - thresh))
        sens = tpr[idx]
        spec = 1 - fpr[idx]
        j = sens + spec - 1
        print(f"  {thresh:<12.2f} {sens:<12.3f} {spec:<12.3f} {j:<12.3f}")
    
    return {
        'auc': roc_auc,
        'optimal_threshold': optimal_threshold,
        'lp_threshold_95sens': lp_threshold,
        'lb_threshold_95spec': lb_threshold
    }

if __name__ == '__main__':
    cascade_dir = sys.argv[1] if len(sys.argv) > 1 else '/home/Ace/analysis/acmg_RERUN_POST_THRESHOLD_FIX'
    
    print(f"🧬 Loading cascade files from: {cascade_dir}")
    df = load_all_cascades(cascade_dir)
    
    if len(df) == 0:
        print("No data loaded!")
        sys.exit(1)
    
    print(f"\n📊 Total variants loaded: {len(df)}")
    
    # Filter to clear ClinVar calls
    clear_df = filter_clear_clinvar(df)
    print(f"📊 Clear P/LP/B/LB variants: {len(clear_df)}")
    
    # Calculate thresholds
    results = find_optimal_thresholds(clear_df)

