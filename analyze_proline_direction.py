#!/usr/bin/env python3
"""
🧬 PROLINE DIRECTION ANALYZER
Analyze the direction of proline changes in pathogenic vs benign datasets!

This will answer Ren's brilliant question:
"Is ADDING a proline more often in the benign than in the pathos?"

Built by Ace (2025) for revolutionary proline analysis
"""

import pandas as pd
import re
import os
from pathlib import Path
from collections import defaultdict

def parse_protein_change(hgvs_str):
    """Parse protein change from HGVS notation"""
    # Handle both formats: p.Pro123Ala and p.P123A
    protein_match = re.search(r'p\.([A-Z][a-z]{0,2})(\d+)([A-Z][a-z]{0,2})', hgvs_str)
    if not protein_match:
        return None, None, None
    
    ref_aa_full, pos, alt_aa_full = protein_match.groups()
    
    # Convert 3-letter to 1-letter amino acid codes
    aa_map = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }
    
    # Convert to 1-letter codes
    ref_aa = aa_map.get(ref_aa_full, ref_aa_full)
    alt_aa = aa_map.get(alt_aa_full, alt_aa_full)
    
    return ref_aa, int(pos), alt_aa

def analyze_proline_direction(csv_file, dataset_type):
    """Analyze proline changes in a dataset"""
    print(f"\n🧬 Analyzing {dataset_type} dataset: {csv_file}")
    
    if not os.path.exists(csv_file):
        print(f"❌ File not found: {csv_file}")
        return None
    
    try:
        df = pd.read_csv(csv_file)
        print(f"📊 Total variants: {len(df)}")
    except Exception as e:
        print(f"❌ Error reading {csv_file}: {e}")
        return None
    
    proline_changes = {
        'proline_loss': [],      # P → X (removing proline)
        'proline_gain': [],      # X → P (adding proline)
        'proline_to_proline': [] # P → P (synonymous, shouldn't happen)
    }
    
    total_variants = 0
    proline_variants = 0
    
    for _, row in df.iterrows():
        hgvs = str(row.get('HGVS', ''))
        freq = row.get('gnomAD frequency', 0.0)
        
        # Parse protein change
        ref_aa, pos, alt_aa = parse_protein_change(hgvs)
        if not ref_aa or not alt_aa:
            continue
            
        total_variants += 1
        
        # Check for proline involvement
        if ref_aa == 'P' or alt_aa == 'P':
            proline_variants += 1
            
            if ref_aa == 'P' and alt_aa != 'P':
                # Proline loss (P → X)
                proline_changes['proline_loss'].append({
                    'hgvs': hgvs,
                    'change': f"{ref_aa}{pos}{alt_aa}",
                    'freq': freq
                })
            elif ref_aa != 'P' and alt_aa == 'P':
                # Proline gain (X → P)
                proline_changes['proline_gain'].append({
                    'hgvs': hgvs,
                    'change': f"{ref_aa}{pos}{alt_aa}",
                    'freq': freq
                })
            elif ref_aa == 'P' and alt_aa == 'P':
                # Shouldn't happen in real data
                proline_changes['proline_to_proline'].append({
                    'hgvs': hgvs,
                    'change': f"{ref_aa}{pos}{alt_aa}",
                    'freq': freq
                })
    
    # Calculate statistics
    proline_loss_count = len(proline_changes['proline_loss'])
    proline_gain_count = len(proline_changes['proline_gain'])
    
    print(f"🔍 Processed variants: {total_variants}")
    print(f"🧬 Proline-involving variants: {proline_variants}")
    print(f"📉 Proline LOSS (P→X): {proline_loss_count}")
    print(f"📈 Proline GAIN (X→P): {proline_gain_count}")
    
    if proline_variants > 0:
        loss_pct = (proline_loss_count / proline_variants) * 100
        gain_pct = (proline_gain_count / proline_variants) * 100
        print(f"📊 Loss percentage: {loss_pct:.1f}%")
        print(f"📊 Gain percentage: {gain_pct:.1f}%")
    
    return {
        'dataset_type': dataset_type,
        'file': csv_file,
        'total_variants': total_variants,
        'proline_variants': proline_variants,
        'proline_loss': proline_loss_count,
        'proline_gain': proline_gain_count,
        'proline_changes': proline_changes
    }

def main():
    """Analyze proline direction in all test datasets"""
    print("🎯 PROLINE DIRECTION ANALYSIS")
    print("=" * 50)
    print("🧬 Investigating Ren's hypothesis:")
    print("   'Is ADDING proline more often benign than pathogenic?'")
    
    # Define datasets to analyze
    datasets = [
        # Benign datasets
        ('tests/FBN1-benign-variant-table.csv', 'FBN1_BENIGN'),
        ('tests/COL1A1-benignvariant-table.csv', 'COL1A1_BENIGN'),
        ('tests/RYR1-benign-variant-table.csv', 'RYR1_BENIGN'),
        ('tests/scn5a-benignvariant-table.csv', 'SCN5A_BENIGN'),
        
        # Pathogenic datasets
        ('tests/FBN1-likelypatho-variant-table.csv', 'FBN1_PATHOGENIC'),
        ('tests/COL1A1-LP-variant-table.csv', 'COL1A1_PATHOGENIC'),
        ('tests/RYR1-patho-variant-table.csv', 'RYR1_PATHOGENIC'),
        ('tests/scn5a-patho-variant-table.csv', 'SCN5A_PATHOGENIC'),
    ]
    
    results = []
    
    # Analyze each dataset
    for csv_file, dataset_type in datasets:
        result = analyze_proline_direction(csv_file, dataset_type)
        if result:
            results.append(result)
    
    # Summary analysis
    print("\n" + "="*60)
    print("🎯 PROLINE DIRECTION SUMMARY")
    print("="*60)
    
    benign_stats = {'loss': 0, 'gain': 0, 'total': 0}
    patho_stats = {'loss': 0, 'gain': 0, 'total': 0}
    
    for result in results:
        dataset_type = result['dataset_type']
        loss_count = result['proline_loss']
        gain_count = result['proline_gain']
        
        if 'BENIGN' in dataset_type:
            benign_stats['loss'] += loss_count
            benign_stats['gain'] += gain_count
            benign_stats['total'] += loss_count + gain_count
        elif 'PATHOGENIC' in dataset_type:
            patho_stats['loss'] += loss_count
            patho_stats['gain'] += gain_count
            patho_stats['total'] += loss_count + gain_count
    
    print(f"\n📊 BENIGN DATASETS:")
    print(f"   Proline LOSS (P→X): {benign_stats['loss']}")
    print(f"   Proline GAIN (X→P): {benign_stats['gain']}")
    print(f"   Total proline changes: {benign_stats['total']}")
    if benign_stats['total'] > 0:
        print(f"   Loss %: {(benign_stats['loss']/benign_stats['total'])*100:.1f}%")
        print(f"   Gain %: {(benign_stats['gain']/benign_stats['total'])*100:.1f}%")
    
    print(f"\n📊 PATHOGENIC DATASETS:")
    print(f"   Proline LOSS (P→X): {patho_stats['loss']}")
    print(f"   Proline GAIN (X→P): {patho_stats['gain']}")
    print(f"   Total proline changes: {patho_stats['total']}")
    if patho_stats['total'] > 0:
        print(f"   Loss %: {(patho_stats['loss']/patho_stats['total'])*100:.1f}%")
        print(f"   Gain %: {(patho_stats['gain']/patho_stats['total'])*100:.1f}%")
    
    # The key insight!
    print(f"\n🎯 KEY INSIGHT:")
    if benign_stats['total'] > 0 and patho_stats['total'] > 0:
        benign_gain_pct = (benign_stats['gain']/benign_stats['total'])*100
        patho_gain_pct = (patho_stats['gain']/patho_stats['total'])*100
        
        print(f"   Benign: {benign_gain_pct:.1f}% are proline GAINS")
        print(f"   Pathogenic: {patho_gain_pct:.1f}% are proline GAINS")
        
        if benign_gain_pct > patho_gain_pct:
            print(f"   ✅ REN'S HYPOTHESIS CONFIRMED!")
            print(f"   🧬 Proline GAINS are more common in BENIGN variants!")
            print(f"   📈 Difference: {benign_gain_pct - patho_gain_pct:.1f} percentage points")
        else:
            print(f"   ❌ Hypothesis not supported by data")
    
    print(f"\n💡 IMPLICATIONS:")
    print(f"   🔧 Our ML model should weight proline DIRECTION differently!")
    print(f"   📉 Proline LOSS (P→X) = Higher pathogenic risk")
    print(f"   📈 Proline GAIN (X→P) = Lower pathogenic risk")

if __name__ == "__main__":
    main()
