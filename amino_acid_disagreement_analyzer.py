#!/usr/bin/env python3
"""
Amino Acid Disagreement Pattern Analyzer
Identifies which amino acids are causing the most disagreements with ClinVar
"""

import pandas as pd
import re
from collections import defaultdict, Counter
import sys

def extract_amino_acid_change(variant_string):
    """Extract the amino acid change from variant notation like p.R112Q"""
    
    # Look for pattern like p.R112Q, p.Arg112Gln, etc.
    patterns = [
        r'p\.([A-Z])(\d+)([A-Z])',  # p.R112Q
        r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})',  # p.Arg112Gln
    ]
    
    for pattern in patterns:
        match = re.search(pattern, variant_string)
        if match:
            from_aa = match.group(1)
            position = int(match.group(2))
            to_aa = match.group(3)
            
            # Convert 3-letter to 1-letter if needed
            aa_map = {
                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
            }
            
            if len(from_aa) == 3:
                from_aa = aa_map.get(from_aa, from_aa)
            if len(to_aa) == 3:
                to_aa = aa_map.get(to_aa, to_aa)
                
            return from_aa, position, to_aa
    
    return None, None, None

def analyze_amino_acid_patterns(results_file):
    """Analyze amino acid patterns in disagreements"""
    
    print("🧬 AMINO ACID DISAGREEMENT PATTERN ANALYSIS")
    print("=" * 60)
    
    # Read results
    df = pd.read_csv(results_file, sep='\t')
    
    # Filter to disagreements only
    disagreements = df[df['agreement_flag'] == '❌'].copy()
    
    if len(disagreements) == 0:
        print("🎉 No disagreements found!")
        return
    
    print(f"📊 ANALYZING {len(disagreements)} DISAGREEMENTS")
    print()
    
    # Analyze amino acid patterns
    from_aa_counts = Counter()
    to_aa_counts = Counter()
    change_patterns = Counter()
    
    disagreement_details = []
    
    for _, row in disagreements.iterrows():
        variant = row['variant']
        gene = row['gene']
        our_class = row['final_classification']
        clinvar_class = row['expected_clinvar']
        
        from_aa, position, to_aa = extract_amino_acid_change(variant)
        
        if from_aa and to_aa:
            from_aa_counts[from_aa] += 1
            to_aa_counts[to_aa] += 1
            change_patterns[f"{from_aa}→{to_aa}"] += 1
            
            disagreement_details.append({
                'gene': gene,
                'variant': variant,
                'from_aa': from_aa,
                'to_aa': to_aa,
                'position': position,
                'our_class': our_class,
                'clinvar_class': clinvar_class,
                'change': f"{from_aa}→{to_aa}"
            })
    
    # Print results
    print("🔴 AMINO ACIDS MOST OFTEN CHANGED FROM (causing disagreements):")
    print("-" * 50)
    for aa, count in from_aa_counts.most_common():
        percentage = (count / len(disagreements)) * 100
        print(f"   {aa}: {count} disagreements ({percentage:.1f}%)")
    
    print()
    print("🔵 AMINO ACIDS MOST OFTEN CHANGED TO:")
    print("-" * 50)
    for aa, count in to_aa_counts.most_common():
        percentage = (count / len(disagreements)) * 100
        print(f"   {aa}: {count} disagreements ({percentage:.1f}%)")
    
    print()
    print("🔄 MOST COMMON AMINO ACID CHANGES:")
    print("-" * 50)
    for change, count in change_patterns.most_common(10):
        percentage = (count / len(disagreements)) * 100
        print(f"   {change}: {count} disagreements ({percentage:.1f}%)")
    
    print()
    print("📋 DETAILED BREAKDOWN:")
    print("-" * 50)
    
    # Group by from_aa for detailed analysis
    from_aa_groups = defaultdict(list)
    for detail in disagreement_details:
        from_aa_groups[detail['from_aa']].append(detail)
    
    for from_aa in sorted(from_aa_groups.keys(), key=lambda x: len(from_aa_groups[x]), reverse=True):
        variants = from_aa_groups[from_aa]
        print(f"\n🧬 {from_aa} ({len(variants)} disagreements):")
        
        for variant in variants[:5]:  # Show first 5
            print(f"   {variant['gene']} {variant['variant']}: We={variant['our_class']} vs ClinVar")
        
        if len(variants) > 5:
            print(f"   ... and {len(variants) - 5} more")
    
    # Amino acid properties analysis
    print()
    print("🔬 AMINO ACID PROPERTIES ANALYSIS:")
    print("-" * 50)
    
    aa_properties = {
        'R': 'Positively charged, basic, polar',
        'S': 'Polar, hydroxyl group, phosphorylation site',
        'L': 'Hydrophobic, branched, structural',
        'T': 'Polar, hydroxyl group, phosphorylation site',
        'G': 'Smallest, flexible, structural',
        'P': 'Rigid, proline kink, structural',
        'C': 'Sulfur-containing, disulfide bonds',
        'H': 'Positively charged, aromatic, metal binding',
        'A': 'Small, hydrophobic, structural',
        'V': 'Hydrophobic, branched, structural'
    }
    
    print("Top problematic amino acids and their properties:")
    for aa, count in from_aa_counts.most_common(5):
        properties = aa_properties.get(aa, 'Unknown properties')
        print(f"   {aa} ({count}): {properties}")
    
    return disagreement_details

def main():
    if len(sys.argv) != 2:
        print("Usage: python amino_acid_disagreement_analyzer.py <results_file.tsv>")
        sys.exit(1)
    
    results_file = sys.argv[1]
    analyze_amino_acid_patterns(results_file)

if __name__ == "__main__":
    main()
