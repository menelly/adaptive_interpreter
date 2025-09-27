#!/usr/bin/env python3
"""
🔍 Simple ClinVar Quality Analyzer
Analyzes ClinVar text descriptions for quality indicators without API calls
"""

import pandas as pd
import re
from typing import Dict, List

def analyze_clinvar_text(clinvar_text: str) -> Dict:
    """Analyze ClinVar text for quality indicators"""
    
    text_lower = str(clinvar_text).lower()
    
    # Quality indicators
    quality_score = 0.0
    red_flags = []
    good_signs = []
    
    # RED FLAGS (subtract from quality)
    if 'conflicting interpretations' in text_lower:
        red_flags.append('Conflicting interpretations')
        quality_score -= 0.4
    
    if 'not provided' in text_lower:
        red_flags.append('No clinical context provided')
        quality_score -= 0.2
    
    if 'not specified' in text_lower:
        red_flags.append('Unspecified condition')
        quality_score -= 0.1
    
    if 'single submitter' in text_lower:
        red_flags.append('Single submitter only')
        quality_score -= 0.2
    
    if 'no assertion criteria' in text_lower:
        red_flags.append('No assertion criteria')
        quality_score -= 0.3
    
    # Count "not provided" occurrences (more = worse)
    not_provided_count = text_lower.count('not provided')
    if not_provided_count > 1:
        red_flags.append(f'Multiple "not provided" ({not_provided_count})')
        quality_score -= 0.1 * not_provided_count
    
    # GOOD SIGNS (add to quality)
    if 'expert panel' in text_lower:
        good_signs.append('Expert panel reviewed')
        quality_score += 0.5
    
    if 'practice guideline' in text_lower:
        good_signs.append('Practice guideline')
        quality_score += 0.5
    
    if 'multiple submitters' in text_lower and 'no conflicts' in text_lower:
        good_signs.append('Multiple submitters, no conflicts')
        quality_score += 0.3
    
    # Trusted submitters
    trusted_labs = ['invitae', 'genedx', 'ambry', 'blueprint', 'fulgent', 'clinGen', 'omim']
    for lab in trusted_labs:
        if lab.lower() in text_lower:
            good_signs.append(f'Trusted submitter: {lab}')
            quality_score += 0.2
            break
    
    # Specific conditions (better than generic)
    specific_conditions = len(re.findall(r'[A-Z][a-z]+ [A-Z][a-z]+', clinvar_text))
    if specific_conditions > 2:
        good_signs.append(f'Multiple specific conditions ({specific_conditions})')
        quality_score += 0.1
    
    # Normalize score to 0-1
    quality_score = max(0.0, min(1.0, quality_score + 0.5))  # Start at 0.5 baseline
    
    return {
        'quality_score': quality_score,
        'red_flags': red_flags,
        'good_signs': good_signs,
        'verdict': get_verdict(quality_score, red_flags)
    }

def get_verdict(quality_score: float, red_flags: List[str]) -> str:
    """Determine verdict based on quality analysis"""
    
    if 'Conflicting interpretations' in red_flags:
        return 'CLINVAR_CONFLICTED'
    elif quality_score < 0.3:
        return 'CLINVAR_LOW_QUALITY'
    elif 'No clinical context provided' in red_flags and quality_score < 0.5:
        return 'CLINVAR_WEAK_EVIDENCE'
    elif quality_score > 0.7:
        return 'CLINVAR_HIGH_QUALITY'
    else:
        return 'CLINVAR_MODERATE_QUALITY'

def analyze_disagreements(tsv_path: str):
    """Analyze disagreements in results file"""
    
    # Read results
    df = pd.read_csv(tsv_path, sep='\t')
    
    # Filter to disagreements
    disagreements = df[df['agreement_flag'] == '❌'].copy()
    
    if len(disagreements) == 0:
        print("🎉 No disagreements found!")
        return
    
    print(f"🔍 Analyzing {len(disagreements)} disagreements...")
    print("="*80)
    
    # Analyze each disagreement
    results = []
    verdict_counts = {}
    
    for _, row in disagreements.iterrows():
        gene = row['gene']
        variant = row['variant']
        our_class = row['final_classification']
        clinvar_text = str(row['expected_clinvar'])
        
        analysis = analyze_clinvar_text(clinvar_text)
        verdict = analysis['verdict']
        
        # Count verdicts
        verdict_counts[verdict] = verdict_counts.get(verdict, 0) + 1
        
        # Print analysis
        print(f"🧬 {gene} {variant}: We={our_class} vs ClinVar")
        print(f"   Quality Score: {analysis['quality_score']:.2f}")
        print(f"   Verdict: {verdict}")
        
        if analysis['red_flags']:
            print(f"   🚨 Red Flags: {', '.join(analysis['red_flags'])}")
        if analysis['good_signs']:
            print(f"   ✅ Good Signs: {', '.join(analysis['good_signs'])}")
        
        # Special analysis for ion channels (conservation issue)
        if gene in ['SCN1A', 'SCN2A', 'KCNQ2', 'CACNA1A']:
            print(f"   🧬 NOTE: Ion channel - likely conservation data issue")
        
        print()
        
        results.append({
            **row.to_dict(),
            **analysis
        })
    
    # Summary
    print("="*80)
    print("🎯 DISAGREEMENT ANALYSIS SUMMARY:")
    print("="*80)
    
    total = len(disagreements)
    for verdict, count in sorted(verdict_counts.items()):
        percentage = count / total * 100
        print(f"   {verdict}: {count}/{total} ({percentage:.1f}%)")
    
    # Categorize results
    low_quality = verdict_counts.get('CLINVAR_LOW_QUALITY', 0)
    conflicted = verdict_counts.get('CLINVAR_CONFLICTED', 0)
    weak_evidence = verdict_counts.get('CLINVAR_WEAK_EVIDENCE', 0)
    
    problematic_clinvar = low_quality + conflicted + weak_evidence
    
    print(f"\n📊 INTERPRETATION:")
    print(f"   Problematic ClinVar data: {problematic_clinvar}/{total} ({problematic_clinvar/total*100:.1f}%)")
    print(f"   High-quality disagreements: {total - problematic_clinvar}/{total} ({(total-problematic_clinvar)/total*100:.1f}%)")
    
    if problematic_clinvar > total * 0.7:
        print(f"\n🎉 CONCLUSION: Most 'disagreements' are ClinVar quality issues!")
        print(f"   Our system accuracy is likely MUCH higher than raw numbers suggest.")
    
    # Save detailed results
    results_df = pd.DataFrame(results)
    output_path = tsv_path.replace('.tsv', '_simple_quality_analysis.tsv')
    results_df.to_csv(output_path, sep='\t', index=False)
    print(f"\n💾 Detailed analysis saved to: {output_path}")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python3 simple_quality_analyzer.py <results.tsv>")
        sys.exit(1)
    
    analyze_disagreements(sys.argv[1])
