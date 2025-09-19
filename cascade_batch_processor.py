#!/usr/bin/env python3
"""
🌊 CASCADE BATCH PROCESSOR - Multi-Analyzer CSV Processing
Process CSV files through the complete DN → LOF → GOF cascade system!

Built by Ace & Nova (2025) for revolutionary batch genetics analysis

Features:
- Processes ClinVar exports and custom CSV files
- Full cascade logic: DN → (LOF + GOF if needed)
- Population frequency filtering
- Comprehensive error handling
- TSV output with all analyzer scores

Usage:
  python3 cascade_batch_processor.py --input variants.csv --output results.tsv
"""

import argparse
import csv
import json
import os
import sys
from typing import Dict, List
from pathlib import Path
import re

from cascade_analyzer import CascadeAnalyzer
from nova_dn.csv_batch_processor import CSVBatchProcessor


class CascadeBatchProcessor:
    """Process CSV files through the complete cascade system with ClinVar comparison"""

    def __init__(self, alphafold_path: str = "/mnt/Arcana/alphafold_human/structures/"):
        self.cascade_analyzer = CascadeAnalyzer(alphafold_path)
        self.csv_processor = CSVBatchProcessor(alphafold_path)  # For HGVS parsing

    def detect_clinvar_expectation(self, filename: str) -> str:
        """Detect expected ClinVar classification from filename"""
        filename_lower = filename.lower()

        if 'benign' in filename_lower or '-lb-' in filename_lower or '-b-' in filename_lower:
            return 'BENIGN'
        elif 'patho' in filename_lower or '-lp-' in filename_lower or '-p-' in filename_lower:
            return 'PATHOGENIC'
        elif 'vus' in filename_lower:
            return 'VUS'
        else:
            return 'UNKNOWN'

    def determine_agreement(self, our_classification: str, expected_clinvar: str) -> tuple:
        """
        Determine agreement status between our prediction and ClinVar expectation

        Returns: (agreement_status, agreement_flag)
        """
        our_class = our_classification.upper()
        expected = expected_clinvar.upper()

        # ✅ AGREEMENT cases
        if (our_class == 'LB' and expected == 'BENIGN') or \
           (our_class in ['VUS', 'VUS-P'] and expected == 'VUS') or \
           (our_class == 'LP' and expected == 'PATHOGENIC'):
            return ('AGREE', '✅')

        # 🎯 BETTER DATA cases (not wrong, we have more confidence!)
        elif (our_class == 'LB' and expected == 'VUS'):
            return ('BETTER_DATA_BENIGN', '🎯')
        elif (our_class in ['LP', 'VUS-P'] and expected == 'VUS'):
            return ('BETTER_DATA_PATHOGENIC', '🎯')

        # ❌ DISAGREEMENT cases (real conflicts)
        elif (our_class == 'LB' and expected == 'PATHOGENIC') or \
             (our_class == 'LP' and expected == 'BENIGN'):
            return ('DISAGREE', '❌')

        # Unknown/unclear cases
        else:
            return ('UNCLEAR', '❓')
        
    def process_csv(self, input_path: str, output_path: str, 
                   freq_threshold: float = 0.01) -> Dict:
        """
        Process CSV file through cascade analysis
        
        Args:
            input_path: Path to input CSV file
            output_path: Path to output TSV file  
            freq_threshold: Skip variants with frequency > threshold
        """
        
        # 🎯 Detect expected ClinVar classification from filename
        expected_clinvar = self.detect_clinvar_expectation(input_path)

        results = []
        stats = {
            'total_variants': 0,
            'processed': 0,
            'skipped_frequency': 0,
            'skipped_frameshift': 0,
            'skipped_synonymous': 0,
            'skipped_intronic': 0,
            'skipped_unparseable': 0,
            'cascade_triggered': 0,
            'dn_only': 0,
            'failed': 0,
            'success': 0,
            # 🎯 ClinVar comparison stats
            'agreements': 0,
            'better_data_benign': 0,
            'better_data_pathogenic': 0,
            'disagreements': 0,
            'unclear': 0
        }

        print(f"🌊 Processing variants through CASCADE SYSTEM")
        print(f"📊 Frequency threshold: {freq_threshold}")
        print(f"🎯 Expected ClinVar classification: {expected_clinvar}")
        print("=" * 60)
        
        try:
            with open(input_path, 'r') as f:
                reader = csv.DictReader(f)
                
                for row in reader:
                    stats['total_variants'] += 1
                    
                    # Parse HGVS using existing CSV processor logic
                    hgvs = row.get('HGVS', '')
                    freq_str = row.get('gnomAD frequency', '0')
                    
                    # Parse frequency
                    try:
                        gnomad_freq = float(freq_str) if freq_str else 0.0
                    except ValueError:
                        gnomad_freq = 0.0
                    
                    # Apply frequency filter
                    if gnomad_freq > freq_threshold:
                        stats['skipped_frequency'] += 1
                        continue
                    
                    # Parse variant using existing logic
                    parsed = self.csv_processor.parse_hgvs(hgvs)
                    if not parsed:
                        stats['skipped_unparseable'] += 1
                        continue

                    # Skip non-missense variants
                    if parsed['type'] != 'missense':
                        skip_reason = parsed.get('skip_reason', 'SKIPPED_OTHER')
                        if skip_reason == 'SKIPPED_FRAMESHIFT':
                            stats['skipped_frameshift'] += 1
                            print(f"⏭️  Skipping {parsed['gene']} {parsed['variant']} (frameshift/nonsense)")
                        elif skip_reason == 'SKIPPED_SYNONYMOUS':
                            stats['skipped_synonymous'] += 1
                            print(f"⏭️  Skipping {parsed['gene']} {parsed['variant']} (synonymous)")
                        elif skip_reason == 'SKIPPED_INTRONIC':
                            stats['skipped_intronic'] += 1
                            print(f"⏭️  Skipping {parsed['gene']} intronic variant")
                        else:
                            stats['skipped_unparseable'] += 1
                            print(f"⏭️  Skipping unparseable variant: {hgvs}")
                        continue
                    
                    gene = parsed['gene']
                    variant = parsed['variant']
                    
                    print(f"\n🧬 Processing {gene} {variant} (freq: {gnomad_freq:.4f})")
                    
                    # 🚀 Run BIOLOGICAL CASCADE analysis (smarter routing!)
                    result = self.cascade_analyzer.analyze_cascade_biological(
                        gene, variant, gnomad_freq, 'missense'  # Default to missense
                    )
                    
                    # Add original HGVS for reference
                    result['hgvs'] = hgvs

                    # 🎯 Add ClinVar comparison if we have a successful result
                    if result['status'] == 'SUCCESS':
                        our_classification = result.get('final_classification', 'UNKNOWN')
                        agreement_status, agreement_flag = self.determine_agreement(our_classification, expected_clinvar)

                        result['expected_clinvar'] = expected_clinvar
                        result['agreement_status'] = agreement_status
                        result['agreement_flag'] = agreement_flag

                        # Update agreement stats
                        if agreement_status == 'AGREE':
                            stats['agreements'] += 1
                        elif agreement_status == 'BETTER_DATA_BENIGN':
                            stats['better_data_benign'] += 1
                        elif agreement_status == 'BETTER_DATA_PATHOGENIC':
                            stats['better_data_pathogenic'] += 1
                        elif agreement_status == 'DISAGREE':
                            stats['disagreements'] += 1
                        else:
                            stats['unclear'] += 1

                        print(f"   🎯 Our: {our_classification} vs Expected: {expected_clinvar} → {agreement_flag}")

                    results.append(result)

                    # Update stats (handle both old and new cascade formats)
                    if result['status'] == 'SUCCESS':
                        stats['success'] += 1

                        # Check analyzers run for new biological routing format
                        analyzers_run = result.get('analyzers_run', [])
                        if len(analyzers_run) > 1:
                            stats['cascade_triggered'] += 1
                        else:
                            stats['dn_only'] += 1
                    else:
                        stats['failed'] += 1
                    
                    stats['processed'] += 1
        
        except Exception as e:
            print(f"❌ Error reading CSV: {e}")
            return {'error': str(e)}
        
        # Write results to TSV
        self.write_results_tsv(results, output_path)

        # 🎯 Calculate and print performance metrics
        self.print_performance_summary(stats, results)

        return {'stats': stats, 'results_file': output_path}

    def print_performance_summary(self, stats: Dict, results: List[Dict]):
        """Print comprehensive performance summary with ClinVar comparison"""

        # Basic processing stats
        print("\n" + "=" * 60)
        print("📊 CASCADE PROCESSING SUMMARY")
        print("=" * 60)
        for key, value in stats.items():
            if not key.startswith(('agreements', 'better_data', 'disagreements', 'unclear')):
                print(f"{key.replace('_', ' ').title()}: {value}")

        # 🎯 ClinVar Performance Analysis
        print("\n" + "🎯" * 20)
        print("🎯 CLINVAR PERFORMANCE ANALYSIS")
        print("🎯" * 20)

        total_analyzed = stats['success']
        if total_analyzed > 0:
            agreements = stats['agreements']
            better_benign = stats['better_data_benign']
            better_pathogenic = stats['better_data_pathogenic']
            disagreements = stats['disagreements']
            unclear = stats['unclear']

            print(f"✅ Perfect Agreements: {agreements}/{total_analyzed} ({agreements/total_analyzed*100:.1f}%)")
            print(f"🎯 Better Data (Benign): {better_benign}/{total_analyzed} ({better_benign/total_analyzed*100:.1f}%)")
            print(f"🎯 Better Data (Pathogenic): {better_pathogenic}/{total_analyzed} ({better_pathogenic/total_analyzed*100:.1f}%)")
            print(f"❌ Disagreements: {disagreements}/{total_analyzed} ({disagreements/total_analyzed*100:.1f}%)")
            print(f"❓ Unclear: {unclear}/{total_analyzed} ({unclear/total_analyzed*100:.1f}%)")

            # Overall success rate (agreements + better data)
            total_success = agreements + better_benign + better_pathogenic
            success_rate = total_success / total_analyzed * 100

            print(f"\n🚀 OVERALL SUCCESS RATE: {total_success}/{total_analyzed} ({success_rate:.1f}%)")
            print(f"   (Includes perfect agreements + cases where we provide better data)")

            # Performance by classification
            print(f"\n📈 PERFORMANCE BREAKDOWN:")
            self.print_classification_breakdown(results)
        else:
            print("⚠️  No successful analyses to evaluate")

    def print_classification_breakdown(self, results: List[Dict]):
        """Print performance breakdown by our classification"""

        # Group results by our classification
        by_classification = {}
        for result in results:
            if result.get('status') == 'SUCCESS':
                our_class = result.get('final_classification', 'UNKNOWN')
                if our_class not in by_classification:
                    by_classification[our_class] = []
                by_classification[our_class].append(result)

        for our_class, class_results in by_classification.items():
            if not class_results:
                continue

            total = len(class_results)
            agreements = len([r for r in class_results if r.get('agreement_flag') == '✅'])
            better_data = len([r for r in class_results if r.get('agreement_flag') == '🎯'])
            disagreements = len([r for r in class_results if r.get('agreement_flag') == '❌'])

            success_rate = (agreements + better_data) / total * 100

            print(f"   {our_class}: {agreements + better_data}/{total} ({success_rate:.1f}%) success")
            print(f"      ✅ {agreements} agreements, 🎯 {better_data} better data, ❌ {disagreements} disagreements")
    
    def write_results_tsv(self, results: List[Dict], output_path: str):
        """Write cascade results to TSV file"""
        if not results:
            print("⚠️  No results to write")
            return
        
        # Define columns for biological cascade output with ClinVar comparison
        columns = [
            'gene', 'variant', 'hgvs', 'gnomad_freq', 'status',
            'routing_strategy', 'routing_confidence', 'routing_rationale',
            'analyzers_run', 'summary',
            'dn_score', 'lof_score', 'gof_score',
            'dn_classification', 'lof_classification', 'gof_classification',
            'final_score', 'final_classification', 'explanation',
            # 🎯 ClinVar comparison columns
            'expected_clinvar', 'agreement_status', 'agreement_flag'
        ]
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
            writer.writeheader()
            
            for result in results:
                row = {}
                routing_info = result.get('routing_info', {})

                for col in columns:
                    if col == 'analyzers_run':
                        row[col] = ', '.join(result.get('analyzers_run', []))
                    elif col.endswith('_score'):
                        analyzer = col.replace('_score', '').upper()
                        # Use plausibility-filtered scores if available, otherwise raw scores
                        filtered_scores = result.get('plausibility_filtered_scores', {})
                        if filtered_scores and analyzer in filtered_scores:
                            row[col] = filtered_scores[analyzer]
                        else:
                            row[col] = result.get('scores', {}).get(analyzer, '')
                    elif col.endswith('_classification'):
                        analyzer = col.replace('_classification', '').upper()
                        row[col] = result.get('classifications', {}).get(analyzer, '')
                    elif col.startswith('routing_'):
                        routing_key = col.replace('routing_', '')
                        row[col] = routing_info.get(routing_key, '')
                    else:
                        row[col] = result.get(col, '')
                
                writer.writerow(row)
        
        print(f"💾 CASCADE results saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="🌊 CASCADE BATCH PROCESSOR - Multi-Analyzer CSV Processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process ClinVar export through full cascade
  python3 cascade_batch_processor.py --input clinvar_variants.csv --output cascade_results.tsv
  
  # Custom frequency threshold
  python3 cascade_batch_processor.py --input variants.csv --output results.tsv --freq-filter 0.001
        """
    )
    
    parser.add_argument('--input', required=True, help='Input CSV file path')
    parser.add_argument('--output', required=True, help='Output TSV file path')
    parser.add_argument('--freq-filter', type=float, default=0.01, 
                       help='Skip variants with gnomAD frequency > threshold (default: 0.01)')
    parser.add_argument('--alphafold-path', default="/mnt/Arcana/alphafold_human/structures/",
                       help='Path to AlphaFold structures directory')
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        print(f"❌ Input file not found: {args.input}")
        return 1
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Process CSV through cascade
    processor = CascadeBatchProcessor(args.alphafold_path)
    result = processor.process_csv(
        args.input, 
        args.output,
        args.freq_filter
    )
    
    if 'error' in result:
        print(f"❌ Processing failed: {result['error']}")
        return 1
    
    print(f"\n✅ CASCADE processing complete! Results in {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
