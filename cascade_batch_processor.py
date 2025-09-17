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

from cascade_analyzer import CascadeAnalyzer
from nova_dn.csv_batch_processor import CSVBatchProcessor


class CascadeBatchProcessor:
    """Process CSV files through the complete cascade system"""
    
    def __init__(self, alphafold_path: str = "/mnt/Arcana/alphafold_human/structures/"):
        self.cascade_analyzer = CascadeAnalyzer(alphafold_path)
        self.csv_processor = CSVBatchProcessor(alphafold_path)  # For HGVS parsing
        
    def process_csv(self, input_path: str, output_path: str, 
                   freq_threshold: float = 0.01) -> Dict:
        """
        Process CSV file through cascade analysis
        
        Args:
            input_path: Path to input CSV file
            output_path: Path to output TSV file  
            freq_threshold: Skip variants with frequency > threshold
        """
        
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
            'success': 0
        }
        
        print(f"🌊 Processing variants through CASCADE SYSTEM")
        print(f"📊 Frequency threshold: {freq_threshold}")
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
        
        # Print summary
        print("\n" + "=" * 60)
        print("📊 CASCADE PROCESSING SUMMARY")
        print("=" * 60)
        for key, value in stats.items():
            print(f"{key.replace('_', ' ').title()}: {value}")
        
        return {'stats': stats, 'results_file': output_path}
    
    def write_results_tsv(self, results: List[Dict], output_path: str):
        """Write cascade results to TSV file"""
        if not results:
            print("⚠️  No results to write")
            return
        
        # Define columns for biological cascade output
        columns = [
            'gene', 'variant', 'hgvs', 'gnomad_freq', 'status',
            'routing_strategy', 'routing_confidence', 'routing_rationale',
            'analyzers_run', 'summary',
            'dn_score', 'lof_score', 'gof_score',
            'dn_classification', 'lof_classification', 'gof_classification',
            'final_score', 'final_classification', 'explanation'
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
