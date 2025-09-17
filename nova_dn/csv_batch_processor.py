#!/usr/bin/env python3
"""
CSV Batch Processor for DN Analysis 🧬⚡
Generic tool for processing variant CSV files through the Nova DN analyzer

Supports:
- ClinVar Miner exports
- Custom variant datasets  
- Population frequency filtering
- DN mechanism filtering
- Real AlphaFold sequences
- Proper error handling

Authors: Ace & Nova (2025)
Usage: python3 -m nova_dn.csv_batch_processor --input variants.csv --output results.tsv
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional, Tuple
from pathlib import Path

from .alphafold_sequence import AlphaFoldSequenceExtractor
from .dn_mechanism_filter_v2 import DNMechanismFilterV2
from .sequence_manager import SequenceManager


class CSVBatchProcessor:
    """Process CSV files of variants through DN analysis pipeline"""
    
    def __init__(self, alphafold_path: str = "/mnt/Arcana/alphafold_human/structures/"):
        self.extractor = AlphaFoldSequenceExtractor(alphafold_path)
        self.dn_filter = DNMechanismFilterV2()
        self.sequence_manager = SequenceManager()
        self.temp_files = []  # Track temp files for cleanup
        
    def parse_hgvs(self, hgvs: str) -> Optional[Dict]:
        """
        Parse HGVS notation to extract gene, transcript, and variant
        
        Examples:
        - NM_000540.3(RYR1):c.1077T>C (p.Ala359=)
        - NM_000546.6(TP53):c.817C>T (p.Arg273His)
        """
        # Extract gene name from parentheses
        gene_match = re.search(r'\(([^)]+)\)', hgvs)
        if not gene_match:
            return None
        gene = gene_match.group(1)
        
        # Look for protein change p.RefPosAlt
        protein_match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs)
        if protein_match:
            ref_long, pos, alt_long = protein_match.groups()
            # Convert 3-letter to 1-letter codes
            aa_map = {
                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
            }
            ref = aa_map.get(ref_long, ref_long)
            alt = aa_map.get(alt_long, alt_long)
            variant = f"p.{ref}{pos}{alt}"
            
            return {
                'gene': gene,
                'variant': variant,
                'position': int(pos),
                'ref': ref,
                'alt': alt,
                'type': 'missense'
            }
        
        # Check for nonsense/frameshift variants p.Arg273Ter or p.Q3485Ter
        nonsense_match = re.search(r'p\.([A-Z][a-z]{0,2})(\d+)Ter', hgvs)
        if nonsense_match:
            aa_code = nonsense_match.group(1)
            position = nonsense_match.group(2)
            return {
                'gene': gene,
                'variant': f'p.{aa_code}{position}Ter',
                'type': 'nonsense',
                'skip_reason': 'SKIPPED_FRAMESHIFT'
            }

        # Check for synonymous variants p.Ala359=
        synonymous_match = re.search(r'p\.([A-Z][a-z]{2})(\d+)=', hgvs)
        if synonymous_match:
            return {
                'gene': gene,
                'variant': f'p.{synonymous_match.group(1)}{synonymous_match.group(2)}=',
                'type': 'synonymous',
                'skip_reason': 'SKIPPED_SYNONYMOUS'
            }

        # Check for intronic variants c.1122+52T>C
        if re.search(r'c\.\d+[+-]\d+', hgvs):
            return {
                'gene': gene,
                'variant': 'intronic',
                'type': 'intronic',
                'skip_reason': 'SKIPPED_INTRONIC'
            }
        
        return None
    
    def get_uniprot_id(self, gene_name: str) -> Optional[str]:
        """Get UniProt ID for gene (using common mappings)"""
        # Common gene -> UniProt mappings
        gene_to_uniprot = {
            'TP53': 'P04637',
            'COL1A1': 'P02452', 
            'FGFR3': 'P22607',
            'VWF': 'P04275',
            'RYR1': 'P21817',
            'FBN1': 'P35555',
            'SCN5A': 'Q14524',
            'KCNQ1': 'P51787',
            'CACNA1I': 'Q9P0X4',
            'BRCA1': 'P38398',
            'BRCA2': 'P51587',
            'MSH2': 'P43246',
            'MLH1': 'P40692',
            'ATM': 'Q13315',
            'PTEN': 'P60484'
        }
        return gene_to_uniprot.get(gene_name)
    
    def process_variant(self, gene: str, variant: str, gnomad_freq: float) -> Dict:
        """Process a single variant through the DN analysis pipeline"""
        
        # Get UniProt ID
        uniprot_id = self.get_uniprot_id(gene)
        if not uniprot_id:
            return {
                'gene': gene,
                'variant': variant,
                'gnomad_freq': gnomad_freq,
                'error': f'No UniProt ID mapping for {gene}',
                'status': 'SKIPPED'
            }
        
        # Check DN likelihood first
        dn_likelihood, dn_evidence = self.dn_filter.assess_dn_likelihood(gene, uniprot_id)
        if dn_likelihood < 0.3:
            return {
                'gene': gene,
                'variant': variant,
                'gnomad_freq': gnomad_freq,
                'dn_likelihood': dn_likelihood,
                'status': 'FILTERED_LOW_DN',
                'reason': 'Gene has low DN likelihood (likely LOF mechanism)'
            }
        
        # Extract variant position for smart sequence selection
        variant_position = None
        if variant.startswith('p.') and len(variant) > 3:
            # Extract position from p.RefPosAlt format
            import re
            pos_match = re.search(r'p\.[A-Z](\d+)[A-Z]', variant)
            if pos_match:
                variant_position = int(pos_match.group(1))

        # Get relevant DN mechanisms for this gene (FINALLY!)
        relevant_mechanisms, mechanism_evidence = self.dn_filter.select_relevant_mechanisms(gene, uniprot_id)
        if not relevant_mechanisms:
            return {
                'gene': gene,
                'variant': variant,
                'gnomad_freq': gnomad_freq,
                'dn_likelihood': dn_likelihood,
                'status': 'FILTERED_NO_DN_MECHANISMS',
                'reason': 'No relevant DN mechanisms for this gene type'
            }

        print(f"🔬 Testing {len(relevant_mechanisms)} relevant mechanisms: {', '.join(relevant_mechanisms)}")

        # Get best sequence (UniProt or AlphaFold based on variant position)
        try:
            sequence, source, temp_fasta_path = self.sequence_manager.get_best_sequence(
                gene, uniprot_id, variant_position
            )
            self.temp_files.append(temp_fasta_path)
            print(f"Using {len(sequence)} residue {source} sequence for {gene} {variant}")

        except Exception as e:
            return {
                'gene': gene,
                'variant': variant,
                'gnomad_freq': gnomad_freq,
                'error': f'Could not get sequence: {e}',
                'status': 'SKIPPED'
            }
        
        # Run DN analyzer with the magic incantation
        cmd = [
            sys.executable, '-m', 'nova_dn.analyzer',
            '--seq-file', temp_fasta_path,
            '--variant', variant,
            '--annotations-json', 'resources/protein_annotations.json',
            '--protein', gene,
            '--json'
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if result.returncode != 0:
                return {
                    'gene': gene,
                    'variant': variant,
                    'gnomad_freq': gnomad_freq,
                    'error': result.stderr.strip(),
                    'status': 'ANALYSIS_FAILED'
                }
            
            # Parse JSON output
            analysis_result = json.loads(result.stdout)
            
            # Add metadata
            analysis_result.update({
                'gene': gene,
                'gnomad_freq': gnomad_freq,
                'dn_likelihood': dn_likelihood,
                'uniprot_id': uniprot_id,
                'status': 'SUCCESS'
            })
            
            return analysis_result
            
        except subprocess.TimeoutExpired:
            return {
                'gene': gene,
                'variant': variant,
                'gnomad_freq': gnomad_freq,
                'error': 'Analysis timeout (>30s)',
                'status': 'TIMEOUT'
            }
        except json.JSONDecodeError as e:
            return {
                'gene': gene,
                'variant': variant,
                'gnomad_freq': gnomad_freq,
                'error': f'JSON decode error: {e}',
                'status': 'JSON_ERROR'
            }
        except Exception as e:
            return {
                'gene': gene,
                'variant': variant,
                'gnomad_freq': gnomad_freq,
                'error': f'Unexpected error: {e}',
                'status': 'ERROR'
            }
    
    def process_csv(self, input_path: str, output_path: str, 
                   freq_threshold: float = 0.01, use_dn_filter: bool = True) -> Dict:
        """
        Process CSV file of variants
        
        Args:
            input_path: Path to input CSV file
            output_path: Path to output TSV file
            freq_threshold: Skip variants with frequency > threshold
            use_dn_filter: Whether to apply DN mechanism filtering
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
            'skipped_dn_filter': 0,
            'failed': 0,
            'success': 0
        }
        
        print(f"🧬 Processing variants from {input_path}")
        print(f"📊 Frequency threshold: {freq_threshold}")
        print(f"🔍 DN filtering: {'ON' if use_dn_filter else 'OFF'}")
        print("=" * 50)
        
        try:
            with open(input_path, 'r') as f:
                reader = csv.DictReader(f)
                
                for row in reader:
                    stats['total_variants'] += 1
                    
                    # Parse HGVS
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
                    
                    # Parse variant
                    parsed = self.parse_hgvs(hgvs)
                    if not parsed:
                        stats['skipped_unparseable'] += 1
                        continue

                    # Skip non-missense variants with specific reasons
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
                    
                    print(f"Processing {gene} {variant} (freq: {gnomad_freq:.4f})")
                    
                    # Process variant
                    result = self.process_variant(gene, variant, gnomad_freq)
                    results.append(result)
                    
                    # Update stats
                    if result['status'] == 'SUCCESS':
                        stats['success'] += 1
                    elif result['status'] == 'FILTERED_LOW_DN':
                        stats['skipped_dn_filter'] += 1
                    else:
                        stats['failed'] += 1
                    
                    stats['processed'] += 1
        
        except Exception as e:
            print(f"❌ Error reading CSV: {e}")
            return {'error': str(e)}
        
        # Write results to TSV
        self.write_results_tsv(results, output_path)
        
        # Print summary
        print("\n" + "=" * 50)
        print("📊 PROCESSING SUMMARY")
        print("=" * 50)
        for key, value in stats.items():
            print(f"{key.replace('_', ' ').title()}: {value}")
        
        # Cleanup temp files
        self.cleanup()
        
        return {'stats': stats, 'results_file': output_path}
    
    def write_results_tsv(self, results: List[Dict], output_path: str):
        """Write results to TSV file"""
        if not results:
            print("⚠️  No results to write")
            return
        
        # Define columns
        columns = [
            'gene', 'variant', 'gnomad_freq', 'status', 'dn_likelihood',
            'interface_poisoning', 'active_site_jamming', 'lattice_disruption', 
            'trafficking_maturation', 'top_mechanism', 'explanation', 'error'
        ]
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
            writer.writeheader()
            
            for result in results:
                row = {}
                for col in columns:
                    if col in ['interface_poisoning', 'active_site_jamming', 
                              'lattice_disruption', 'trafficking_maturation']:
                        # Extract mechanism scores
                        scores = result.get('mechanism_scores', {})
                        row[col] = scores.get(col, '')
                    else:
                        row[col] = result.get(col, '')
                
                writer.writerow(row)
        
        print(f"💾 Results saved to {output_path}")
    
    def cleanup(self):
        """Clean up temporary files"""
        for temp_file in self.temp_files:
            try:
                os.unlink(temp_file)
            except:
                pass
        self.temp_files = []


def main():
    parser = argparse.ArgumentParser(
        description="CSV Batch Processor for DN Analysis 🧬⚡",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process ClinVar export with default settings
  python3 -m nova_dn.csv_batch_processor --input clinvar_benign.csv --output results.tsv
  
  # Custom frequency threshold and no DN filtering
  python3 -m nova_dn.csv_batch_processor --input variants.csv --output results.tsv --freq-filter 0.001 --no-dn-filter
        """
    )
    
    parser.add_argument('--input', required=True, help='Input CSV file path')
    parser.add_argument('--output', required=True, help='Output TSV file path')
    parser.add_argument('--freq-filter', type=float, default=0.01, 
                       help='Skip variants with gnomAD frequency > threshold (default: 0.01)')
    parser.add_argument('--no-dn-filter', action='store_true',
                       help='Disable DN mechanism filtering')
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
    
    # Process CSV
    processor = CSVBatchProcessor(args.alphafold_path)
    result = processor.process_csv(
        args.input, 
        args.output,
        args.freq_filter,
        not args.no_dn_filter
    )
    
    if 'error' in result:
        print(f"❌ Processing failed: {result['error']}")
        return 1
    
    print(f"\n✅ Processing complete! Results in {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
