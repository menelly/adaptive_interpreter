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

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.gnomad_frequency_fetcher import GnomADFrequencyFetcher
from cascade.cascade_analyzer import CascadeAnalyzer
from nova_dn.csv_batch_processor import CSVBatchProcessor
from nova_dn.mixed_mechanism_resolver import UnifiedMechanismResolver


class CascadeBatchProcessor:
    """Process CSV files through the complete cascade system with ClinVar comparison"""

    def __init__(self, alphafold_path: str = "/mnt/Arcana/alphafold_human/structures/",
                 override_family: str = None, conservative_mode: bool = False):
        self.cascade_analyzer = CascadeAnalyzer(alphafold_path, override_family, conservative_mode)
        self.csv_processor = CSVBatchProcessor(alphafold_path)  # For HGVS parsing

        # 🧬 NOVA'S UNIFIED MECHANISM RESOLVER
        self.unified_resolver = UnifiedMechanismResolver()

        # 🔥 REAL FREQUENCY FETCHER (No more hardcoding!)
        self.frequency_fetcher = GnomADFrequencyFetcher()
        self.frequency_fetcher.load_cache()  # Load cached frequencies

        # Gene to chromosome mapping (common genes)
        self.gene_to_chromosome = {
            'FBN1': 'chr15',
            'KCNMA1': 'chr10',
            'TP53': 'chr17',
            'BRCA1': 'chr17',
            'BRCA2': 'chr13',
            'COL1A1': 'chr17',
            'COL1A2': 'chr7',
            'CFTR': 'chr7',
            'SCN5A': 'chr3',
            'CACNA1C': 'chr12',
            'KCNH2': 'chr7',
            'MYO7A': 'chr11',
            'RYR1': 'chr19',
            'TTN': 'chr2',
            'MYBPC3': 'chr11',
            'MYH7': 'chr14',
            'LDLR': 'chr19',
            'APOB': 'chr2',
            'PCSK9': 'chr1'
        }

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

    def extract_clinical_significance(self, clinvar_string: str) -> str:
        """
        Extract the actual clinical significance from messy ClinVar strings using regex

        🔥 REN'S BRILLIANT IDEA: Look for keywords, ignore the disease names!
        """
        import re

        if not clinvar_string:
            return 'UNKNOWN'

        # Convert to uppercase for matching
        text = clinvar_string.upper()

        # 🎯 PRIORITY ORDER: More specific matches first
        patterns = [
            (r'\bLIKELY\s+PATHOGENIC\b', 'LIKELY_PATHOGENIC'),
            (r'\bLIKELY\s+BENIGN\b', 'LIKELY_BENIGN'),
            (r'\bPATHOGENIC\b', 'PATHOGENIC'),
            (r'\bBENIGN\b', 'BENIGN'),
            (r'\bUNCERTAIN\s+SIGNIFICANCE\b', 'VUS'),
            (r'\bVUS\b', 'VUS'),
            (r'\bNOT\s+PROVIDED\b', 'UNKNOWN'),
            (r'\bCONFLICTING\s+INTERPRETATIONS\b', 'CONFLICTING'),
        ]

        # Find the first match (most specific)
        for pattern, classification in patterns:
            if re.search(pattern, text):
                return classification

        return 'UNKNOWN'

    def determine_agreement(self, our_classification: str, expected_clinvar: str) -> tuple:
        """
        Determine agreement status between our prediction and ClinVar expectation

        Returns: (agreement_status, agreement_flag)
        """
        our_class = our_classification.upper()

        # 🔥 EXTRACT ACTUAL CLINICAL SIGNIFICANCE from messy ClinVar string
        expected = self.extract_clinical_significance(expected_clinvar)

        # 🧬 REN'S BRILLIANT INSIGHT: 1 step away = same biological conclusion!

        # Define biological families
        benign_family = ['B', 'LB']
        pathogenic_family = ['LP', 'P']
        uncertain_family = ['VUS', 'VUS-P']

        clinvar_benign = ['BENIGN', 'LIKELY_BENIGN']
        clinvar_pathogenic = ['PATHOGENIC', 'LIKELY_PATHOGENIC']
        clinvar_uncertain = ['VUS']

        # ✅ AGREEMENT cases (same biological family)
        if (our_class in benign_family and expected in clinvar_benign) or \
           (our_class in pathogenic_family and expected in clinvar_pathogenic) or \
           (our_class in uncertain_family and expected in clinvar_uncertain):
            return ('AGREE', '✅')

        # 🎯 BETTER DATA cases (moving OUTWARD from uncertainty - GOOD!)
        elif (our_class in benign_family and expected in clinvar_uncertain):
            return ('BETTER_DATA_BENIGN', '🎯')  # VUS → LB/B (found benign evidence!)
        elif (our_class in pathogenic_family and expected in clinvar_uncertain):
            return ('BETTER_DATA_PATHOGENIC', '🎯')  # VUS → LP/P (found pathogenic evidence!)

        # 🩺 CLINICAL CORRELATION NEEDED (one step off - math vs clinical judgment)
        elif ((our_class == 'VUS-P' and expected in clinvar_pathogenic) or
              (our_class in pathogenic_family and expected == 'VUS-P') or
              (our_class == 'VUS' and expected in clinvar_benign) or
              (our_class in benign_family and expected == 'VUS')):
            return ('CLINICAL_CORRELATION', '🩺')  # Math is iffy but clinical judgment may differ

        # ❌ REAL DISAGREEMENTS - Only major flips across the neutral line (B/LB ↔ LP/P)
        elif (our_class in benign_family and expected in clinvar_pathogenic) or \
             (our_class in pathogenic_family and expected in clinvar_benign):
            return ('DISAGREE', '❌')  # MAJOR FLIP - True disagreement!

        # 📊 CORRELATIONS - Adjacent classifications (refinements, not errors)
        # Moving to/from VUS is correlation, not disagreement
        elif (our_class in uncertain_family and expected in clinvar_benign):
            return ('CLINICAL_CORRELATION', '🩺')  # LB/B → VUS (being more cautious)
        elif (our_class in uncertain_family and expected in clinvar_pathogenic):
            return ('CLINICAL_CORRELATION', '🩺')  # LP/P → VUS (being more cautious)

        # Handle conflicting interpretations as unclear
        elif expected == 'CONFLICTING':
            return ('UNCLEAR', '❓')

        # Unknown/unclear cases
        else:
            return ('UNCLEAR', '❓')

    def run_nova_analysis(self, gene: str, variant: str, cascade_result: Dict) -> Dict:
        """
        🧬 RUN NOVA'S UNIFIED MECHANISM RESOLVER

        Extract sequence and position info from cascade result and run
        Nova's Revolutionary Framework for comparison.
        """
        # Get sequence from cascade result or fetch it
        sequence = cascade_result.get('sequence', '')
        uniprot_id = cascade_result.get('uniprot_id', '')

        if not sequence and uniprot_id:
            # Extract position for sequence selection
            import re
            pos_match = re.search(r'p\.[A-Z](\d+)[A-Z]', variant)
            variant_position = int(pos_match.group(1)) if pos_match else None

            try:
                sequence, source, temp_fasta_path = self.cascade_analyzer.sequence_manager.get_best_sequence(
                    gene, uniprot_id, variant_position
                )
                # Sequence obtained successfully
            except Exception as e:
                raise ValueError(f"Could not get sequence for {gene}: {e}")
        elif not sequence:
            raise ValueError(f"No sequence or UniProt ID available for {gene}")

        # Parse variant (e.g., "p.P840L" or "P840L" -> pos=840, ref="P", alt="L")
        import re

        # Remove "p." prefix if present
        clean_variant = variant.replace('p.', '') if variant.startswith('p.') else variant

        match = re.match(r'([A-Z])(\d+)([A-Z*])', clean_variant)
        if not match:
            raise ValueError(f"Cannot parse variant: {variant} (cleaned: {clean_variant})")

        ref, pos_str, alt = match.groups()
        pos1 = int(pos_str)

        # Build context from cascade result
        context = {
            'gene_name': gene,
            'gene_family': self.get_gene_family(gene),
            'domains': cascade_result.get('domains', []),
            'uniprot_id': cascade_result.get('uniprot_id', ''),
        }

        # Determine variant type (simplified for now)
        variant_type = "nonsense" if alt == "*" else "missense"

        # Running Nova's analysis

        # Run Nova's Revolutionary Framework!
        nova_result = self.unified_resolver.resolve_mechanisms(
            sequence, pos1, ref, alt, context, variant_type
        )

        return nova_result

    def get_gene_family(self, gene: str) -> str:
        """
        🧬 GET GENE FAMILY FOR NOVA'S ANALYSIS WITH INHERITANCE AWARENESS

        Map gene names to protein families for Nova's specialized analyzers.
        Now includes inheritance pattern detection for proper mechanism weighting!
        """
        # 🧬 INHERITANCE-AWARE GENE CLASSIFICATION
        inheritance_patterns = {
            # Autosomal Recessive genes (LOF boost!)
            'CFTR': 'AUTOSOMAL_RECESSIVE',  # Cystic fibrosis - LOF mechanism!
            'FKRP': 'AUTOSOMAL_RECESSIVE',  # Muscular dystrophy - LOF mechanism!
            'ABCA4': 'AUTOSOMAL_RECESSIVE', # Stargardt disease - LOF mechanism!

            # Autosomal Dominant structural genes (DN boost!)
            'COL1A1': 'COLLAGEN_FIBRILLAR',  # OI - DN mechanism!
            'COL1A2': 'COLLAGEN_FIBRILLAR',  # OI - DN mechanism!
            'COL3A1': 'COLLAGEN_FIBRILLAR',  # EDS - DN mechanism!
            'FBN1': 'FIBRILLIN',            # Marfan - DN mechanism!

            # Ion channels (mixed mechanisms)
            'SCN1A': 'ION_CHANNEL',
            'KCNQ2': 'ION_CHANNEL',
            'RYR1': 'ION_CHANNEL',
            'KCNMA1': 'ION_CHANNEL',
            'SCN5A': 'ION_CHANNEL',
            'CACNA1C': 'ION_CHANNEL',

            # Tumor suppressors
            'TP53': 'TUMOR_SUPPRESSOR',
            'BRCA1': 'TUMOR_SUPPRESSOR',
            'BRCA2': 'TUMOR_SUPPRESSOR',
        }

        return inheritance_patterns.get(gene, 'GENERAL')

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
            'clinical_correlation': 0,
            'disagreements': 0,
            'unclear': 0
        }

        print(f"🌊 Processing variants through CASCADE SYSTEM")
        print(f"📊 Frequency threshold: {freq_threshold}")
        print(f"🎯 Expected ClinVar classification: {expected_clinvar}")
        print("=" * 60)
        
        try:
            with open(input_path, 'r') as f:
                # Detect delimiter (TSV vs CSV)
                delimiter = '\t' if input_path.endswith('.tsv') else ','
                reader = csv.DictReader(f, delimiter=delimiter)
                delimiter_name = 'TAB' if delimiter == '\t' else 'COMMA'
                print(f"📄 Using delimiter: {delimiter_name}")
                
                for row in reader:
                    stats['total_variants'] += 1

                    # 🧬 NEW: Handle Ren's TSV format with coordinates!
                    if 'Chrpos' in row and 'AA Chg' in row:
                        # Ren's format: Clinical significance and condition | Chrpos | Variation Name | AA Chg
                        clinical_sig = row.get('Clinical significance and condition', '')
                        chrpos = row.get('Chrpos', '')
                        variation_name = row.get('Variation Name', '')
                        aa_chg = row.get('AA Chg', '')

                        # Extract gene from variation name (e.g., "NM_000138.5(FBN1):c.5938G>C")
                        import re
                        gene_match = re.search(r'\(([A-Z0-9]+)\)', variation_name)
                        if not gene_match:
                            stats['skipped_unparseable'] += 1
                            print(f"⏭️  Could not extract gene from: {variation_name}")
                            continue

                        gene = gene_match.group(1)

                        # Skip if no protein change (synonymous or intronic)
                        if not aa_chg or aa_chg.strip() == '':
                            stats['skipped_synonymous'] += 1
                            print(f"⏭️  Skipping {gene} (no protein change)")
                            continue

                        # Convert AA change format (E1980Q -> p.Glu1980Gln)
                        if not aa_chg.startswith('p.'):
                            # Convert single letter to p. format
                            aa_match = re.match(r'([A-Z])(\d+)([A-Z])', aa_chg)
                            if aa_match:
                                variant = f"p.{aa_match.group(1)}{aa_match.group(2)}{aa_match.group(3)}"
                            else:
                                stats['skipped_unparseable'] += 1
                                print(f"⏭️  Could not parse AA change: {aa_chg}")
                                continue
                        else:
                            variant = aa_chg

                        # Parse coordinates (format: "48,444,640(-)" -> need to add chromosome)
                        if chrpos:
                            # Remove commas and parentheses
                            pos_clean = chrpos.replace(',', '').replace('(', '').replace(')', '').replace('-', '').replace('+', '')
                            try:
                                pos = int(pos_clean)
                                # Get chromosome from mapping
                                chrom = self.gene_to_chromosome.get(gene, 'chr1')  # Default to chr1 if unknown

                                coord_key = f"{gene}_{variant}"

                                # 🔥 ADD TO CONSERVATION SYSTEM!
                                if hasattr(self.cascade_analyzer, '_get_conservation_multiplier'):
                                    # Update the known coordinates in the cascade analyzer
                                    if not hasattr(self.cascade_analyzer, '_temp_coordinates'):
                                        self.cascade_analyzer._temp_coordinates = {}
                                    self.cascade_analyzer._temp_coordinates[coord_key] = (chrom, pos)
                                    print(f"🎯 Added coordinates for {gene} {variant}: {chrom}:{pos}")
                            except ValueError:
                                print(f"⚠️  Could not parse position from: {chrpos}")
                                continue

                        hgvs = variation_name  # For reference

                        # 🔥 GET REAL FREQUENCY (NO FAKE ESTIMATES!)
                        clinical_significance = row.get('Clinical significance and condition', '')
                        print(f"🔍 Looking up real gnomAD frequency...")
                        freq_result = self.frequency_fetcher.get_frequency_from_clinvar_row(chrpos, variation_name, "")  # No clinical significance = no fake estimates
                        gnomad_freq = freq_result['frequency']

                        if freq_result['error']:
                            print(f"⚠️ Frequency fetch issue: {freq_result['error']} (using {gnomad_freq})")
                        else:
                            print(f"✅ Real gnomAD frequency: {gnomad_freq:.6f} (source: {freq_result['source']})")

                    else:
                        # Original format: HGVS | gnomAD frequency
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
                    
                    # 🚀 Run BIOLOGICAL CASCADE analysis (clean single pipeline)
                    result = self.cascade_analyzer.analyze_cascade_biological(
                        gene, variant, gnomad_freq, 'missense'  # Default to missense
                    )
                    result['hgvs'] = hgvs

                    # 🎯 Add ClinVar comparison if we have a successful result
                    if result['status'] == 'SUCCESS':
                        our_classification = result.get('final_classification', 'UNKNOWN')

                        # Use actual clinical significance from TSV if available
                        if 'Chrpos' in row:
                            expected_clinvar = clinical_sig  # From the TSV row

                        agreement_status, agreement_flag = self.determine_agreement(our_classification, expected_clinvar)

                        result['expected_clinvar'] = expected_clinvar
                        result['agreement_status'] = agreement_status
                        result['agreement_flag'] = agreement_flag

                        # Add conservation info if available
                        if hasattr(self.cascade_analyzer, '_temp_coordinates'):
                            result['coordinates'] = chrpos if 'Chrpos' in row else ''

                        # Update agreement stats
                        if agreement_status == 'AGREE':
                            stats['agreements'] += 1
                        elif agreement_status == 'BETTER_DATA_BENIGN':
                            stats['better_data_benign'] += 1
                        elif agreement_status == 'BETTER_DATA_PATHOGENIC':
                            stats['better_data_pathogenic'] += 1
                        elif agreement_status == 'CLINICAL_CORRELATION':
                            stats['clinical_correlation'] += 1
                        elif agreement_status == 'DISAGREE':
                            stats['disagreements'] += 1
                        else:
                            stats['unclear'] += 1

                        # 🔥 QUICK FIX: Show parsed clinical significance instead of raw messy string
                        parsed_expected = self.extract_clinical_significance(expected_clinvar)
                        print(f"   🎯 Our: {our_classification} vs Expected: {parsed_expected} → {agreement_flag}")

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

        # 🧬 Print human-friendly summary
        self.print_human_friendly_summary(results)

        # 💾 Save frequency cache for future runs
        self.frequency_fetcher.save_cache()

        return {'stats': stats, 'results_file': output_path}

    def print_performance_summary(self, stats: Dict, results: List[Dict]):
        """Print comprehensive performance summary with ClinVar comparison"""

        # Basic processing stats
        print("\n" + "=" * 60)
        print("📊 CASCADE PROCESSING SUMMARY")
        print("=" * 60)
        for key, value in stats.items():
            if not key.startswith(('agreements', 'better_data', 'clinical_correlation', 'disagreements', 'unclear')):
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
            clinical_correlation = stats['clinical_correlation']
            disagreements = stats['disagreements']
            unclear = stats['unclear']

            print(f"✅ Perfect Agreements: {agreements}/{total_analyzed} ({agreements/total_analyzed*100:.1f}%)")
            print(f"🎯 Better Data (Benign): {better_benign}/{total_analyzed} ({better_benign/total_analyzed*100:.1f}%)")
            print(f"🎯 Better Data (Pathogenic): {better_pathogenic}/{total_analyzed} ({better_pathogenic/total_analyzed*100:.1f}%)")
            print(f"🩺 Clinical Correlation Needed: {clinical_correlation}/{total_analyzed} ({clinical_correlation/total_analyzed*100:.1f}%)")
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

            # 📝 Create disagree lists for review
            self.create_disagree_lists(results, output_path)
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

    def create_disagree_lists(self, results: List[Dict], output_path: str):
        """Create separate files for disagreements and clinical correlation cases"""

        # Separate results by agreement status
        disagreements = []
        clinical_correlation = []

        for result in results:
            if result.get('status') == 'SUCCESS':
                agreement_flag = result.get('agreement_flag', '')
                if agreement_flag == '❌':
                    disagreements.append(result)
                elif agreement_flag == '🩺':
                    clinical_correlation.append(result)

        # Write disagreements file
        if disagreements:
            disagree_path = output_path.replace('.tsv', '_DISAGREEMENTS.tsv')
            self.write_results_tsv(disagreements, disagree_path)
            print(f"\n❌ DISAGREEMENTS written to: {disagree_path}")
            print(f"   {len(disagreements)} variants need investigation")

        # Write clinical correlation file
        if clinical_correlation:
            clinical_path = output_path.replace('.tsv', '_CLINICAL_CORRELATION.tsv')
            self.write_results_tsv(clinical_correlation, clinical_path)
            print(f"\n🩺 CLINICAL CORRELATION cases written to: {clinical_path}")
            print(f"   {len(clinical_correlation)} variants: math vs clinical judgment")

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
                    elif col == 'final_score' or col == 'final_classification':
                        # 🔥 FIX: Handle final_score and final_classification directly from result dict
                        row[col] = result.get(col, '')
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

    def print_human_friendly_summary(self, results: List[Dict]):
        """Print clean human-readable summary of results"""

        print("\n" + "🧬" * 30)
        print("🧬 HUMAN-FRIENDLY VARIANT ANALYSIS RESULTS")
        print("🧬" * 30)
        print(f"{'Gene':8} {'Variant':12} | {'Mechanism Scores':25} | {'Result':15} | ClinVar")
        print("-" * 80)

        for result in results:
            if result.get('status') != 'SUCCESS':
                continue

            gene = result.get('gene', 'Unknown')
            variant = result.get('variant', 'Unknown')

            # Get scores
            dn_score = result.get('dn_score', 0.0)
            lof_score = result.get('lof_score', 0.0)
            gof_score = result.get('gof_score', 0.0)

            # Format scores
            score_parts = []
            if dn_score > 0:
                score_parts.append(f"DN:{dn_score:.2f}")
            if lof_score > 0:
                score_parts.append(f"LOF:{lof_score:.2f}")
            if gof_score > 0:
                score_parts.append(f"GOF:{gof_score:.2f}")

            score_summary = " | ".join(score_parts) if score_parts else "No scores"

            # Get final result
            final_class = result.get('final_classification', 'UNKNOWN')
            final_score = result.get('final_score', 0.0)

            # Add indicators
            indicators = ""
            if 'hotspot_info' in result and result['hotspot_info']:
                hotspot = result['hotspot_info']
                indicators += f" 🔥 {hotspot.get('hotspot_type', 'hotspot')}"

            if result.get('conservation_multiplier_applied', 1.0) > 1.0:
                mult = result['conservation_multiplier_applied']
                indicators += f" 🧬 {mult:.1f}x"

            # ClinVar comparison
            clinvar = result.get('expected_clinvar', 'Unknown')
            agreement = result.get('agreement_flag', '')

            # Format line
            result_col = f"{final_class:6} ({final_score:.3f}){indicators}"
            print(f"{gene:8} {variant:12} | {score_summary:25} | {result_col:15} | {clinvar} {agreement}")

        print("-" * 80)
        print("🔥 = Hotspot detected | 🧬 = Conservation boost applied")
        print("B=Benign, LB=Likely Benign, VUS=Uncertain, VUS-P=VUS-Pathogenic, LP=Likely Pathogenic, P=Pathogenic")
        print("✅ = Agreement | 🎯 = Better data | 🩺 = Clinical correlation | ❌ = Disagreement")


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
    parser.add_argument('--override-family', type=str,
                       help='Override gene family classification (e.g., TUMOR_SUPPRESSOR, ION_CHANNEL, COLLAGEN_FIBRILLAR)')
    parser.add_argument('--conservative-mode', action='store_true',
                       help='Enable conservative classification: downgrade P/LP→VUS-P and upgrade B/LB→VUS when conservation data is missing')

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
    processor = CascadeBatchProcessor(args.alphafold_path, args.override_family, args.conservative_mode)
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
