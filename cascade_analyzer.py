#!/usr/bin/env python3
"""
🌊 CASCADE ANALYZER - Multi-Mechanism Pathogenicity Analysis
The REAL cascade system that coordinates DN, LOF, and GOF analyzers!

Built by Ace & Nova (2025) for revolutionary genetics analysis

Logic:
1. Run DN analysis first (fast, mechanism-aware)
2. IF DN < 0.3 AND frequency < 0.1% THEN try LOF + GOF
3. Return combined results: "DN:0.2(LB) LOF:0.85(LP) GOF:0.1(LB) FINAL:LP"

This is the missing piece that makes the cascade system complete!
"""

import sys
import os
import json
import tempfile
from typing import Dict, List, Optional, Tuple
from pathlib import Path

# Add paths for all analyzers
sys.path.append(str(Path(__file__).parent / "nova_dn"))
sys.path.append(str(Path(__file__).parent.parent / "caller" / "phase1" / "code"))

# Import all three analyzers
from nova_dn.analyzer import NovaDNAnalyzer
from analyzers.lof_analyzer import LOFAnalyzer
from analyzers.gof_variant_analyzer import GOFVariantAnalyzer
from nova_dn.alphafold_sequence import AlphaFoldSequenceExtractor
from nova_dn.sequence_manager import SequenceManager
from biological_router import BiologicalRouter


class CascadeAnalyzer:
    """Coordinates DN, LOF, and GOF analyzers for comprehensive pathogenicity analysis"""
    
    def __init__(self, alphafold_path: str = "/mnt/Arcana/alphafold_human/structures/"):
        self.dn_analyzer = NovaDNAnalyzer(use_smart_filtering=True)
        self.lof_analyzer = LOFAnalyzer(offline_mode=True)
        self.gof_analyzer = GOFVariantAnalyzer(offline_mode=True)
        self.sequence_manager = SequenceManager()
        self.biological_router = BiologicalRouter()  # 🧬 NEW: Smart routing!
        self.temp_files = []
        
        # Gene -> UniProt mappings (expanded for biological routing)
        self.gene_to_uniprot = {
            'TP53': 'P04637', 'COL1A1': 'P02452', 'FGFR3': 'P22607',
            'VWF': 'P04275', 'RYR1': 'P21817', 'FBN1': 'P35555',
            'SCN5A': 'Q14524', 'KCNQ1': 'P51787', 'BRCA1': 'P38398',
            'BRCA2': 'P51587', 'MSH2': 'P43246', 'MLH1': 'P40692',
            'ATM': 'Q13315', 'PTEN': 'P60484',
            # New genes for biological routing
            'CFTR': 'P13569', 'G6PD': 'P11413', 'HEXA': 'P06865',
            'PAH': 'P00439', 'LDLR': 'P01130', 'HBB': 'P68871',
            'RET': 'P07949', 'MET': 'P08581', 'EGFR': 'P00533',
            'COL1A2': 'P08123', 'COL3A1': 'P02461',
            # Tumor suppressors for GO-based routing
            'APC': 'P25054', 'NF1': 'P21359', 'RB1': 'P06400',
            # Additional genes for testing
            'BMPR2': 'Q13873',
            # Hearing loss genes - use canonical isoforms
            'MYO7A': 'Q13402'  # Canonical full-length MYO7A (2215 residues)
        }
    
    def interpret_score(self, score: float) -> str:
        """Convert numeric score to clinical classification"""
        if score >= 0.8:
            return "LP"  # Likely Pathogenic
        elif score >= 0.5:
            return "VUS-P"  # VUS favor pathogenic
        elif score >= 0.3:
            return "VUS"  # Uncertain significance
        else:
            return "LB"  # Likely Benign
    
    def analyze_cascade_biological(self, gene: str, variant: str, gnomad_freq: float = 0.0,
                                  variant_type: str = 'missense', sequence: Optional[str] = None) -> Dict:
        """
        🧬 NEW: Biologically-guided cascade analysis
        Uses BiologicalRouter to determine which analyzers to run

        Args:
            gene: Gene symbol (e.g., 'TP53')
            variant: Variant in p.RefPosAlt format (e.g., 'p.R273H')
            gnomad_freq: Population frequency (0.0-1.0)
            variant_type: Type of variant ('missense', 'frameshift', 'nonsense', etc.)
            sequence: Optional protein sequence (will fetch if not provided)

        Returns:
            Comprehensive biologically-guided analysis results
        """

        # Get UniProt ID - let BiologicalRouter handle the lookup
        # (BiologicalRouter already uses UniversalContext for gene annotation)
        uniprot_id = self.gene_to_uniprot.get(gene)  # Try hardcoded first for speed

        # 🧬 STEP 1: Biological Routing Decision
        routing_result = self.biological_router.route_variant(gene, variant, variant_type, uniprot_id)
        analyzers_to_run = routing_result['analyzers_to_run']

        print(f"🧬 BIOLOGICAL ROUTING for {gene} {variant}")
        print(f"   Strategy: {routing_result['routing_strategy']}")
        print(f"   Confidence: {routing_result['confidence']:.2f}")
        print(f"   Analyzers: {', '.join(analyzers_to_run)}")
        print(f"   Rationale: {routing_result['rationale']}")

        # Get sequence if not provided
        if not sequence:
            try:
                # Extract position for smart sequence selection
                import re
                pos_match = re.search(r'p\.[A-Z](\d+)[A-Z]', variant)
                variant_position = int(pos_match.group(1)) if pos_match else None

                sequence, source, temp_fasta_path = self.sequence_manager.get_best_sequence(
                    gene, uniprot_id, variant_position
                )
                self.temp_files.append(temp_fasta_path)
                print(f"Using {len(sequence)} residue {source} sequence for {gene}")

            except Exception as e:
                return {
                    'error': f'Could not get sequence: {e}',
                    'gene': gene,
                    'variant': variant,
                    'status': 'FAILED'
                }

        results = {
            'gene': gene,
            'variant': variant,
            'variant_type': variant_type,
            'gnomad_freq': gnomad_freq,
            'uniprot_id': uniprot_id,
            'routing_info': routing_result,
            'analyzers_run': analyzers_to_run,
            'scores': {},
            'classifications': {},
            'final_classification': None,
            'final_score': 0.0,
            'explanation': '',
            'status': 'SUCCESS'
        }

        # 🧬 STEP 2: Run Selected Analyzers
        analyzer_results = {}

        if 'DN' in analyzers_to_run:
            print(f"🧬 Running DN analysis...")
            analyzer_results['DN'] = self._run_dn_analysis(gene, variant, sequence, uniprot_id)

        if 'LOF' in analyzers_to_run:
            print(f"🔬 Running LOF analysis...")
            analyzer_results['LOF'] = self._run_lof_analysis(gene, variant, sequence, uniprot_id)

        if 'GOF' in analyzers_to_run:
            print(f"🔥 Running GOF analysis...")
            analyzer_results['GOF'] = self._run_gof_analysis(gene, variant, sequence, uniprot_id)

        # Process results
        for analyzer, result in analyzer_results.items():
            if result['success']:
                results['scores'][analyzer] = result['score']
                results['classifications'][analyzer] = self.interpret_score(result['score'])
                results[f'{analyzer.lower()}_details'] = result['details']
                print(f"   {analyzer} score: {result['score']:.3f} ({results['classifications'][analyzer]})")
            else:
                results['scores'][analyzer] = 0.0
                results['classifications'][analyzer] = 'ERROR'
                print(f"   {analyzer} analysis failed: {result['error']}")

        # 🧬 STEP 3: Final Classification with Primary/Backup Logic
        valid_scores = {k: v for k, v in results['scores'].items() if v > 0}

        if valid_scores:
            # ALWAYS check for synergy first - it takes precedence!
            import math
            primary_analyzer = routing_result.get('primary_analyzer')

            # Get all valid scores with their names
            valid_score_list = [(name, score) for name, score in valid_scores.items() if score > 0]
            valid_score_list.sort(key=lambda x: x[1], reverse=True)  # Sort by score (highest first)

            # Check for mixed mechanism synergy
            synergistic_score = 0
            synergy_explanation = ""

            if len(valid_score_list) >= 2:
                    top_2_scores = [valid_score_list[0][1], valid_score_list[1][1]]
                    top_2_names = [valid_score_list[0][0], valid_score_list[1][0]]

                    # NOVA'S V2 SYNERGY SYSTEM! 🚀
                    synergy_result = self.calculate_synergy_score_v2(
                        {top_2_names[0]: top_2_scores[0], top_2_names[1]: top_2_scores[1]},
                        gene_family=self.get_gene_family(gene)
                    )

                    synergistic_score = synergy_result['synergy_score']
                    synergy_explanation = synergy_result['explanation']

            # Determine final score: use synergy if it's higher than single mechanism
            max_analyzer = max(valid_scores.items(), key=lambda x: x[1])[0]
            max_single_score = valid_scores[max_analyzer]

            if synergistic_score > 0 and synergistic_score > max_single_score:
                # Synergy wins - mixed mechanism is more dangerous
                results['final_score'] = synergistic_score
                results['final_classification'] = self.interpret_score(synergistic_score)
                results['explanation'] = f"Mixed mechanism synergy (biological routing): {synergy_explanation}"
                results['synergy_used'] = True
                results['synergy_details'] = synergy_explanation
            elif primary_analyzer and primary_analyzer in valid_scores and valid_scores[primary_analyzer] == max_single_score:
                # Primary analyzer wins ONLY if it's also the highest score
                results['final_score'] = valid_scores[primary_analyzer]
                results['final_classification'] = results['classifications'][primary_analyzer]
                results['explanation'] = f"Primary {primary_analyzer} score (biological confidence: {routing_result['confidence']:.2f})"

                # Add backup information
                backup_scores = {k: v for k, v in valid_scores.items() if k != primary_analyzer}
                if backup_scores:
                    backup_info = []
                    for analyzer, score in backup_scores.items():
                        classification = results['classifications'][analyzer]
                        backup_info.append(f"{analyzer}:{score:.2f}({classification})")
                    results['backup_scores'] = backup_scores
                    results['backup_summary'] = f"Backup scores: {', '.join(backup_info)}"
                    results['explanation'] += f" | Backup: {', '.join(backup_info)}"
                results['synergy_used'] = False
            else:
                # Highest score wins - don't let "primary" override better evidence
                results['final_score'] = max_single_score
                results['final_classification'] = results['classifications'][max_analyzer]
                results['explanation'] = f"Highest score from {max_analyzer} analyzer (biological routing)"
                results['synergy_used'] = False
        else:
            results['final_score'] = 0.0
            results['final_classification'] = 'ERROR'
            results['explanation'] = "All analyzers failed"

        # Create summary string with primary/backup indication
        primary_analyzer = routing_result.get('primary_analyzer')
        summary_parts = []

        for analyzer in ['DN', 'LOF', 'GOF']:
            if analyzer in results['scores']:
                score = results['scores'][analyzer]
                classification = results['classifications'][analyzer]
                if classification not in ['ERROR']:
                    # Mark primary analyzer
                    if analyzer == primary_analyzer:
                        summary_parts.append(f"*{analyzer}:{score:.2f}({classification})")  # * indicates primary
                    else:
                        summary_parts.append(f"{analyzer}:{score:.2f}({classification})")

        results['summary'] = f"{' '.join(summary_parts)} FINAL:{results['final_classification']}"

        # Add routing confidence to summary for transparency
        if routing_result.get('confidence'):
            results['summary'] += f" [Confidence:{routing_result['confidence']:.2f}]"

        print(f"🎯 BIOLOGICAL RESULT: {results['summary']}")

        # Cleanup
        self.cleanup()

        return results

    def analyze_cascade(self, gene: str, variant: str, gnomad_freq: float = 0.0,
                       sequence: Optional[str] = None) -> Dict:
        """
        Run cascade analysis: DN → (LOF + GOF if needed)
        
        Args:
            gene: Gene symbol (e.g., 'TP53')
            variant: Variant in p.RefPosAlt format (e.g., 'p.R273H')
            gnomad_freq: Population frequency (0.0-1.0)
            sequence: Optional protein sequence (will fetch if not provided)
            
        Returns:
            Comprehensive cascade analysis results
        """
        
        # Get UniProt ID
        uniprot_id = self.gene_to_uniprot.get(gene)
        if not uniprot_id:
            return {
                'error': f'No UniProt mapping for gene {gene}',
                'gene': gene,
                'variant': variant,
                'status': 'FAILED'
            }
        
        # Get sequence if not provided
        if not sequence:
            try:
                # Extract position for smart sequence selection
                import re
                pos_match = re.search(r'p\.[A-Z](\d+)[A-Z]', variant)
                variant_position = int(pos_match.group(1)) if pos_match else None
                
                sequence, source, temp_fasta_path = self.sequence_manager.get_best_sequence(
                    gene, uniprot_id, variant_position
                )
                self.temp_files.append(temp_fasta_path)
                print(f"Using {len(sequence)} residue {source} sequence for {gene}")
                
            except Exception as e:
                return {
                    'error': f'Could not get sequence: {e}',
                    'gene': gene,
                    'variant': variant,
                    'status': 'FAILED'
                }
        
        results = {
            'gene': gene,
            'variant': variant,
            'gnomad_freq': gnomad_freq,
            'uniprot_id': uniprot_id,
            'cascade_triggered': False,
            'analyzers_run': ['DN'],
            'scores': {},
            'classifications': {},
            'final_classification': None,
            'final_score': 0.0,
            'explanation': '',
            'status': 'SUCCESS'
        }
        
        # STEP 1: Run DN Analysis
        print(f"🧬 Running DN analysis for {gene} {variant}")
        try:
            # Load annotations for DN analyzer
            annotations_path = Path(__file__).parent / "resources" / "protein_annotations.json"
            context = None
            if annotations_path.exists():
                from nova_dn.context import load_annotations_json, build_position_context
                try:
                    anns = load_annotations_json(str(annotations_path))
                    pos_match = re.search(r'p\.[A-Z](\d+)[A-Z]', variant)
                    if pos_match:
                        pos = int(pos_match.group(1))
                        context = build_position_context(anns, gene, pos)
                except:
                    pass
            
            dn_result = self.dn_analyzer.analyze(
                sequence, variant, context, gene, uniprot_id
            )
            
            # Extract DN score (use top mechanism score)
            dn_score = dn_result['mechanism_scores'][dn_result['top_mechanism']]
            results['scores']['DN'] = dn_score
            results['classifications']['DN'] = self.interpret_score(dn_score)
            results['dn_details'] = dn_result
            
            print(f"   DN score: {dn_score:.3f} ({results['classifications']['DN']})")
            
        except Exception as e:
            return {
                'error': f'DN analysis failed: {e}',
                'gene': gene,
                'variant': variant,
                'status': 'FAILED'
            }
        
        # STEP 2: Cascade Decision Logic
        cascade_threshold = 0.3
        freq_threshold = 0.001  # 0.1%
        
        if dn_score < cascade_threshold and gnomad_freq < freq_threshold:
            print(f"🌊 CASCADE TRIGGERED: DN={dn_score:.3f} < {cascade_threshold}, freq={gnomad_freq:.4f} < {freq_threshold}")
            results['cascade_triggered'] = True
            results['analyzers_run'].extend(['LOF', 'GOF'])
            
            # STEP 3: Run LOF Analysis
            print(f"🔬 Running LOF analysis...")
            try:
                lof_result = self.lof_analyzer.analyze_lof(
                    variant.replace('p.', ''), sequence, uniprot_id, gene
                )
                lof_score = lof_result.get('lof_score', 0.0)
                results['scores']['LOF'] = lof_score
                results['classifications']['LOF'] = self.interpret_score(lof_score)
                results['lof_details'] = lof_result
                
                print(f"   LOF score: {lof_score:.3f} ({results['classifications']['LOF']})")
                
            except Exception as e:
                print(f"   LOF analysis failed: {e}")
                results['scores']['LOF'] = 0.0
                results['classifications']['LOF'] = 'ERROR'
            
            # STEP 4: Run GOF Analysis
            print(f"🔥 Running GOF analysis...")
            try:
                gof_result = self.gof_analyzer.analyze_gof(
                    variant.replace('p.', ''), sequence, uniprot_id
                )
                gof_score = gof_result.get('gof_score', 0.0)
                results['scores']['GOF'] = gof_score
                results['classifications']['GOF'] = self.interpret_score(gof_score)
                results['gof_details'] = gof_result
                
                print(f"   GOF score: {gof_score:.3f} ({results['classifications']['GOF']})")
                
            except Exception as e:
                print(f"   GOF analysis failed: {e}")
                results['scores']['GOF'] = 0.0
                results['classifications']['GOF'] = 'ERROR'
        
        else:
            print(f"🚫 CASCADE NOT TRIGGERED: DN={dn_score:.3f} >= {cascade_threshold} OR freq={gnomad_freq:.4f} >= {freq_threshold}")
            # Use DN result only
            results['scores']['LOF'] = 0.0
            results['scores']['GOF'] = 0.0
            results['classifications']['LOF'] = 'NOT_RUN'
            results['classifications']['GOF'] = 'NOT_RUN'
        
        # STEP 5: Determine Final Classification
        all_scores = [score for score in results['scores'].values() if score > 0]
        if all_scores:
            results['final_score'] = max(all_scores)
            # Find which analyzer gave the max score
            max_analyzer = max(results['scores'].items(), key=lambda x: x[1])[0]
            results['final_classification'] = results['classifications'][max_analyzer]
            results['explanation'] = f"Highest score from {max_analyzer} analyzer"
        else:
            results['final_score'] = dn_score
            results['final_classification'] = results['classifications']['DN']
            results['explanation'] = "DN analysis only"
        
        # Create summary string
        summary_parts = []
        for analyzer in ['DN', 'LOF', 'GOF']:
            if analyzer in results['scores']:
                score = results['scores'][analyzer]
                classification = results['classifications'][analyzer]
                if classification not in ['NOT_RUN', 'ERROR']:
                    summary_parts.append(f"{analyzer}:{score:.2f}({classification})")
        
        results['summary'] = f"{' '.join(summary_parts)} FINAL:{results['final_classification']}"
        
        print(f"🎯 FINAL RESULT: {results['summary']}")
        
        # Cleanup
        self.cleanup()
        
        return results
    
    def _run_dn_analysis(self, gene: str, variant: str, sequence: str, uniprot_id: str) -> Dict:
        """Run DN analysis and return standardized result"""
        try:
            # Load annotations for DN analyzer
            annotations_path = Path(__file__).parent / "resources" / "protein_annotations.json"
            context = None
            if annotations_path.exists():
                from nova_dn.context import load_annotations_json, build_position_context
                try:
                    anns = load_annotations_json(str(annotations_path))
                    import re
                    pos_match = re.search(r'p\.[A-Z](\d+)[A-Z]', variant)
                    if pos_match:
                        pos = int(pos_match.group(1))
                        context = build_position_context(anns, gene, pos)
                except:
                    pass

            dn_result = self.dn_analyzer.analyze(sequence, variant, context, gene, uniprot_id)
            dn_score = dn_result['mechanism_scores'][dn_result['top_mechanism']]

            return {
                'success': True,
                'score': dn_score,
                'details': dn_result,
                'error': None
            }
        except Exception as e:
            return {
                'success': False,
                'score': 0.0,
                'details': None,
                'error': str(e)
            }

    def _run_lof_analysis(self, gene: str, variant: str, sequence: str, uniprot_id: str) -> Dict:
        """Run LOF analysis and return standardized result"""
        try:
            lof_result = self.lof_analyzer.analyze_lof(
                variant.replace('p.', ''), sequence, uniprot_id, gene
            )
            lof_score = lof_result.get('lof_score', 0.0)

            return {
                'success': True,
                'score': lof_score,
                'details': lof_result,
                'error': None
            }
        except Exception as e:
            return {
                'success': False,
                'score': 0.0,
                'details': None,
                'error': str(e)
            }

    def _run_gof_analysis(self, gene: str, variant: str, sequence: str, uniprot_id: str) -> Dict:
        """Run GOF analysis and return standardized result"""
        try:
            gof_result = self.gof_analyzer.analyze_gof(
                variant.replace('p.', ''), sequence, uniprot_id
            )
            gof_score = gof_result.get('gof_score', 0.0)

            return {
                'success': True,
                'score': gof_score,
                'details': gof_result,
                'error': None
            }
        except Exception as e:
            return {
                'success': False,
                'score': 0.0,
                'details': None,
                'error': str(e)
            }

    def calculate_synergy_score_v2(self, mechanism_scores, gene_family=None):
        """
        🧬 NOVA'S V2 SYNERGY ALGORITHM 🧬
        Calculate synergistic score for mixed-mechanism variants with tiered thresholds,
        contextual weighting, and improved biological rationale
        """
        import math

        # Get top 2 scoring mechanisms
        valid_scores = [(name, score) for name, score in mechanism_scores.items() if score > 0]
        valid_scores.sort(key=lambda x: x[1], reverse=True)  # Highest first

        if len(valid_scores) < 2:
            return {'synergy_score': 0, 'synergy_used': False, 'explanation': 'Need 2+ mechanisms for synergy'}

        top_2_names = [valid_scores[0][0], valid_scores[1][0]]
        top_2_scores = [valid_scores[0][1], valid_scores[1][1]]

        # Tiered thresholds for strength of evidence
        min_score = min(top_2_scores)
        if min_score < 0.3:
            return {'synergy_score': 0, 'synergy_used': False, 'explanation': 'Scores too low for synergy'}
        elif min_score >= 0.7:
            tier = 'strong'
            synergy_boost = 0.3
        elif min_score >= 0.5:
            tier = 'moderate'
            synergy_boost = 0.2
        else:
            tier = 'weak'
            synergy_boost = 0.1

        # Biological plausibility check
        mechanism_pair = tuple(sorted(top_2_names))
        plausible = False
        rationale = ''

        if ('LOF' in mechanism_pair and 'DN' in mechanism_pair):
            plausible = True
            rationale = 'Protein both unstable and poisons complex (classic dominant negative).'
        elif ('DN' in mechanism_pair and 'GOF' in mechanism_pair):
            plausible = True
            rationale = 'Protein is hyperactive AND disrupts wild-type partners.'
        elif ('LOF' in mechanism_pair and 'GOF' in mechanism_pair):
            plausible = False
            rationale = 'Unusual: catalytically dead but also hyperactive. Possible only in special cases (e.g. dimerization sequestering).'

        # Contextual adjustment
        context_multiplier = 1.0
        if gene_family:
            gf = gene_family.lower()
            if gf == 'collagen' and ('DN' in mechanism_pair and 'LOF' in mechanism_pair):
                context_multiplier = 1.1  # boost collagen DN+LOF
            if gf == 'kinase' and ('DN' in mechanism_pair and 'GOF' in mechanism_pair):
                context_multiplier = 0.9  # penalize unless strong

        if not plausible and 'LOF' in mechanism_pair and 'GOF' in mechanism_pair:
            return {
                'synergy_score': min(sum(top_2_scores)/2, 0.5), # keep low but non-zero
                'synergy_used': True,
                'tier': tier,
                'biological_rationale': rationale,
                'explanation': f"LOF+GOF flagged as biologically unlikely; downweighted score assigned for caution."
            }

        # Balance weighting
        balance_factor = min(top_2_scores) / max(top_2_scores)
        synergy_multiplier = 1.0 + (balance_factor * synergy_boost * context_multiplier)

        # Base synergistic score
        base_synergistic_score = math.sqrt(top_2_scores[0]**2 + top_2_scores[1]**2)

        # Final score, normalized to max 1.0
        synergistic_score = min(base_synergistic_score * synergy_multiplier, 1.0)

        return {
            'synergy_score': synergistic_score,
            'synergy_used': True,
            'tier': tier,
            'balance_factor': balance_factor,
            'synergy_multiplier': synergy_multiplier,
            'base_score': base_synergistic_score,
            'biological_rationale': rationale,
            'explanation': f"{tier.capitalize()} mixed mechanism synergy: {top_2_names[0]}+{top_2_names[1]} "
                           f"= sqrt({top_2_scores[0]:.2f}² + {top_2_scores[1]:.2f}²) * {synergy_multiplier:.3f} "
                           f"(balance {balance_factor:.2f}, context {context_multiplier:.2f}) "
                           f"= {synergistic_score:.3f}"
        }

    def get_gene_family(self, gene):
        """Determine gene family for contextual synergy weighting"""
        # Simple gene family classification based on gene name patterns
        gene_upper = gene.upper()

        if gene_upper.startswith('COL') or gene_upper in ['FBN1', 'FBN2', 'TGFBR1', 'TGFBR2']:
            return 'collagen'  # Structural proteins including fibrillin
        elif gene_upper.endswith('K') or 'KINASE' in gene_upper or gene_upper in ['ATM', 'BRCA1', 'BRCA2']:
            return 'kinase'
        elif 'SCN' in gene_upper or 'KCNQ' in gene_upper or 'CACNA' in gene_upper:
            return 'channel'
        else:
            return None

    def cleanup(self):
        """Clean up temporary files"""
        for temp_file in self.temp_files:
            try:
                os.unlink(temp_file)
            except:
                pass
        self.temp_files = []


def main():
    """CLI interface for cascade analysis"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="🌊 CASCADE ANALYZER - Multi-Mechanism Pathogenicity Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze single variant
  python3 cascade_analyzer.py --gene TP53 --variant p.R273H --freq 0.0001
  
  # High frequency variant (cascade won't trigger)
  python3 cascade_analyzer.py --gene TP53 --variant p.P72R --freq 0.25
        """
    )
    
    parser.add_argument('--gene', required=True, help='Gene symbol (e.g., TP53)')
    parser.add_argument('--variant', required=True, help='Variant (e.g., p.R273H)')
    parser.add_argument('--freq', type=float, default=0.0, help='gnomAD frequency (default: 0.0)')
    parser.add_argument('--json', action='store_true', help='Output JSON format')
    
    args = parser.parse_args()
    
    # Run cascade analysis
    analyzer = CascadeAnalyzer()
    result = analyzer.analyze_cascade(args.gene, args.variant, args.freq)
    
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        if result['status'] == 'SUCCESS':
            print(f"\n🎯 CASCADE ANALYSIS COMPLETE")
            print(f"Gene: {result['gene']}")
            print(f"Variant: {result['variant']}")
            print(f"Frequency: {result['gnomad_freq']:.4f}")
            print(f"Cascade Triggered: {result['cascade_triggered']}")
            print(f"Analyzers Run: {', '.join(result['analyzers_run'])}")
            print(f"Summary: {result['summary']}")
            print(f"Final: {result['final_score']:.3f} ({result['final_classification']})")
        else:
            print(f"❌ Analysis failed: {result.get('error', 'Unknown error')}")
    
    return 0 if result['status'] == 'SUCCESS' else 1


if __name__ == "__main__":
    sys.exit(main())
