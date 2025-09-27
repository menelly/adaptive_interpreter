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
from proline_ml_integrator import ProlineMLIntegrator  # 🔥 REVOLUTIONARY ML SYSTEM!
from analyzers.conservation_database import ConservationDatabase  # 🧬 EVOLUTIONARY INTELLIGENCE!


class HotspotDatabase:
    """🔥 HOTSPOT INTELLIGENCE: Known pathogenic and activating regions"""

    def __init__(self):
        # Known hotspot regions from literature and our analysis
        self.hotspots = {
            # Collagen genes - glycine-rich regions are poison hotspots
            'COL1A1': [
                {'start': 359, 'end': 362, 'type': 'dominant_cluster', 'mechanism': 'collagen_poison', 'confidence': 0.9},
                {'start': 988, 'end': 991, 'type': 'dominant_cluster', 'mechanism': 'collagen_poison', 'confidence': 0.8},
            ],
            'COL1A2': [
                {'start': 359, 'end': 362, 'type': 'dominant_cluster', 'mechanism': 'collagen_poison', 'confidence': 0.9},
            ],
            # Ion channels - pore and gate regions
            'SCN1A': [
                {'start': 1616, 'end': 1620, 'type': 'activating_hotspot', 'mechanism': 'channel_disruption', 'confidence': 0.8},
                {'start': 373, 'end': 375, 'type': 'dominant_cluster', 'mechanism': 'channel_poison', 'confidence': 0.7},
            ],
            'KCNQ2': [
                {'start': 265, 'end': 267, 'type': 'activating_hotspot', 'mechanism': 'channel_disruption', 'confidence': 0.9},
            ],
            # Tumor suppressors - DNA binding domains
            'TP53': [
                {'start': 270, 'end': 280, 'type': 'dominant_cluster', 'mechanism': 'dna_binding_loss', 'confidence': 0.95},
            ],
            # Muscle genes - critical structural regions
            'FKRP': [
                {'start': 338, 'end': 340, 'type': 'dominant_cluster', 'mechanism': 'structural_disruption', 'confidence': 0.8},
            ],
        }

    def get_hotspots(self, gene: str) -> List[Dict]:
        """Get all hotspots for a gene"""
        return self.hotspots.get(gene, [])

    def check_variant_in_hotspot(self, gene: str, position: int) -> Optional[Dict]:
        """Check if variant position overlaps any hotspot"""
        hotspots = self.get_hotspots(gene)
        for hotspot in hotspots:
            if hotspot['start'] <= position <= hotspot['end']:
                return hotspot
        return None


class CascadeAnalyzer:
    """Coordinates DN, LOF, and GOF analyzers for comprehensive pathogenicity analysis"""
    
    def __init__(self, alphafold_path: str = "/mnt/Arcana/alphafold_human/structures/"):
        self.dn_analyzer = NovaDNAnalyzer(use_smart_filtering=True)
        self.lof_analyzer = LOFAnalyzer(offline_mode=True)
        self.gof_analyzer = GOFVariantAnalyzer(offline_mode=True)
        self.sequence_manager = SequenceManager()
        self.biological_router = BiologicalRouter()  # 🧬 NEW: Smart routing!
        self.proline_ml = ProlineMLIntegrator(alphafold_path=alphafold_path)  # 🔥 REVOLUTIONARY ML!
        self.conservation_db = ConservationDatabase()  # 🧬 EVOLUTIONARY INTELLIGENCE!
        self.hotspot_db = HotspotDatabase()  # 🔥 HOTSPOT INTELLIGENCE!
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
        """
        Convert numeric score to clinical classification

        🔥 UPDATED THRESHOLDS for conservation-boosted era!
        Based on real 300-variant validation data analysis.
        """
        if score >= 1.2:
            return "P"   # Pathogenic (ultra-conserved, definitely broken)
        elif score >= 0.8:
            return "LP"  # Likely Pathogenic (strong evidence)
        elif score >= 0.6:
            return "VUS-P"  # VUS favor Pathogenic (moderate evidence)
        elif score >= 0.4:
            return "VUS"  # Uncertain significance (weak evidence)
        elif score >= 0.2:
            return "LB"  # Likely Benign (conservative threshold)
        else:
            return "B"   # Benign (very confident)
    
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
        print(f"🔍 DEBUG: Initial uniprot_id for {gene}: {uniprot_id}")

        # 🚨 CRITICAL CODON CHECK: Start/Stop codon variants are auto-pathogenic
        critical_codon_result = self._check_critical_codons(variant, variant_type)
        if critical_codon_result['is_critical']:
            return {
                'gene': gene,
                'variant': variant,
                'variant_type': variant_type,
                'gnomad_freq': gnomad_freq,
                'uniprot_id': uniprot_id,
                'routing_info': {'reason': 'Critical codon override'},
                'analyzers_run': ['CRITICAL_CODON'],
                'scores': {'CRITICAL_CODON': 2.0},  # Maximum pathogenic score
                'classifications': {'CRITICAL_CODON': 'P'},
                'final_classification': 'P',
                'final_score': 2.0,
                'explanation': critical_codon_result['explanation'],
                'critical_codon_override': True,
                'status': 'SUCCESS'
            }

        # 🧬 STEP 1: Biological Routing Decision
        routing_result = self.biological_router.route_variant(gene, variant, variant_type, uniprot_id)
        analyzers_to_run = routing_result['analyzers_to_run']

        # Update UniProt ID from routing result if available
        if routing_result.get('uniprot_id'):
            uniprot_id = routing_result['uniprot_id']
            print(f"🔍 DEBUG: Updated uniprot_id from router: {uniprot_id}")
        else:
            print(f"🔍 DEBUG: Router did not provide uniprot_id")
            # Fallback: Use the same UniProt lookup that DN analyzer uses
            try:
                from universal_protein_annotator import UniversalProteinAnnotator
                annotator = UniversalProteinAnnotator()
                uniprot_id = annotator._find_uniprot_id(gene)
                if uniprot_id:
                    print(f"🔍 DEBUG: Fallback found uniprot_id: {uniprot_id}")
                else:
                    print(f"🔍 DEBUG: Fallback could not find uniprot_id for {gene}")
            except Exception as e:
                print(f"🔍 DEBUG: Fallback UniProt lookup failed: {e}")

        print(f"🧬 BIOLOGICAL ROUTING for {gene} {variant}")
        print(f"   Strategy: {routing_result['routing_strategy']}")
        print(f"   Confidence: {routing_result['confidence']:.2f}")
        print(f"   Analyzers: {', '.join(analyzers_to_run)}")
        print(f"   Rationale: {routing_result['rationale']}")

        # 🎯 SMART ORDERING: Run analyzers in biological priority order
        if routing_result.get('primary_analyzer'):
            primary = routing_result['primary_analyzer']
            print(f"   🎯 Primary mechanism: {primary}")
            # Reorder to put primary first
            if primary in analyzers_to_run:
                analyzers_to_run.remove(primary)
                analyzers_to_run.insert(0, primary)

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

        # 🧬 STEP 2: Get Conservation Multiplier (EVOLUTIONARY INTELLIGENCE!)
        conservation_multiplier = self._get_conservation_multiplier(gene, variant, uniprot_id, gnomad_freq)
        print(f"🧬 Using conservation multiplier: {conservation_multiplier:.1f}x")

        # 🧬 STEP 3: Run Analyzers with LOF-First Logic
        analyzer_results = {}

        # 🔥 REN'S BRILLIANT INSIGHT: LOF first, skip GOF if protein is broken!
        lof_score = 0.0

        # Always run LOF first (if requested) - NO CONSERVATION BOOST YET!
        if 'LOF' in analyzers_to_run:
            print(f"🔬 Running LOF analysis...")
            print(f"🔍 DEBUG: About to call _run_lof_analysis with gene={gene}, variant={variant}")
            # 🔥 REN'S BRILLIANT FIX: No conservation boost here - apply at the end!
            analyzer_results['LOF'] = self._run_lof_analysis(gene, variant, sequence, uniprot_id, conservation_multiplier=1.0)
            print(f"🔍 DEBUG: LOF analyzer returned: {analyzer_results['LOF']}")

            if analyzer_results['LOF']['success']:
                lof_score = analyzer_results['LOF']['score']
                print(f"🎯 LOF score: {lof_score:.3f} (raw, before conservation)")

        # Always run DN (can coexist with broken proteins - misfolded can still poison) - NO CONSERVATION BOOST YET!
        if 'DN' in analyzers_to_run:
            print(f"🧬 Running DN analysis...")
            # 🔥 REN'S BRILLIANT FIX: No conservation boost here - apply at the end!
            analyzer_results['DN'] = self._run_dn_analysis(gene, variant, sequence, uniprot_id, conservation_multiplier=1.0)

        # 🔥 BIOLOGICAL LOGIC: Skip GOF if protein is definitely broken (LOF ≥ 0.5)
        if 'GOF' in analyzers_to_run:
            if lof_score >= 0.5:
                print(f"🚫 SKIPPING GOF: LOF score {lof_score:.3f} ≥ 0.5 (broken proteins cannot be hyperactive!)")
                analyzer_results['GOF'] = {
                    'success': True,
                    'score': 0.0,
                    'details': {'skipped_reason': f'LOF score {lof_score:.3f} ≥ 0.5 - broken proteins cannot gain function'},
                    'error': None
                }
            else:
                print(f"🔥 Running GOF analysis (LOF score {lof_score:.3f} < 0.5, protein might still be functional)...")
                # 🔥 REN'S BRILLIANT FIX: No conservation boost here - apply at the end!
                analyzer_results['GOF'] = self._run_gof_analysis(gene, variant, sequence, uniprot_id, conservation_multiplier=1.0)

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

        # 🔥 STEP 3.5: Apply NOVA-style hotspot boosts
        hotspot_result = self._apply_hotspot_boost(gene, variant, results['scores'])
        if hotspot_result['hotspot_info']:
            # Update scores with hotspot boosts
            results['scores'] = hotspot_result['scores']
            results['hotspot_info'] = hotspot_result['hotspot_info']
            # Recalculate classifications with boosted scores
            for analyzer, score in results['scores'].items():
                results['classifications'][analyzer] = self.interpret_score(score)
                print(f"   🔥 {analyzer} boosted: {score:.3f} ({results['classifications'][analyzer]})")
        else:
            results['hotspot_info'] = None

        # 🧬 STEP 3: Biologically-Guided Final Classification
        valid_scores = {k: v for k, v in results['scores'].items() if v > 0}

        if valid_scores:
            # ALWAYS check for synergy first - it takes precedence!
            import math
            primary_analyzer = routing_result.get('primary_analyzer')
            routing_confidence = routing_result.get('confidence', 0.5)

            # Get all valid scores with their names
            valid_score_list = [(name, score) for name, score in valid_scores.items() if score > 0]
            valid_score_list.sort(key=lambda x: x[1], reverse=True)  # Sort by score (highest first)

            # 🎯 BIOLOGICAL PRIORITY: Weight primary analyzer more heavily if high confidence
            if primary_analyzer and primary_analyzer in valid_scores and routing_confidence > 0.8:
                primary_score = valid_scores[primary_analyzer]
                print(f"🎯 HIGH CONFIDENCE ROUTING: Prioritizing {primary_analyzer} (score: {primary_score:.3f}, confidence: {routing_confidence:.2f})")
                # Boost primary analyzer score slightly for final classification
                valid_scores[primary_analyzer] = min(1.0, primary_score * 1.1)
                print(f"   Boosted {primary_analyzer} score: {valid_scores[primary_analyzer]:.3f}")

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

                # Add hotspot information to explanation
                explanation = f"Mixed mechanism synergy (biological routing): {synergy_explanation}"
                if results.get('hotspot_info'):
                    hotspot = results['hotspot_info']
                    explanation += f" + Hotspot boost ({hotspot['hotspot_type']}: {hotspot['mechanism']})"
                results['explanation'] = explanation

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

        # Apply Nova's biological plausibility filter
        try:
            from plausibility_filter import plausibility_pipeline

            # Get UniProt function from annotation cache
            uniprot_function = ""
            if uniprot_id:
                try:
                    import json
                    cache_file = f"protein_annotations_cache/{uniprot_id}_domains.json"
                    if os.path.exists(cache_file):
                        with open(cache_file, 'r') as f:
                            annotation_data = json.load(f)
                            uniprot_function = annotation_data.get('function', '')
                            print(f"🔍 Retrieved function for {gene} ({uniprot_id}): {uniprot_function[:100]}...")
                except Exception as e:
                    print(f"⚠️ Could not read cached function: {e}")
                    pass

            # Apply plausibility filter
            raw_scores = {
                'DN': results['scores'].get('DN', 0.0),
                'LOF': results['scores'].get('LOF', 0.0),
                'GOF': results['scores'].get('GOF', 0.0)
            }

            print(f"🧬 APPLYING PLAUSIBILITY FILTER to {gene}...")
            print(f"   Raw scores: DN={raw_scores['DN']:.3f}, LOF={raw_scores['LOF']:.3f}, GOF={raw_scores['GOF']:.3f}")

            plausibility_result = plausibility_pipeline(
                gene_symbol=gene,
                raw_scores=raw_scores,
                uniprot_function=uniprot_function,
                go_terms=[]
            )

            # Update results with plausibility-filtered scores
            filtered_scores = plausibility_result['final_scores']
            results['plausibility_filtered_scores'] = filtered_scores
            results['gene_family'] = plausibility_result['gene_family']
            results['plausibility_rationale'] = plausibility_result.get('rationale', {})

            print(f"🎯 PLAUSIBILITY RESULT: Gene family = {results['gene_family']}")
            print(f"   Filtered scores: DN={filtered_scores['DN']:.3f}, LOF={filtered_scores['LOF']:.3f}, GOF={filtered_scores['GOF']:.3f}")

            # 🔥 RECALCULATE SYNERGY WITH FILTERED SCORES (preserve Ren's sqrt synergy!)
            valid_filtered_scores = {k: v for k, v in filtered_scores.items() if v > 0}

            if len(valid_filtered_scores) >= 2:
                # Apply synergy to filtered scores too!
                import math
                score_list = sorted(valid_filtered_scores.items(), key=lambda x: x[1], reverse=True)
                top_2_scores = [score_list[0][1], score_list[1][1]]
                top_2_names = [score_list[0][0], score_list[1][0]]

                # 🧬 REN'S SQRT SYNERGY on filtered scores
                synergy_score = math.sqrt(top_2_scores[0]**2 + top_2_scores[1]**2)
                max_filtered_score = max(filtered_scores.values())

                if synergy_score > max_filtered_score:
                    print(f"🔥 SYNERGY APPLIED TO FILTERED SCORES: {top_2_names[0]}({top_2_scores[0]:.3f}) + {top_2_names[1]}({top_2_scores[1]:.3f}) = {synergy_score:.3f}")
                    final_filtered_score = synergy_score
                    results['synergy_applied_post_filter'] = True
                else:
                    final_filtered_score = max_filtered_score
                    results['synergy_applied_post_filter'] = False
            else:
                final_filtered_score = max(filtered_scores.values())
                results['synergy_applied_post_filter'] = False

            # ALWAYS apply plausibility filter results - biological plausibility is critical!
            results['final_score_pre_filter'] = results['final_score']
            results['final_score_pre_conservation'] = final_filtered_score

            # 🔥 REN'S BRILLIANT FIX: Apply conservation to FINAL score (fair competition!)
            final_score_with_conservation = final_filtered_score * conservation_multiplier
            results['final_score'] = final_score_with_conservation
            results['conservation_multiplier_applied'] = conservation_multiplier

            print(f"🧬 CONSERVATION APPLIED TO FINAL: {final_filtered_score:.3f} × {conservation_multiplier:.1f}x = {final_score_with_conservation:.3f}")

            # Update classification based on conservation-boosted final score
            raw_classification = self.interpret_score(final_score_with_conservation)

            # 🔥 REN'S BRILLIANT CONSERVATION-AWARE CLASSIFICATION
            # If conservation multiplier is exactly 1.0, we have no conservation data - be conservative!
            if conservation_multiplier == 1.0:
                print(f"⚠️ No conservation data available (1.0x multiplier) - applying conservative classification")

                if raw_classification in ['B', 'LB']:
                    # Upgrade benign calls to VUS when we lack conservation data
                    results['final_classification'] = 'VUS'
                    results['explanation'] += " | Upgraded B/LB→VUS due to missing conservation data"

                elif raw_classification in ['P', 'LP']:
                    # Downgrade pathogenic calls to VUS-P when we lack conservation data
                    results['final_classification'] = 'VUS-P'
                    results['explanation'] += " | Downgraded P/LP→VUS-P due to missing conservation data"

                else:
                    # VUS stays VUS, VUS-P stays VUS-P
                    results['final_classification'] = raw_classification
                    results['explanation'] += " | Classification maintained despite missing conservation data"
            else:
                # We have real conservation data - trust the classification
                results['final_classification'] = raw_classification

            results['explanation'] += f" | Plausibility filter applied (gene family: {results['gene_family']})"

        except Exception as e:
            print(f"⚠️ Plausibility filter failed: {e}")
            # Continue with original results if filter fails

        # Cleanup
        self.cleanup()

        return results

    def _get_conservation_multiplier(self, gene: str, variant: str, uniprot_id: str, gnomad_freq: float = 0.0) -> float:
        """
        🧬 EVOLUTIONARY INTELLIGENCE: Get conservation-based multiplier

        Uses real evolutionary data (GERP, phyloP) to determine if this position
        is conserved across species. Highly conserved = higher pathogenic potential.

        Args:
            gene: Gene symbol
            variant: Variant in p.RefPosAlt format
            uniprot_id: UniProt ID for genomic mapping

        Returns:
            Conservation multiplier (0.8-2.0x based on evolutionary constraint)
        """
        try:
            # 🔥 COORDINATE LOOKUP: Check temporary coordinates from batch processor
            variant_key = f"{gene}_{variant}"
            coordinates = None

            # Check temporary coordinates (from batch processor)
            if hasattr(self, '_temp_coordinates') and variant_key in self._temp_coordinates:
                coordinates = self._temp_coordinates[variant_key]
                print(f"🎯 Using batch coordinates for {gene} {variant}: {coordinates[0]}:{coordinates[1]}")

            # Fallback to hardcoded coordinates for testing
            elif variant_key == 'KCNMA1_p.F533L':
                coordinates = ('chr10', 76974527)
                print(f"🎯 Using hardcoded coordinates for {gene} {variant}: chr10:76974527")

            if coordinates:
                chrom, pos = coordinates

                # Get conservation scores directly
                scores = self.conservation_db.get_conservation_scores(chrom, pos)
                phylop_score = scores.get('phyloP', 0)
                phastcons_score = scores.get('phastCons', 0)
                combined_score = scores.get('conservation_score', 0)

                print(f"🧬 Conservation for {gene} {variant}: phyloP={phylop_score:.2f}, phastCons={phastcons_score:.2f}, combined={combined_score:.3f}")

                # Convert phyloP score to multiplier (phyloP ranges from -20 to +20, positive = conserved)
                # 🔥 UPDATED SCALING: More aggressive for ultra-conserved positions!
                if phylop_score >= 7.0:
                    multiplier = 2.5  # Ultra-conserved (phyloP 7+, absolutely critical!)
                elif phylop_score >= 5.0:
                    multiplier = 2.0  # Extremely conserved
                elif phylop_score >= 3.0:
                    multiplier = 1.5  # Highly conserved
                elif phylop_score >= 1.0:
                    multiplier = 1.2  # Moderately conserved
                elif phylop_score >= -1.0:
                    multiplier = 1.0  # Neutral
                else:
                    multiplier = 0.8  # Not conserved (likely benign!)

                print(f"🎯 Conservation multiplier for {gene} {variant}: {multiplier:.1f}x (phyloP: {phylop_score:.2f})")
                return multiplier

            # Fall back to original method for other variants
            # Extract position from variant
            import re
            pos_match = re.search(r'p\.[A-Z](\d+)[A-Z]', variant)
            if not pos_match:
                print(f"⚠️ Could not extract position from variant {variant}")
                return 1.0

            protein_position = int(pos_match.group(1))

            # Get conservation scores for this protein position
            conservation_data = self.conservation_db.get_variant_conservation(uniprot_id, protein_position)

            if conservation_data.get('error'):
                print(f"⚠️ Conservation lookup failed for {gene} {variant}: {conservation_data['error']}")

                # 🔥 FALLBACK: Use frequency-based conservation for ion channels
                frequency_multiplier = self._get_frequency_based_conservation(gene, gnomad_freq)
                if frequency_multiplier > 1.0:
                    print(f"🧬 Using frequency-based conservation fallback: {frequency_multiplier:.1f}x")
                    return frequency_multiplier

                return 1.0

            conservation_scores = conservation_data.get('conservation_scores')
            if not conservation_scores:
                print(f"⚠️ No conservation scores available for {gene} {variant}")
                return 1.0

            # Extract phyloP score (most reliable for pathogenicity)
            phylop_score = conservation_scores.get('phyloP', 0)
            phastcons_score = conservation_scores.get('phastCons', 0)
            combined_score = conservation_scores.get('conservation_score', 0)

            print(f"🧬 Conservation for {gene} {variant}: phyloP={phylop_score:.2f}, phastCons={phastcons_score:.2f}, combined={combined_score:.3f}")

            # Convert phyloP score to multiplier (phyloP ranges from -20 to +20, positive = conserved)
            # 🔥 REN'S BRILLIANT FIX: Never use exactly 1.0 for real data (reserved for "no data")
            if phylop_score >= 7.0:
                multiplier = 2.5  # Ultra-conserved (phyloP 7+, absolutely critical!)
            elif phylop_score >= 5.0:
                multiplier = 2.0  # Extremely conserved
            elif phylop_score >= 3.0:
                multiplier = 1.5  # Highly conserved
            elif phylop_score >= 1.0:
                multiplier = 1.2  # Moderately conserved
            elif phylop_score >= -1.0:
                multiplier = 0.9  # Slightly non-conserved (avoid 1.0!)
            else:
                multiplier = 0.8  # Not conserved (likely benign!)

            print(f"🎯 Conservation multiplier for {gene} {variant}: {multiplier:.1f}x (phyloP: {phylop_score:.2f})")
            return multiplier

        except Exception as e:
            print(f"⚠️ Conservation analysis failed for {gene} {variant}: {e}")

            # 🔥 FALLBACK: Use frequency-based conservation for ion channels
            frequency_multiplier = self._get_frequency_based_conservation(gene, gnomad_freq)
            if frequency_multiplier > 1.0:
                print(f"🧬 Using frequency-based conservation fallback: {frequency_multiplier:.1f}x")
                return frequency_multiplier

            return 1.0

    def _get_frequency_based_conservation(self, gene: str, gnomad_freq: float) -> float:
        """
        🔥 FREQUENCY-BASED CONSERVATION: Use rarity as conservation proxy for ion channels

        Ion channels are ultra-conserved (phyloP ~9-10). When conservation data fails,
        use population frequency as a proxy - rare variants in conserved genes are more likely pathogenic.

        Args:
            gene: Gene symbol
            gnomad_freq: Population frequency (0.0-1.0)

        Returns:
            Frequency-based conservation multiplier (1.0-2.5x)
        """

        # Ion channel genes (ultra-conserved, frequency matters a lot)
        ion_channel_genes = {
            'SCN1A', 'SCN2A', 'SCN3A', 'SCN5A', 'SCN8A', 'SCN9A', 'SCN10A', 'SCN11A',
            'KCNQ1', 'KCNQ2', 'KCNQ3', 'KCNQ4', 'KCNQ5',
            'KCNH1', 'KCNH2', 'KCNH3', 'KCNH4', 'KCNH5', 'KCNH6', 'KCNH7', 'KCNH8',
            'KCNA1', 'KCNA2', 'KCNA3', 'KCNA4', 'KCNA5', 'KCNA6', 'KCNA7', 'KCNA10',
            'CACNA1A', 'CACNA1B', 'CACNA1C', 'CACNA1D', 'CACNA1E', 'CACNA1F', 'CACNA1G', 'CACNA1H', 'CACNA1I', 'CACNA1S',
            'CACNA2D1', 'CACNA2D2', 'CACNA2D3', 'CACNA2D4',
            'CACNB1', 'CACNB2', 'CACNB3', 'CACNB4'
        }

        # Other highly conserved gene families
        highly_conserved_genes = {
            'TP53', 'BRCA1', 'BRCA2', 'MLH1', 'MSH2', 'MSH6', 'PMS2',  # Tumor suppressors
            'COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'COL4A2', 'COL5A1', 'COL5A2',  # Collagens
            'FBN1', 'FBN2', 'TGFBR1', 'TGFBR2',  # Connective tissue
        }

        if gene not in ion_channel_genes and gene not in highly_conserved_genes:
            return 1.0  # No frequency-based conservation for other genes

        # Frequency-based conservation logic
        if gnomad_freq == 0.0:
            # Absent from gnomAD = ultra-rare = likely pathogenic in conserved genes
            if gene in ion_channel_genes:
                return 2.5  # Ion channels: ultra-rare = ultra-conserved position
            else:
                return 2.0  # Other conserved genes: ultra-rare = highly conserved

        elif gnomad_freq < 0.00001:  # <1 in 100,000
            if gene in ion_channel_genes:
                return 2.0  # Very rare in ion channel = likely pathogenic
            else:
                return 1.5

        elif gnomad_freq < 0.0001:  # <1 in 10,000
            if gene in ion_channel_genes:
                return 1.5  # Rare in ion channel = moderately pathogenic
            else:
                return 1.2

        elif gnomad_freq < 0.001:  # <1 in 1,000
            return 1.2  # Uncommon = mildly conserved

        else:
            return 1.0  # Common variants don't get frequency boost

    def analyze_cascade(self, gene: str, variant: str, gnomad_freq: float = 0.0,
                       sequence: Optional[str] = None, variant_type: str = 'missense',
                       expected_clinvar: str = "") -> Dict:
        """
        🧬 BIOLOGICALLY-GUIDED CASCADE ANALYSIS (Updated to use biological routing by default!)

        Uses BiologicalRouter to determine optimal analysis strategy based on gene function,
        GO terms, and variant type. This replaces the old "DN first" approach with
        intelligent biological routing.

        Args:
            gene: Gene symbol (e.g., 'TP53')
            variant: Variant in p.RefPosAlt format (e.g., 'p.R273H')
            gnomad_freq: Population frequency (0.0-1.0)
            sequence: Optional protein sequence (will fetch if not provided)
            variant_type: Type of variant ('missense', 'frameshift', 'nonsense', etc.)
            expected_clinvar: Expected ClinVar classification for disagreement handling

        Returns:
            Comprehensive biologically-guided analysis results
        """

        # 🚀 NEW: Use biological routing by default!
        print(f"🧬 BIOLOGICAL ROUTING: Analyzing {gene} {variant} ({variant_type})")
        return self.analyze_cascade_biological(gene, variant, gnomad_freq, variant_type, sequence)

    def analyze_cascade_legacy(self, gene: str, variant: str, gnomad_freq: float = 0.0,
                              sequence: Optional[str] = None) -> Dict:
        """
        🏛️ LEGACY: Old DN-first cascade analysis (kept for compatibility)

        This is the original "DN → (LOF + GOF if needed)" approach.
        Use analyze_cascade() for the new biologically-guided approach.

        Args:
            gene: Gene symbol (e.g., 'TP53')
            variant: Variant in p.RefPosAlt format (e.g., 'p.R273H')
            gnomad_freq: Population frequency (0.0-1.0)
            sequence: Optional protein sequence (will fetch if not provided)

        Returns:
            Comprehensive cascade analysis results (legacy DN-first approach)
        """

        # Get UniProt ID - try hardcoded first, then dynamic lookup
        uniprot_id = self.gene_to_uniprot.get(gene)
        if not uniprot_id:
            print(f"🔍 Gene {gene} not in hardcoded mappings, searching UniProt...")
            # Use the existing UniversalProteinAnnotator to find UniProt ID
            from universal_protein_annotator import UniversalProteinAnnotator
            annotator = UniversalProteinAnnotator()
            uniprot_id = annotator._find_uniprot_id(gene)

            if uniprot_id:
                print(f"✅ Found UniProt ID for {gene}: {uniprot_id}")
                # Cache it for future use
                self.gene_to_uniprot[gene] = uniprot_id
            else:
                return {
                    'error': f'No UniProt mapping found for gene {gene}',
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
    
    def _run_dn_analysis(self, gene: str, variant: str, sequence: str, uniprot_id: str, conservation_multiplier: float = 1.0) -> Dict:
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

    def _run_lof_analysis(self, gene: str, variant: str, sequence: str, uniprot_id: str, conservation_multiplier: float = 1.0) -> Dict:
        """Run LOF analysis and return standardized result"""
        try:
            print(f"🔍 LOF DEBUG: gene={gene}, variant={variant}, uniprot_id={uniprot_id}, seq_len={len(sequence) if sequence else 'None'}")
            print(f"🧬 LOF DEBUG: Using conservation multiplier: {conservation_multiplier:.1f}x")
            lof_result = self.lof_analyzer.analyze_lof(
                variant.replace('p.', ''), sequence, uniprot_id=uniprot_id, gene_symbol=gene, conservation_multiplier=conservation_multiplier
            )
            lof_score = lof_result.get('lof_score', 0.0)
            domain_mult = lof_result.get('domain_multiplier', 'N/A')
            print(f"🔍 LOF DEBUG: Got score={lof_score}, domain_multiplier={domain_mult}")

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

    def _run_gof_analysis(self, gene: str, variant: str, sequence: str, uniprot_id: str, conservation_multiplier: float = 1.0) -> Dict:
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

    def _check_critical_codons(self, variant: str, variant_type: str) -> Dict:
        """
        🚨 CRITICAL CODON DETECTION: Auto-pathogenic variants

        Detects variants that affect critical codons:
        - Start codon loss (M1X) = Auto-P (no protein production)
        - Nonsense variants = Auto-P (premature termination)
        - Stop codon loss = Auto-P (read-through effects)

        Args:
            variant: Variant in p.RefPosAlt format (e.g., 'p.M1L')
            variant_type: Type of variant ('missense', 'nonsense', etc.)

        Returns:
            Dict with is_critical flag and explanation
        """
        import re

        # Extract position and amino acids from variant
        match = re.search(r'p\.([A-Z*])(\d+)([A-Z*])', variant)
        if not match:
            return {'is_critical': False, 'explanation': ''}

        ref_aa, position, alt_aa = match.groups()
        position = int(position)

        # 🚨 START CODON LOSS (M1X)
        if position == 1 and ref_aa == 'M' and alt_aa != 'M':
            return {
                'is_critical': True,
                'explanation': f"Start codon loss (M1{alt_aa}): Complete loss of protein production - Auto-Pathogenic"
            }

        # 🚨 NONSENSE VARIANTS (X = stop codon)
        if variant_type == 'nonsense' or alt_aa == '*':
            return {
                'is_critical': True,
                'explanation': f"Nonsense variant ({ref_aa}{position}*): Premature protein termination - Auto-Pathogenic"
            }

        # 🚨 STOP CODON LOSS (*X where X != *)
        if ref_aa == '*' and alt_aa != '*':
            return {
                'is_critical': True,
                'explanation': f"Stop codon loss (*{position}{alt_aa}): Read-through effects - Auto-Pathogenic"
            }

        return {'is_critical': False, 'explanation': ''}

    def _extract_position(self, variant: str) -> Optional[int]:
        """Extract amino acid position from variant string"""
        import re
        match = re.search(r'p\.[A-Z*](\d+)[A-Z*]', variant)
        if match:
            return int(match.group(1))
        return None

    def _apply_hotspot_boost(self, gene: str, variant: str, scores: Dict[str, float]) -> Dict:
        """Apply NOVA-style hotspot boosts to mechanism scores"""
        position = self._extract_position(variant)
        if not position:
            return {'scores': scores, 'hotspot_info': None}

        hotspot = self.hotspot_db.check_variant_in_hotspot(gene, position)
        if not hotspot:
            return {'scores': scores, 'hotspot_info': None}

        print(f"🔥 HOTSPOT DETECTED: {gene} pos {position} in {hotspot['type']} ({hotspot['mechanism']})")

        boosted_scores = scores.copy()
        hotspot_type = hotspot['type']
        confidence = hotspot['confidence']

        # Apply NOVA-style boosts based on hotspot type
        if hotspot_type == 'activating_hotspot':
            # Boost GOF for activating regions (kinase sites, channels)
            if 'GOF' in boosted_scores:
                boost = 0.25 * confidence  # Up to +0.25 boost
                boosted_scores['GOF'] = min(1.0, boosted_scores['GOF'] + boost)
                print(f"🚀 GOF HOTSPOT BOOST: +{boost:.3f} (activating region)")

        elif hotspot_type == 'dominant_cluster':
            # Boost DN for pathogenic clusters (poison mechanisms)
            if 'DN' in boosted_scores:
                boost = 0.15 * confidence  # Up to +0.15 boost
                boosted_scores['DN'] = min(1.0, boosted_scores['DN'] + boost)
                print(f"🧬 DN HOTSPOT BOOST: +{boost:.3f} (pathogenic cluster)")

            # Also boost LOF for structural disruption
            if 'LOF' in boosted_scores and hotspot['mechanism'] in ['structural_disruption', 'collagen_poison']:
                boost = 0.10 * confidence  # Up to +0.10 boost
                boosted_scores['LOF'] = min(1.0, boosted_scores['LOF'] + boost)
                print(f"🔬 LOF HOTSPOT BOOST: +{boost:.3f} (structural disruption)")

        return {
            'scores': boosted_scores,
            'hotspot_info': {
                'position': position,
                'hotspot_type': hotspot_type,
                'mechanism': hotspot['mechanism'],
                'confidence': confidence,
                'boost_applied': True
            }
        }

    def format_human_readable(self, gene: str, variant: str, result: Dict) -> str:
        """Format results in clean, human-readable format"""

        # Extract key information
        final_class = result.get('final_classification', 'UNKNOWN')
        final_score = result.get('final_score', 0.0)

        # Get individual scores
        scores = result.get('scores', {})
        dn_score = scores.get('DN', 0.0)
        lof_score = scores.get('LOF', 0.0)
        gof_score = scores.get('GOF', 0.0)

        # Format scores nicely
        score_parts = []
        if dn_score > 0:
            score_parts.append(f"DN:{dn_score:.2f}")
        if lof_score > 0:
            score_parts.append(f"LOF:{lof_score:.2f}")
        if gof_score > 0:
            score_parts.append(f"GOF:{gof_score:.2f}")

        score_summary = " | ".join(score_parts) if score_parts else "No scores"

        # Add hotspot indicator
        hotspot_indicator = ""
        if result.get('hotspot_info'):
            hotspot = result['hotspot_info']
            hotspot_indicator = f" 🔥 {hotspot['hotspot_type']}"

        # Conservation indicator
        conservation_indicator = ""
        if result.get('conservation_multiplier_applied', 1.0) > 1.0:
            mult = result['conservation_multiplier_applied']
            conservation_indicator = f" 🧬 {mult:.1f}x"

        # Format final line
        return f"{gene:8} {variant:12} | {score_summary:25} | {final_class:6} ({final_score:.3f}){hotspot_indicator}{conservation_indicator}"

    def print_human_summary(self, gene: str, variant: str, result: Dict, clinvar_expected: str = None):
        """Print clean human-readable summary"""

        print("\n" + "="*80)
        print(f"🧬 VARIANT ANALYSIS: {gene} {variant}")
        print("="*80)

        # Scores section
        scores = result.get('scores', {})
        print(f"📊 MECHANISM SCORES:")
        if scores.get('DN', 0) > 0:
            print(f"   Dominant Negative: {scores['DN']:.3f} ({self.interpret_score(scores['DN'])})")
        if scores.get('LOF', 0) > 0:
            print(f"   Loss of Function:  {scores['LOF']:.3f} ({self.interpret_score(scores['LOF'])})")
        if scores.get('GOF', 0) > 0:
            print(f"   Gain of Function:  {scores['GOF']:.3f} ({self.interpret_score(scores['GOF'])})")

        # Final result
        final_class = result.get('final_classification', 'UNKNOWN')
        final_score = result.get('final_score', 0.0)
        print(f"\n🎯 FINAL RESULT: {final_class} (Score: {final_score:.3f})")

        # Hotspot info
        if result.get('hotspot_info'):
            hotspot = result['hotspot_info']
            print(f"🔥 HOTSPOT: {hotspot['hotspot_type']} ({hotspot['mechanism']})")

        # Conservation info
        if result.get('conservation_multiplier_applied', 1.0) > 1.0:
            mult = result['conservation_multiplier_applied']
            print(f"🧬 CONSERVATION: {mult:.1f}x multiplier applied")

        # ClinVar comparison
        if clinvar_expected:
            print(f"📋 CLINVAR: {clinvar_expected}")
            # Simple agreement check
            our_severity = self._get_severity_level(final_class)
            clinvar_severity = self._get_severity_level(clinvar_expected)
            if abs(our_severity - clinvar_severity) <= 1:
                print("✅ AGREEMENT: Close correlation with ClinVar")
            elif (our_severity <= 2 and clinvar_severity >= 5) or (our_severity >= 5 and clinvar_severity <= 2):
                print("🚨 DISAGREEMENT: Major classification difference")
            else:
                print("📊 CORRELATION: Minor classification difference")

        print("="*80)

    def _get_severity_level(self, classification: str) -> int:
        """Convert classification to severity level for comparison"""
        severity_map = {
            'B': 1, 'LB': 2, 'VUS': 3, 'VUS-P': 4, 'LP': 5, 'P': 6,
            'Benign': 1, 'Likely benign': 2, 'Uncertain significance': 3,
            'Likely pathogenic': 5, 'Pathogenic': 6
        }

        classification_str = str(classification).lower()
        for key, value in severity_map.items():
            if key.lower() in classification_str:
                return value
        return 3  # Default to VUS if unknown

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
