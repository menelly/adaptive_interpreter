#!/usr/bin/env python3
"""
🔥 GAIN OF FUNCTION VARIANT ANALYZER
Built by Ace for revolutionary GOF variant analysis

This analyzer determines if a specific variant causes gain-of-function
by analyzing structural and biochemical impacts that lead to:
- Constitutive activation
- Loss of autoinhibition  
- Increased binding affinity
- Degradation resistance
- Enhanced dimerization

NO HARDCODED GENES - Pure mathematical analysis!
"""

from typing import Dict, Any, List, Tuple
import re
import math
import logging
import sys
import os
from .smart_protein_analyzer import SmartProteinAnalyzer
from .conservation_database import ConservationDatabase
from AdaptiveInterpreter.utils.proline_ml_integrator import ProlineMLIntegrator
from AdaptiveInterpreter.utils.gly_cys_ml_integrator import FamilyAwareGlyCysIntegrator as GlyCysMLIntegrator
# from AdaptiveInterpreter.utils.active_site_scanner import ActiveSiteScanner # TODO: Restore this once the file is found
# from AdaptiveInterpreter.utils.motif_detector import MotifDetector # TODO: Restore this once the file is found

# 🎯 Add domain awareness system
try:
    from AdaptiveInterpreter.data_processing.universal_protein_annotator import UniversalProteinAnnotator
    DOMAIN_AWARENESS_AVAILABLE = True
except ImportError:
    DOMAIN_AWARENESS_AVAILABLE = False

logger = logging.getLogger(__name__)

class GOFVariantAnalyzer:
    """Analyze gain of function potential for specific variants"""
    
    def __init__(self, offline_mode=False):
        self.name = "GOFVariantAnalyzer"
        self.offline_mode = offline_mode
        self.smart_analyzer = SmartProteinAnalyzer(offline_mode=offline_mode)
        self.conservation_db = ConservationDatabase() if not offline_mode else None

        # 🎯 Initialize domain awareness system
        self.protein_annotator = UniversalProteinAnnotator() if DOMAIN_AWARENESS_AVAILABLE and not offline_mode else None

        if offline_mode:
            logger.info("🔥 GOF Analyzer initialized in OFFLINE MODE - skipping conservation analysis")
        
        # Grantham distance matrix - CRITICAL for all mechanisms!
        self.grantham_matrix = {
            ('A', 'R'): 112, ('A', 'N'): 111, ('A', 'D'): 126, ('A', 'C'): 195,
            ('A', 'Q'): 91, ('A', 'E'): 107, ('A', 'G'): 60, ('A', 'H'): 86,
            ('A', 'I'): 94, ('A', 'L'): 96, ('A', 'K'): 106, ('A', 'M'): 84,
            ('A', 'F'): 113, ('A', 'P'): 27, ('A', 'S'): 99, ('A', 'T'): 58,
            ('A', 'W'): 148, ('A', 'Y'): 112, ('A', 'V'): 64,
            # Add more as needed - this is a subset for now
        }
        
        # Amino acid properties for GOF analysis
        self.aa_properties = {
            'G': {'size': 1, 'charge': 0, 'hydrophobic': False, 'flexibility': 'high', 'stability': 'low'},
            'A': {'size': 2, 'charge': 0, 'hydrophobic': True, 'flexibility': 'medium', 'stability': 'medium'},
            'V': {'size': 3, 'charge': 0, 'hydrophobic': True, 'flexibility': 'low', 'stability': 'high'},
            'L': {'size': 4, 'charge': 0, 'hydrophobic': True, 'flexibility': 'low', 'stability': 'high'},
            'I': {'size': 4, 'charge': 0, 'hydrophobic': True, 'flexibility': 'low', 'stability': 'high'},
            'M': {'size': 4, 'charge': 0, 'hydrophobic': True, 'flexibility': 'medium', 'stability': 'medium'},
            'F': {'size': 5, 'charge': 0, 'hydrophobic': True, 'flexibility': 'low', 'stability': 'high'},
            'W': {'size': 6, 'charge': 0, 'hydrophobic': True, 'flexibility': 'low', 'stability': 'high'},
            'P': {'size': 3, 'charge': 0, 'hydrophobic': False, 'flexibility': 'rigid', 'stability': 'high'},
            'S': {'size': 2, 'charge': 0, 'hydrophobic': False, 'flexibility': 'high', 'stability': 'low'},
            'T': {'size': 3, 'charge': 0, 'hydrophobic': False, 'flexibility': 'medium', 'stability': 'medium'},
            'C': {'size': 2, 'charge': 0, 'hydrophobic': False, 'flexibility': 'medium', 'stability': 'medium'},
            'Y': {'size': 5, 'charge': 0, 'hydrophobic': False, 'flexibility': 'medium', 'stability': 'high'},
            'N': {'size': 3, 'charge': 0, 'hydrophobic': False, 'flexibility': 'high', 'stability': 'low'},
            'Q': {'size': 4, 'charge': 0, 'hydrophobic': False, 'flexibility': 'high', 'stability': 'low'},
            'D': {'size': 3, 'charge': -1, 'hydrophobic': False, 'flexibility': 'high', 'stability': 'low'},
            'E': {'size': 4, 'charge': -1, 'hydrophobic': False, 'flexibility': 'high', 'stability': 'low'},
            'K': {'size': 4, 'charge': 1, 'hydrophobic': False, 'flexibility': 'high', 'stability': 'low'},
            'R': {'size': 5, 'charge': 1, 'hydrophobic': False, 'flexibility': 'high', 'stability': 'low'},
            'H': {'size': 4, 'charge': 0.5, 'hydrophobic': False, 'flexibility': 'high', 'stability': 'medium'}
        }
        
        # GOF mechanism signatures
        self.gof_signatures = {
            'constitutive_activation': {
                'charge_disruption_weight': 0.3,
                'flexibility_increase_weight': 0.2,
                'hydrophobic_disruption_weight': 0.25,
                'size_change_weight': 0.25
            },
            'increased_binding_affinity': {
                'charge_enhancement_weight': 0.4,
                'hydrophobic_enhancement_weight': 0.3,
                'size_optimization_weight': 0.3
            },
            'degradation_resistance': {
                'stability_increase_weight': 0.4,
                'flexibility_decrease_weight': 0.3,
                'hydrophobic_increase_weight': 0.3
            },
            'autoinhibition_loss': {
                'flexibility_increase_weight': 0.35,
                'charge_disruption_weight': 0.35,
                'size_disruption_weight': 0.3
            }
        }
    
    def get_grantham_distance(self, aa1: str, aa2: str) -> float:
        """Get Grantham distance between two amino acids"""
        if aa1 == aa2:
            return 0.0
        
        # Try both directions in the matrix
        key = (aa1, aa2)
        reverse_key = (aa2, aa1)
        
        if key in self.grantham_matrix:
            return self.grantham_matrix[key]
        elif reverse_key in self.grantham_matrix:
            return self.grantham_matrix[reverse_key]
        else:
            # Fallback calculation if not in matrix
            return self._calculate_grantham_fallback(aa1, aa2)
    
    def _calculate_grantham_fallback(self, aa1: str, aa2: str) -> float:
        """Fallback Grantham calculation for missing pairs"""
        if aa1 not in self.aa_properties or aa2 not in self.aa_properties:
            return 100.0  # Default moderate distance
        
        prop1 = self.aa_properties[aa1]
        prop2 = self.aa_properties[aa2]
        
        # Simplified Grantham-like calculation
        size_diff = abs(prop1['size'] - prop2['size']) * 20
        charge_diff = abs(prop1['charge'] - prop2['charge']) * 50
        
        return min(size_diff + charge_diff, 215)  # Cap at max Grantham distance
    
    def analyze_gof(self, mutation: str, sequence: str, uniprot_id: str = None, **kwargs) -> Dict[str, Any]:
        """Main GOF analysis with domain awareness wrapper"""
        # Get the gene symbol from kwargs for domain awareness
        gene_symbol = kwargs.get('gene_symbol', '')

        # Run the core GOF analysis
        result = self._analyze_gof_core(mutation, sequence, uniprot_id, **kwargs)

        # Apply domain awareness to the final result
        if 'gof_score' in result and 'error' not in result:
            try:
                # Parse position from mutation for domain awareness
                match = re.match(r'([A-Z])(\d+)([A-Z])', mutation)
                if match and uniprot_id:
                    position = int(match.group(2))
                    original_score = result['gof_score']
                    domain_aware_score, domain_multiplier = self._apply_domain_awareness_to_gof_score(
                        original_score, position, uniprot_id, gene_symbol
                    )

                    # 🔥 REN'S BRILLIANT FIX: ADD ML PROLINE PANIC TO GOF TOO!
                    ml_proline_multiplier = 1.0
                    ref_aa = match.group(1)
                    alt_aa = match.group(3)

                    if ref_aa == 'P' or alt_aa == 'P':  # Proline substitution detected!
                        try:
                            # Import and use our revolutionary ML system
                            from AdaptiveInterpreter.utils.proline_ml_integrator import get_ml_proline_multiplier
                            variant_str = f"p.{mutation}"
                            ml_proline_multiplier = get_ml_proline_multiplier(gene_symbol, variant_str)
                            logger.info(f"🔥 GOF ML PROLINE: {gene_symbol} {variant_str} -> ML multiplier = {ml_proline_multiplier:.3f}")

                        except Exception as e:
                            logger.warning(f"⚠️ GOF ML proline system error: {e}")
                            # Fallback to conservative multiplier for GOF
                            if ref_aa == 'P':  # Proline loss - less likely to cause GOF
                                ml_proline_multiplier = 0.8  # Reduce GOF score
                                logger.info(f"🔥 GOF PROLINE FALLBACK: {gene_symbol} {variant_str} -> fallback multiplier = {ml_proline_multiplier:.3f}")

                    # Apply both domain awareness AND ML proline panic
                    final_gof_score = domain_aware_score * ml_proline_multiplier

                    # Update the result with all multipliers
                    result['gof_score'] = min(final_gof_score, 1.0)
                    result['domain_multiplier'] = domain_multiplier
                    result['ml_proline_multiplier'] = ml_proline_multiplier
                    result['original_gof_score'] = original_score

                    # Update prediction based on new score
                    result['prediction'] = ('GOF_LIKELY' if result['gof_score'] > 0.6 else
                                          'GOF_POSSIBLE' if result['gof_score'] > 0.3 else 'GOF_UNLIKELY')

            except Exception as e:
                logger.warning(f"⚠️ Domain awareness application failed: {e}")
                result['domain_multiplier'] = 1.0
                result['ml_proline_multiplier'] = 1.0

        return result

    def _analyze_gof_core(self, mutation: str, sequence: str, uniprot_id: str = None, **kwargs) -> Dict[str, Any]:
        """
        TRIPLE-GATED GOF ANALYSIS - Revolutionary efficiency!

        GATE 1: Grantham screening (is this change worth analyzing?)
        GATE 2: Math mechanism screening (does math suggest GOF potential?)
        GATE 3: Structural modeling (full FASTA + AlphaFold analysis)

        Args:
            mutation: Mutation string (e.g., "R175H")
            sequence: Protein sequence
            uniprot_id: UniProt ID for additional analysis
            **kwargs: Additional parameters

        Returns:
            Dict with GOF analysis results
        """
        try:
            # Parse mutation
            match = re.match(r'([A-Z])(\d+)([A-Z])', mutation)
            if not match:
                return {'error': f'Invalid mutation format: {mutation}', 'gof_score': 0.0}

            original_aa, position_str, mutant_aa = match.groups()
            position = int(position_str)

            # Validate position
            if position < 1 or position > len(sequence):
                return {'error': f'Position {position} out of range for sequence length {len(sequence)}', 'gof_score': 0.0}

            # 🎯 NOVA'S EARLY MOTIF DETECTION - Check for canonical GOF variants first!
            # try:
            #     from AdaptiveInterpreter.utils.motif_detector import detect_regulatory_context
            #
            #     motif_context = detect_regulatory_context(sequence, original_aa, mutant_aa, position)
            #     motif_score = motif_context.get('score', 0.0)
            #
            #     # If Nova's motif detection gives a high score (0.8+), this is likely a canonical GOF
            #     # Return immediately with high confidence
            #     if motif_score >= 0.8:
            #         logger.info(f"🎯 CANONICAL GOF DETECTED: {motif_context.get('explanation', 'Unknown')}")
            #         return {
            #             'gof_score': motif_score,
            #             'prediction': 'GOF_LIKELY',
            #             'confidence': 0.95,
            #             'analysis_level': 'NOVA_MOTIF_DETECTION',
            #             'motif_context': motif_context,
            #             'reason': f"Canonical GOF variant detected by motif analysis: {motif_context.get('explanation', 'Unknown')}"
            #         }
            #
            # except Exception as e:
            #     logger.warning(f"⚠️ Early motif detection failed: {e}")
            #     # Continue with normal analysis

            # Check if original amino acid matches sequence
            sequence_mismatch = False
            if sequence[position - 1] != original_aa:
                sequence_mismatch = True
                logger.warning(f"⚠️ Sequence mismatch: Expected {original_aa} at position {position}, found {sequence[position - 1]}")
                logger.warning(f"⚠️ Likely transcript/isoform difference - proceeding with mathematical analysis only")

            # REVOLUTIONARY APPROACH: SKIP GRANTHAM ENTIRELY FOR GOF!
            # Small changes can cause MASSIVE regulatory disruption!
            grantham_distance = self.get_grantham_distance(original_aa, mutant_aa)
            logger.info(f"🎯 Grantham distance: {grantham_distance:.1f} (but we don't care for GOF!)")

            # STRUCTURAL CONFIDENCE ANALYSIS - Using AlphaFold pLDDT scores!
            conservation_data = None
            conservation_multiplier = 1.0

            if uniprot_id:
                try:
                    # 🛡️ BULLETPROOF CONFIDENCE ANALYSIS with smart fallbacks!
                    confidence_result = self._get_structural_confidence_with_fallbacks(uniprot_id, position)

                    conservation_multiplier = confidence_result['multiplier']
                    conservation_data = confidence_result

                    logger.info(f"🧬 {confidence_result['source'].upper()}: {confidence_result['score']:.1f} → {conservation_multiplier:.2f}x GOF multiplier")
                    logger.info(f"   Rationale: {confidence_result['rationale']}")

                except Exception as e:
                    logger.warning(f"⚠️ All confidence analysis failed: {e}")
                    # Ultimate fallback - always works!
                    conservation_multiplier = 1.0
                    conservation_data = {
                        'score': 50.0,
                        'source': 'emergency_fallback',
                        'confidence': 'none',
                        'multiplier': 1.0,
                        'rationale': 'Emergency fallback - no confidence data available'
                    }
                    logger.info(f"🆘 EMERGENCY FALLBACK: Using 1.0x multiplier")
            else:
                logger.debug(f"🔥 No UniProt ID provided - using base GOF analysis")
                conservation_multiplier = 1.0
                conservation_data = {
                    'score': 50.0,
                    'source': 'no_uniprot',
                    'confidence': 'none',
                    'multiplier': 1.0,
                    'rationale': 'No UniProt ID provided'
                }

            # GATE 1: REGULATORY DISRUPTION SCREENING - Enhanced with conservation!
            regulatory_disruption_score = self._calculate_regulatory_disruption_potential(
                original_aa, mutant_aa, position, sequence
            )

            # Apply conservation enhancement to regulatory disruption!
            conservation_enhanced_disruption = regulatory_disruption_score * conservation_multiplier
            conservation_enhanced_disruption = min(conservation_enhanced_disruption, 1.0)  # Cap at 1.0

            logger.info(f"🎯 Base regulatory disruption: {regulatory_disruption_score:.3f}")
            logger.info(f"🎯 Conservation-enhanced disruption: {conservation_enhanced_disruption:.3f}")

            # Use conservation-enhanced score for filtering
            if conservation_enhanced_disruption < 0.1:
                # No regulatory disruption potential - truly boring change
                return {
                    'mutation': mutation,
                    'grantham_distance': grantham_distance,
                    'regulatory_disruption_score': regulatory_disruption_score,
                    'conservation_enhanced_disruption': conservation_enhanced_disruption,
                    'conservation_multiplier': conservation_multiplier,
                    'conservation_data': conservation_data,
                    'gof_score': 0.0,
                    'prediction': 'GOF_UNLIKELY',
                    'confidence': 0.9,
                    'analysis_level': 'GATE_1_CONSERVATION_REGULATORY_FILTERED',
                    'reason': f'No conservation-enhanced regulatory disruption ({conservation_enhanced_disruption:.3f})'
                }

            # GATE 2: CONSERVATION-ENHANCED REGULATORY MECHANISM ANALYSIS!
            regulatory_gof_scores = self._run_regulatory_gof_analysis(original_aa, mutant_aa, position, sequence)
            regulatory_overall_score = self._calculate_regulatory_gof_score(regulatory_gof_scores, conservation_enhanced_disruption)

            logger.info(f"🎯 Regulatory GOF score: {regulatory_overall_score:.3f}")

            if regulatory_overall_score < 0.2:
                # Low regulatory GOF potential - return regulatory results
                analysis_level = 'GATE_2_REGULATORY_SCREENING'
                if sequence_mismatch:
                    analysis_level += '_SEQUENCE_MISMATCH'

                return {
                    'mutation': mutation,
                    'grantham_distance': grantham_distance,
                    'regulatory_disruption_score': regulatory_disruption_score,
                    'conservation_enhanced_disruption': conservation_enhanced_disruption,
                    'conservation_multiplier': conservation_multiplier,
                    'conservation_data': conservation_data,
                    'gof_mechanisms': regulatory_gof_scores,
                    'gof_score': regulatory_overall_score,
                    'prediction': 'GOF_UNLIKELY',
                    'confidence': 0.8,
                    'analysis_level': analysis_level,
                    'reason': f'Conservation-enhanced regulatory score ({regulatory_overall_score:.3f}) below GOF threshold',
                    'sequence_mismatch': sequence_mismatch
                }

            # GATE 3: ENHANCED REGULATORY ANALYSIS - Skip structural if sequence mismatch!
            if sequence_mismatch:
                # Sequence mismatch - return enhanced regulatory results
                logger.info(f"🔬 Sequence mismatch detected - using enhanced regulatory analysis")
                enhanced_result = self._run_enhanced_regulatory_analysis(
                    original_aa, mutant_aa, position, conservation_enhanced_disruption,
                    regulatory_gof_scores, conservation_multiplier, conservation_data
                )
                return enhanced_result
            else:
                # Normal structural + regulatory analysis
                structural_result = self._run_structural_regulatory_analysis(
                    original_aa, mutant_aa, position, sequence, uniprot_id,
                    conservation_enhanced_disruption, regulatory_gof_scores,
                    conservation_multiplier, conservation_data
                )
                return structural_result

        except Exception as e:
            return {'error': f'GOF analysis failed: {str(e)}', 'gof_score': 0.0}

    def _analyze_constitutive_activation(self, original_aa: str, mutant_aa: str, position: int,
                                       sequence: str, grantham_distance: float) -> float:
        """Analyze potential for constitutive activation"""
        if original_aa not in self.aa_properties or mutant_aa not in self.aa_properties:
            return 0.0

        orig_props = self.aa_properties[original_aa]
        mut_props = self.aa_properties[mutant_aa]

        score = 0.0

        # Charge disruption (breaks regulatory interactions)
        charge_change = abs(mut_props['charge'] - orig_props['charge'])
        if charge_change > 0:
            score += charge_change * self.gof_signatures['constitutive_activation']['charge_disruption_weight']

        # Flexibility increase (disrupts autoinhibitory conformations)
        flexibility_map = {'low': 1, 'medium': 2, 'high': 3, 'rigid': 0}
        orig_flex = flexibility_map.get(orig_props['flexibility'], 2)
        mut_flex = flexibility_map.get(mut_props['flexibility'], 2)

        if mut_flex > orig_flex:
            flexibility_increase = (mut_flex - orig_flex) / 3.0
            score += flexibility_increase * self.gof_signatures['constitutive_activation']['flexibility_increase_weight']

        # Hydrophobic disruption (breaks hydrophobic regulatory patches)
        if orig_props['hydrophobic'] and not mut_props['hydrophobic']:
            score += self.gof_signatures['constitutive_activation']['hydrophobic_disruption_weight']

        # Size change impact
        size_change = abs(mut_props['size'] - orig_props['size'])
        if size_change > 1:
            score += (size_change / 5.0) * self.gof_signatures['constitutive_activation']['size_change_weight']

        # Grantham distance amplification
        grantham_factor = min(grantham_distance / 100.0, 1.5)  # Cap at 1.5x

        return min(score * grantham_factor, 1.0)

    def _analyze_binding_affinity(self, original_aa: str, mutant_aa: str, position: int,
                                sequence: str, grantham_distance: float) -> float:
        """Analyze potential for increased binding affinity"""
        if original_aa not in self.aa_properties or mutant_aa not in self.aa_properties:
            return 0.0

        orig_props = self.aa_properties[original_aa]
        mut_props = self.aa_properties[mutant_aa]

        score = 0.0

        # Charge enhancement (stronger ionic interactions)
        charge_enhancement = abs(mut_props['charge']) - abs(orig_props['charge'])
        if charge_enhancement > 0:
            score += charge_enhancement * self.gof_signatures['increased_binding_affinity']['charge_enhancement_weight']

        # Hydrophobic enhancement (stronger hydrophobic interactions)
        if not orig_props['hydrophobic'] and mut_props['hydrophobic']:
            score += self.gof_signatures['increased_binding_affinity']['hydrophobic_enhancement_weight']

        # Size optimization (better fit in binding pockets)
        size_change = mut_props['size'] - orig_props['size']
        if 1 <= size_change <= 2:  # Optimal size increase
            score += self.gof_signatures['increased_binding_affinity']['size_optimization_weight']

        # Grantham distance amplification
        grantham_factor = min(grantham_distance / 120.0, 1.3)  # Different scaling for binding

        return min(score * grantham_factor, 1.0)

    def _analyze_degradation_resistance(self, original_aa: str, mutant_aa: str, position: int,
                                      sequence: str, grantham_distance: float) -> float:
        """Analyze potential for degradation resistance"""
        if original_aa not in self.aa_properties or mutant_aa not in self.aa_properties:
            return 0.0

        orig_props = self.aa_properties[original_aa]
        mut_props = self.aa_properties[mutant_aa]

        score = 0.0

        # Stability increase (harder to degrade)
        stability_map = {'low': 1, 'medium': 2, 'high': 3}
        orig_stab = stability_map.get(orig_props['stability'], 2)
        mut_stab = stability_map.get(mut_props['stability'], 2)

        if mut_stab > orig_stab:
            stability_increase = (mut_stab - orig_stab) / 2.0
            score += stability_increase * self.gof_signatures['degradation_resistance']['stability_increase_weight']

        # Flexibility decrease (more rigid, harder to unfold)
        flexibility_map = {'low': 1, 'medium': 2, 'high': 3, 'rigid': 0}
        orig_flex = flexibility_map.get(orig_props['flexibility'], 2)
        mut_flex = flexibility_map.get(mut_props['flexibility'], 2)

        if mut_flex < orig_flex and mut_flex > 0:
            flexibility_decrease = (orig_flex - mut_flex) / 3.0
            score += flexibility_decrease * self.gof_signatures['degradation_resistance']['flexibility_decrease_weight']

        # Hydrophobic increase (more stable core)
        if not orig_props['hydrophobic'] and mut_props['hydrophobic']:
            score += self.gof_signatures['degradation_resistance']['hydrophobic_increase_weight']

        # Grantham distance amplification
        grantham_factor = min(grantham_distance / 80.0, 1.4)  # Different scaling for stability

        return min(score * grantham_factor, 1.0)

    def _analyze_autoinhibition_loss(self, original_aa: str, mutant_aa: str, position: int,
                                   sequence: str, grantham_distance: float) -> float:
        """Analyze potential for autoinhibition loss"""
        if original_aa not in self.aa_properties or mutant_aa not in self.aa_properties:
            return 0.0

        orig_props = self.aa_properties[original_aa]
        mut_props = self.aa_properties[mutant_aa]

        score = 0.0

        # Flexibility increase (disrupts autoinhibitory conformations)
        flexibility_map = {'low': 1, 'medium': 2, 'high': 3, 'rigid': 0}
        orig_flex = flexibility_map.get(orig_props['flexibility'], 2)
        mut_flex = flexibility_map.get(mut_props['flexibility'], 2)

        if mut_flex > orig_flex:
            flexibility_increase = (mut_flex - orig_flex) / 3.0
            score += flexibility_increase * self.gof_signatures['autoinhibition_loss']['flexibility_increase_weight']

        # Charge disruption (breaks autoinhibitory salt bridges)
        charge_change = abs(mut_props['charge'] - orig_props['charge'])
        if charge_change > 0:
            score += charge_change * self.gof_signatures['autoinhibition_loss']['charge_disruption_weight']

        # Size disruption (breaks autoinhibitory packing)
        size_change = abs(mut_props['size'] - orig_props['size'])
        if size_change > 1:
            score += (size_change / 5.0) * self.gof_signatures['autoinhibition_loss']['size_disruption_weight']

        # Grantham distance amplification
        grantham_factor = min(grantham_distance / 90.0, 1.4)

        return min(score * grantham_factor, 1.0)

    def _calculate_overall_gof_score(self, gof_scores: Dict[str, float], grantham_distance: float) -> float:
        """Calculate overall GOF score from mechanism scores"""
        # Weight the mechanisms based on their general importance
        mechanism_weights = {
            'constitutive_activation': 0.3,
            'increased_binding_affinity': 0.25,
            'degradation_resistance': 0.2,
            'autoinhibition_loss': 0.25
        }

        weighted_score = 0.0
        for mechanism, score in gof_scores.items():
            if mechanism in mechanism_weights:
                weighted_score += score * mechanism_weights[mechanism]

        # Apply Grantham distance scaling - higher distances more likely to cause GOF
        grantham_scaling = min(grantham_distance / 150.0, 1.2)  # Cap at 1.2x boost

        final_score = weighted_score * grantham_scaling

        return min(final_score, 1.0)  # Cap at 1.0

    def _run_math_gof_screening(self, original_aa: str, mutant_aa: str, position: int,
                               sequence: str, grantham_distance: float) -> Dict[str, float]:
        """
        GATE 2: REVOLUTIONARY REGULATORY CONTEXT-AWARE GOF SCREENING!

        Combines traditional mechanism analysis with regulatory context disruption
        """
        # Traditional mechanism scores
        traditional_scores = {}

        # 1. Constitutive Activation Analysis
        traditional_scores['constitutive_activation'] = self._analyze_constitutive_activation(
            original_aa, mutant_aa, position, sequence, grantham_distance
        )

        # 2. Increased Binding Affinity Analysis
        traditional_scores['increased_binding_affinity'] = self._analyze_binding_affinity(
            original_aa, mutant_aa, position, sequence, grantham_distance
        )

        # 3. Degradation Resistance Analysis
        traditional_scores['degradation_resistance'] = self._analyze_degradation_resistance(
            original_aa, mutant_aa, position, sequence, grantham_distance
        )

        # 4. Autoinhibition Loss Analysis
        traditional_scores['autoinhibition_loss'] = self._analyze_autoinhibition_loss(
            original_aa, mutant_aa, position, sequence, grantham_distance
        )

        # REVOLUTIONARY ADDITION: Regulatory Context Analysis!
        context_scores = self._analyze_regulatory_context_disruption(
            original_aa, mutant_aa, position, sequence
        )

        # Integrate context scores with traditional scores
        enhanced_scores = {}

        for mechanism in traditional_scores:
            base_score = traditional_scores[mechanism]

            # Apply regulatory context enhancement
            context_multiplier = 1.0

            # Phosphorylation disruption enhances all mechanisms
            if context_scores['phosphorylation_disruption'] > 0.5:
                context_multiplier *= 1.5  # Major boost for phospho site loss!
                logger.info(f"🎯 Phospho disruption boosting {mechanism}: {base_score:.3f} -> {base_score * context_multiplier:.3f}")

            # Flexibility disruption especially enhances constitutive activation
            if mechanism == 'constitutive_activation' and context_scores['flexibility_regulatory_disruption'] > 0.5:
                context_multiplier *= 1.3
                logger.info(f"🎯 Flexibility disruption boosting constitutive activation")

            # Charge disruption enhances binding affinity and autoinhibition loss
            if mechanism in ['increased_binding_affinity', 'autoinhibition_loss'] and context_scores['charge_regulatory_disruption'] > 0.3:
                context_multiplier *= 1.2
                logger.info(f"🎯 Charge disruption boosting {mechanism}")

            # Apply context enhancement
            enhanced_scores[mechanism] = min(base_score * context_multiplier, 1.0)

        # Add context scores as separate mechanisms for transparency
        enhanced_scores.update({
            f"context_{key}": value for key, value in context_scores.items()
        })

        return enhanced_scores

    def _run_structural_gof_analysis(self, original_aa: str, mutant_aa: str, position: int,
                                   sequence: str, uniprot_id: str, grantham_distance: float,
                                   math_scores: Dict[str, float]) -> Dict[str, Any]:
        """
        GATE 3: Full structural GOF analysis with FASTA + AlphaFold
        Only runs on high-potential variants that passed Gates 1 & 2
        """
        try:
            # Start with math scores as baseline
            structural_scores = math_scores.copy()

            # TODO: Add structural analysis here
            # This is where we would:
            # 1. Create FASTA mutation pairs
            # 2. Analyze AlphaFold structures
            # 3. Check regulatory domain context
            # 4. Refine scores based on structural impact

            # For now, enhance math scores with structural context hints
            structural_enhancement = self._estimate_structural_context(position, sequence, uniprot_id)

            # Apply structural enhancement to math scores
            for mechanism in structural_scores:
                structural_scores[mechanism] *= structural_enhancement
                structural_scores[mechanism] = min(structural_scores[mechanism], 1.0)  # Cap at 1.0

            # Calculate enhanced overall score
            overall_score = self._calculate_overall_gof_score(structural_scores, grantham_distance)

            return {
                'mutation': f"{original_aa}{position}{mutant_aa}",
                'grantham_distance': grantham_distance,
                'gof_mechanisms': structural_scores,
                'gof_score': overall_score,
                'prediction': 'GOF_LIKELY' if overall_score > 0.6 else 'GOF_POSSIBLE' if overall_score > 0.3 else 'GOF_UNLIKELY',
                'confidence': min(grantham_distance / 215.0, 1.0),
                'analysis_level': 'GATE_3_STRUCTURAL_MODELING',
                'structural_enhancement': structural_enhancement,
                'math_scores': math_scores  # Include original math scores for comparison
            }

        except Exception as e:
            # Fallback to math scores if structural analysis fails
            overall_score = self._calculate_overall_gof_score(math_scores, grantham_distance)
            return {
                'mutation': f"{original_aa}{position}{mutant_aa}",
                'grantham_distance': grantham_distance,
                'gof_mechanisms': math_scores,
                'gof_score': overall_score,
                'prediction': 'GOF_LIKELY' if overall_score > 0.6 else 'GOF_POSSIBLE' if overall_score > 0.3 else 'GOF_UNLIKELY',
                'confidence': min(grantham_distance / 215.0, 1.0),
                'analysis_level': 'GATE_3_STRUCTURAL_FALLBACK',
                'error': f'Structural analysis failed: {str(e)}'
            }

    def _estimate_structural_context(self, position: int, sequence: str, uniprot_id: str) -> float:
        """
        Estimate structural context for GOF enhancement
        TODO: Replace with real AlphaFold analysis
        """
        # For now, return a modest enhancement factor
        # In the future, this would analyze:
        # - Is position in a regulatory domain?
        # - Is position near allosteric sites?
        # - Is position in autoinhibitory regions?
        # - AlphaFold confidence at this position

        # Simple heuristic for now - positions in middle of protein more likely regulatory
        protein_length = len(sequence)
        relative_position = position / protein_length

        if 0.2 < relative_position < 0.8:
            return 1.2  # Modest enhancement for middle regions
        else:
            return 1.0  # No enhancement for terminal regions

    def _run_enhanced_math_analysis(self, original_aa: str, mutant_aa: str, position: int,
                                   grantham_distance: float, math_scores: Dict[str, float]) -> Dict[str, Any]:
        """
        Enhanced mathematical analysis for sequence mismatch cases
        Uses more sophisticated mathematical modeling when structural analysis isn't possible
        """
        try:
            # Start with base math scores
            enhanced_scores = math_scores.copy()

            # Apply enhanced mathematical modeling for sequence mismatch cases
            # This compensates for not having structural analysis

            # 1. Boost scores for high Grantham distances (likely more impactful)
            grantham_boost = min(grantham_distance / 150.0, 1.3)  # Up to 30% boost

            # 2. Apply amino acid change type analysis
            change_type_multiplier = self._analyze_change_type_impact(original_aa, mutant_aa)

            # 3. Enhance each mechanism score
            for mechanism in enhanced_scores:
                enhanced_scores[mechanism] *= grantham_boost * change_type_multiplier
                enhanced_scores[mechanism] = min(enhanced_scores[mechanism], 1.0)  # Cap at 1.0

            # Calculate enhanced overall score
            overall_score = self._calculate_overall_gof_score(enhanced_scores, grantham_distance)

            return {
                'mutation': f"{original_aa}{position}{mutant_aa}",
                'grantham_distance': grantham_distance,
                'gof_mechanisms': enhanced_scores,
                'gof_score': overall_score,
                'prediction': 'GOF_LIKELY' if overall_score > 0.6 else 'GOF_POSSIBLE' if overall_score > 0.3 else 'GOF_UNLIKELY',
                'confidence': min(grantham_distance / 215.0, 1.0),
                'analysis_level': 'GATE_3_ENHANCED_MATH_SEQUENCE_MISMATCH',
                'sequence_mismatch': True,
                'grantham_boost': grantham_boost,
                'change_type_multiplier': change_type_multiplier,
                'math_scores': math_scores  # Include original math scores for comparison
            }

        except Exception as e:
            # Fallback to basic math scores
            overall_score = self._calculate_overall_gof_score(math_scores, grantham_distance)
            return {
                'mutation': f"{original_aa}{position}{mutant_aa}",
                'grantham_distance': grantham_distance,
                'gof_mechanisms': math_scores,
                'gof_score': overall_score,
                'prediction': 'GOF_LIKELY' if overall_score > 0.6 else 'GOF_POSSIBLE' if overall_score > 0.3 else 'GOF_UNLIKELY',
                'confidence': min(grantham_distance / 215.0, 1.0),
                'analysis_level': 'GATE_3_ENHANCED_MATH_FALLBACK',
                'sequence_mismatch': True,
                'error': f'Enhanced math analysis failed: {str(e)}'
            }

    def _analyze_change_type_impact(self, original_aa: str, mutant_aa: str) -> float:
        """
        Analyze the type of amino acid change for GOF impact
        Returns a multiplier based on change characteristics
        """
        if original_aa not in self.aa_properties or mutant_aa not in self.aa_properties:
            return 1.0

        orig_props = self.aa_properties[original_aa]
        mut_props = self.aa_properties[mutant_aa]

        multiplier = 1.0

        # Charge changes are highly impactful for GOF
        charge_change = abs(mut_props['charge'] - orig_props['charge'])
        if charge_change > 0:
            multiplier *= (1.0 + charge_change * 0.3)  # Up to 60% boost for double charge change

        # Hydrophobic to polar changes (or vice versa) are impactful
        if orig_props['hydrophobic'] != mut_props['hydrophobic']:
            multiplier *= 1.2  # 20% boost

        # Size changes matter
        size_change = abs(mut_props['size'] - orig_props['size'])
        if size_change > 2:
            multiplier *= 1.15  # 15% boost for large size changes

        return min(multiplier, 1.5)  # Cap at 50% total boost

    def _analyze_regulatory_context_disruption(self, original_aa: str, mutant_aa: str, position: int,
                                             sequence: str) -> Dict[str, float]:
        """
        REVOLUTIONARY REGULATORY CONTEXT ANALYSIS!

        Analyzes how this specific change disrupts regulatory mechanisms
        NO HARDCODING - uses universal regulatory patterns!
        """
        context_scores = {
            'phosphorylation_disruption': 0.0,
            'charge_regulatory_disruption': 0.0,
            'flexibility_regulatory_disruption': 0.0,
            'allosteric_disruption': 0.0,
            'binding_interface_disruption': 0.0
        }

        # 1. PHOSPHORYLATION SITE DISRUPTION - Universal brake pedal detection!
        context_scores['phosphorylation_disruption'] = self._analyze_phosphorylation_disruption(
            original_aa, mutant_aa, position, sequence
        )

        # 2. CHARGE DISRUPTION IN REGULATORY CONTEXTS
        context_scores['charge_regulatory_disruption'] = self._analyze_charge_regulatory_disruption(
            original_aa, mutant_aa, position, sequence
        )

        # 3. FLEXIBILITY DISRUPTION IN REGULATORY REGIONS
        context_scores['flexibility_regulatory_disruption'] = self._analyze_flexibility_regulatory_disruption(
            original_aa, mutant_aa, position, sequence
        )

        # 4. ALLOSTERIC NETWORK DISRUPTION
        context_scores['allosteric_disruption'] = self._analyze_allosteric_disruption(
            original_aa, mutant_aa, position, sequence
        )

        # 5. BINDING INTERFACE DISRUPTION
        context_scores['binding_interface_disruption'] = self._analyze_binding_interface_disruption(
            original_aa, mutant_aa, position, sequence
        )

        return context_scores

    def _analyze_phosphorylation_disruption(self, original_aa: str, mutant_aa: str,
                                          position: int, sequence: str) -> float:
        """
        UNIVERSAL PHOSPHORYLATION BRAKE PEDAL DETECTION!

        Detects disruption of phosphorylation sites - the ultimate GOF mechanism!
        """
        score = 0.0

        # Check if we're losing a phosphorylation site (BRAKE PEDAL REMOVAL!)
        if original_aa in ['S', 'T', 'Y'] and mutant_aa not in ['S', 'T', 'Y']:
            # Losing a potential phosphorylation site!

            # Check for kinase consensus sequences around this position
            window_start = max(0, position - 6)
            window_end = min(len(sequence), position + 5)
            local_sequence = sequence[window_start:window_end]

            # Universal kinase consensus patterns (NO HARDCODING!)
            kinase_score = self._detect_kinase_consensus(local_sequence, position - window_start)

            if kinase_score > 0.5:
                score = 0.9  # VERY HIGH - losing a real phosphorylation site!
                logger.info(f"🎯 PHOSPHO BRAKE PEDAL REMOVAL: {original_aa}{position}{mutant_aa} in kinase consensus!")
            elif kinase_score > 0.2:
                score = 0.6  # HIGH - losing a potential phosphorylation site
                logger.info(f"🎯 Potential phospho site loss: {original_aa}{position}{mutant_aa}")
            else:
                score = 0.3  # MODERATE - losing phospho potential

        # Check if we're gaining a phosphorylation site (NEW BRAKE PEDAL!)
        elif original_aa not in ['S', 'T', 'Y'] and mutant_aa in ['S', 'T', 'Y']:
            # This could create a new regulatory site - usually less GOF potential
            score = 0.1  # Low GOF potential (new brakes usually reduce GOF)

        return score

    def _detect_kinase_consensus(self, local_sequence: str, target_pos: int) -> float:
        """
        Universal kinase consensus detection - NO HARDCODING!

        Detects common kinase targeting patterns around phosphorylation sites
        """
        if target_pos >= len(local_sequence):
            return 0.0

        consensus_score = 0.0

        # PKA consensus: R/K-R/K-X-S/T
        if target_pos >= 3:
            upstream = local_sequence[target_pos-3:target_pos]
            if len(upstream) >= 2:
                basic_count = sum(1 for aa in upstream[-2:] if aa in 'RK')
                if basic_count >= 1:
                    consensus_score += 0.4
                if basic_count >= 2:
                    consensus_score += 0.3

        # CK2 consensus: S/T-X-X-E/D or E/D-X-X-S/T
        if target_pos >= 3 and target_pos < len(local_sequence) - 3:
            downstream = local_sequence[target_pos+1:target_pos+4]
            upstream = local_sequence[target_pos-3:target_pos]

            # Check downstream acidic
            if len(downstream) >= 2 and any(aa in 'ED' for aa in downstream):
                consensus_score += 0.3

            # Check upstream acidic
            if len(upstream) >= 2 and any(aa in 'ED' for aa in upstream):
                consensus_score += 0.3

        # Proline-directed kinases: S/T-P
        if target_pos < len(local_sequence) - 1:
            if local_sequence[target_pos + 1] == 'P':
                consensus_score += 0.4

        return min(consensus_score, 1.0)

    def _analyze_charge_regulatory_disruption(self, original_aa: str, mutant_aa: str,
                                            position: int, sequence: str) -> float:
        """
        Analyze charge changes in regulatory contexts
        """
        if original_aa not in self.aa_properties or mutant_aa not in self.aa_properties:
            return 0.0

        orig_charge = self.aa_properties[original_aa]['charge']
        mut_charge = self.aa_properties[mutant_aa]['charge']
        charge_change = abs(mut_charge - orig_charge)

        if charge_change == 0:
            return 0.0

        # Analyze local charge environment
        window_start = max(0, position - 10)
        window_end = min(len(sequence), position + 10)
        local_sequence = sequence[window_start:window_end]

        # Count charged residues nearby
        charged_nearby = sum(1 for aa in local_sequence if aa in 'DEKR')
        charge_density = charged_nearby / len(local_sequence)

        # High charge density = likely regulatory region
        if charge_density > 0.3:
            score = charge_change * 0.7  # High impact in charge-rich regions
            logger.info(f"🎯 Charge disruption in regulatory region: {original_aa}{position}{mutant_aa}")
        else:
            score = charge_change * 0.3  # Lower impact in charge-poor regions

        return min(score, 1.0)

    def _analyze_flexibility_regulatory_disruption(self, original_aa: str, mutant_aa: str,
                                                 position: int, sequence: str) -> float:
        """
        Analyze flexibility changes in regulatory contexts
        Gly->anything in hinge regions = conformational locking!
        """
        if original_aa not in self.aa_properties or mutant_aa not in self.aa_properties:
            return 0.0

        orig_flex = self.aa_properties[original_aa]['flexibility']
        mut_flex = self.aa_properties[mutant_aa]['flexibility']

        flexibility_map = {'low': 1, 'medium': 2, 'high': 3, 'rigid': 0}
        orig_flex_score = flexibility_map.get(orig_flex, 2)
        mut_flex_score = flexibility_map.get(mut_flex, 2)

        flexibility_change = orig_flex_score - mut_flex_score  # Positive = losing flexibility

        if flexibility_change <= 0:
            return 0.0  # Gaining flexibility rarely causes GOF

        # Special case: Glycine loss (ultimate flexibility loss!)
        if original_aa == 'G':
            # Glycine provides unique conformational freedom
            # Losing it can lock conformations - MAJOR GOF potential!

            # Check if we're in a potential hinge region
            hinge_score = self._detect_hinge_region(position, sequence)

            if hinge_score > 0.5:
                score = 0.8  # VERY HIGH - Gly loss in hinge region!
                logger.info(f"CONFORMATIONAL LOCK: Gly{position}{mutant_aa} in hinge region!")
            else:
                score = 0.5  # HIGH - Gly loss anywhere is significant
                logger.info(f"Flexibility loss: Gly{position}{mutant_aa}")

        # Proline introduction (rigidity introduction)
        elif mutant_aa == 'P':
            score = flexibility_change * 0.4  # Proline can lock conformations
            logger.info(f"Rigidity introduction: {original_aa}{position}Pro")

        else:
            score = flexibility_change * 0.2  # General flexibility loss

        return min(score, 1.0)

    def _detect_hinge_region(self, position: int, sequence: str) -> float:
        """
        Detect potential hinge regions using universal patterns
        """
        # Hinge regions often have:
        # 1. High glycine content
        # 2. Low hydrophobic content
        # 3. Mixed charge patterns

        window_start = max(0, position - 5)
        window_end = min(len(sequence), position + 6)
        local_sequence = sequence[window_start:window_end]

        if len(local_sequence) < 3:
            return 0.0

        # Count glycines (hinge indicators)
        gly_count = local_sequence.count('G')
        gly_density = gly_count / len(local_sequence)

        # Count hydrophobic residues (structural indicators)
        hydrophobic_count = sum(1 for aa in local_sequence if aa in 'AILMFWV')
        hydrophobic_density = hydrophobic_count / len(local_sequence)

        # Hinge score: high gly, low hydrophobic
        hinge_score = gly_density * 0.7 + (1 - hydrophobic_density) * 0.3

        return min(hinge_score, 1.0)

    def _analyze_allosteric_disruption(self, original_aa: str, mutant_aa: str,
                                     position: int, sequence: str) -> float:
        """
        Analyze potential allosteric network disruption
        """
        # Allosteric networks often involve:
        # 1. Conserved hydrophobic patches
        # 2. Salt bridge networks
        # 3. Aromatic stacking interactions

        if original_aa not in self.aa_properties or mutant_aa not in self.aa_properties:
            return 0.0

        score = 0.0

        # Aromatic residue changes (often allosteric)
        if original_aa in 'FWY' and mutant_aa not in 'FWY':
            score += 0.4  # Losing aromatic interactions
            logger.info(f"Aromatic loss: {original_aa}{position}{mutant_aa}")
        elif original_aa not in 'FWY' and mutant_aa in 'FWY':
            score += 0.2  # Gaining aromatic interactions (less common GOF)

        # Hydrophobic patch disruption
        orig_hydrophobic = self.aa_properties[original_aa]['hydrophobic']
        mut_hydrophobic = self.aa_properties[mutant_aa]['hydrophobic']

        if orig_hydrophobic and not mut_hydrophobic:
            # Disrupting hydrophobic patch
            hydrophobic_context = self._analyze_hydrophobic_context(position, sequence)
            score += hydrophobic_context * 0.3

        return min(score, 1.0)

    def _analyze_hydrophobic_context(self, position: int, sequence: str) -> float:
        """
        Analyze hydrophobic context around position
        """
        window_start = max(0, position - 3)
        window_end = min(len(sequence), position + 4)
        local_sequence = sequence[window_start:window_end]

        hydrophobic_count = sum(1 for aa in local_sequence if aa in 'AILMFWV')
        return hydrophobic_count / len(local_sequence)

    def _analyze_binding_interface_disruption(self, original_aa: str, mutant_aa: str,
                                            position: int, sequence: str) -> float:
        """
        Analyze potential protein-protein interface disruption
        """
        # Interface residues often:
        # 1. Are surface exposed
        # 2. Have specific charge/hydrophobic patterns
        # 3. Are in beta-sheets or loops

        if original_aa not in self.aa_properties or mutant_aa not in self.aa_properties:
            return 0.0

        # Simple interface prediction based on amino acid properties
        interface_score = 0.0

        # Charge changes at interfaces are highly disruptive
        orig_charge = self.aa_properties[original_aa]['charge']
        mut_charge = self.aa_properties[mutant_aa]['charge']
        charge_change = abs(mut_charge - orig_charge)

        if charge_change > 0:
            interface_score += charge_change * 0.4

        # Size changes at interfaces
        orig_size = self.aa_properties[original_aa]['size']
        mut_size = self.aa_properties[mutant_aa]['size']
        size_change = abs(mut_size - orig_size)

        if size_change > 1:
            interface_score += (size_change / 5.0) * 0.3

        return min(interface_score, 1.0)

    def _calculate_regulatory_disruption_potential(self, original_aa: str, mutant_aa: str,
                                                 position: int, sequence: str) -> float:
        """
        REVOLUTIONARY REGULATORY DISRUPTION POTENTIAL CALCULATOR!

        This is the ONLY gate that matters for GOF - forget Grantham!
        Calculates how likely this change is to disrupt regulatory mechanisms
        """
        disruption_score = 0.0

        # NOVA'S MOTIF DETECTION SYSTEM - Universal regulatory motif detection!
        # try:
        #     from AdaptiveInterpreter.utils.motif_detector import detect_regulatory_context
        #
        #     motif_context = detect_regulatory_context(sequence, original_aa, mutant_aa, position)
        #     motif_score = motif_context.get('score', 0.0)
        #
        #     if motif_score > disruption_score:
        #         disruption_score = motif_score
        #         logger.info(f"NOVA'S MOTIF DETECTION: {motif_context.get('explanation', 'Unknown')}")
        #
        #         # If Nova's motif detection gives a high score (0.8+), this is likely a canonical GOF
        #         # Skip the complex mechanism analysis and return the motif score directly
        #         if motif_score >= 0.8:
        #             logger.info(f"HIGH MOTIF SCORE ({motif_score:.3f}) - Canonical GOF detected!")
        #             return motif_score
        #
        # except Exception as e:
        #     logger.warning(f"Motif detection failed: {e}")
        #     # Continue with original logic as fallback

        # 1. PHOSPHORYLATION SITE DISRUPTION - The ultimate GOF trigger!
        phospho_disruption = self._analyze_phosphorylation_disruption(original_aa, mutant_aa, position, sequence)
        if phospho_disruption > 0.5:
            disruption_score = max(disruption_score, 0.9)  # MASSIVE disruption potential!
            logger.info(f"PHOSPHO BRAKE PEDAL DISRUPTION: {original_aa}{position}{mutant_aa}")
        elif phospho_disruption > 0.2:
            disruption_score = max(disruption_score, 0.6)  # High disruption potential

        # 2. GLYCINE LOSS - Conformational lock potential
        if original_aa == 'G':
            hinge_score = self._detect_hinge_region(position, sequence)
            if hinge_score > 0.5:
                disruption_score = max(disruption_score, 0.8)  # Very high - Gly loss in hinge!
                logger.info(f"CONFORMATIONAL LOCK: Gly{position}{mutant_aa} in hinge region!")
            else:
                disruption_score = max(disruption_score, 0.4)  # Moderate - Gly loss anywhere

        # 3. CHARGE CHANGES IN REGULATORY CONTEXTS
        charge_disruption = self._analyze_charge_regulatory_disruption(original_aa, mutant_aa, position, sequence)
        if charge_disruption > 0.5:
            disruption_score = max(disruption_score, 0.7)  # High regulatory disruption
            logger.info(f"CHARGE REGULATORY DISRUPTION: {original_aa}{position}{mutant_aa}")
        elif charge_disruption > 0.3:
            disruption_score = max(disruption_score, 0.4)  # Moderate disruption

        # 4. PROLINE INTRODUCTION - Rigidity introduction
        if mutant_aa == 'P':
            disruption_score = max(disruption_score, 0.5)  # Moderate - can lock conformations
            logger.info(f"RIGIDITY INTRODUCTION: {original_aa}{position}Pro")

        # 5. AROMATIC CHANGES - Allosteric disruption
        if (original_aa in 'FWY' and mutant_aa not in 'FWY') or (original_aa not in 'FWY' and mutant_aa in 'FWY'):
            disruption_score = max(disruption_score, 0.3)  # Moderate allosteric potential
            logger.info(f"AROMATIC CHANGE: {original_aa}{position}{mutant_aa}")

        # 6. CYSTEINE CHANGES - Disulfide bond disruption
        if (original_aa == 'C' and mutant_aa != 'C') or (original_aa != 'C' and mutant_aa == 'C'):
            disruption_score = max(disruption_score, 0.4)  # Moderate - structural/regulatory
            logger.info(f"CYSTEINE CHANGE: {original_aa}{position}{mutant_aa}")

        # 7. HYDROPHOBIC PATCH DISRUPTION
        if original_aa in 'AILMFWV' and mutant_aa not in 'AILMFWV':
            hydrophobic_context = self._analyze_hydrophobic_context(position, sequence)
            if hydrophobic_context > 0.6:
                disruption_score = max(disruption_score, 0.3)  # Moderate hydrophobic patch disruption

        logger.info(f"Final regulatory disruption score: {disruption_score:.3f}")
        return disruption_score

    def _run_regulatory_gof_analysis(self, original_aa: str, mutant_aa: str, position: int,
                                   sequence: str) -> Dict[str, float]:
        """
        REVOLUTIONARY REGULATORY-FOCUSED GOF ANALYSIS!

        Analyzes GOF mechanisms through the lens of regulatory disruption
        """
        # Get all regulatory context scores
        context_scores = self._analyze_regulatory_context_disruption(
            original_aa, mutant_aa, position, sequence
        )

        # Convert regulatory disruptions into GOF mechanism scores
        gof_scores = {}

        # 1. Phosphorylation disruption -> Constitutive Activation + Autoinhibition Loss
        phospho_score = context_scores['phosphorylation_disruption']
        gof_scores['constitutive_activation'] = phospho_score * 0.9  # Phospho loss = always on
        gof_scores['autoinhibition_loss'] = phospho_score * 0.95     # Phospho loss = brake removal

        # 2. Flexibility disruption -> Constitutive Activation
        flexibility_score = context_scores['flexibility_regulatory_disruption']
        gof_scores['constitutive_activation'] = max(
            gof_scores.get('constitutive_activation', 0),
            flexibility_score * 0.8
        )

        # 3. Charge disruption -> Binding Affinity + Autoinhibition Loss
        charge_score = context_scores['charge_regulatory_disruption']
        gof_scores['increased_binding_affinity'] = charge_score * 0.7
        gof_scores['autoinhibition_loss'] = max(
            gof_scores.get('autoinhibition_loss', 0),
            charge_score * 0.6
        )

        # 4. Allosteric disruption -> All mechanisms (moderate)
        allosteric_score = context_scores['allosteric_disruption']
        for mechanism in ['constitutive_activation', 'increased_binding_affinity', 'autoinhibition_loss']:
            gof_scores[mechanism] = max(
                gof_scores.get(mechanism, 0),
                allosteric_score * 0.4
            )

        # 5. Binding interface disruption -> Binding Affinity
        interface_score = context_scores['binding_interface_disruption']
        gof_scores['increased_binding_affinity'] = max(
            gof_scores.get('increased_binding_affinity', 0),
            interface_score * 0.6
        )

        # 6. Degradation resistance (less common for regulatory changes)
        gof_scores['degradation_resistance'] = 0.0  # Regulatory changes rarely affect stability

        # Add context scores for transparency
        gof_scores.update({
            f"context_{key}": value for key, value in context_scores.items()
        })

        return gof_scores

    def _calculate_regulatory_gof_score(self, gof_scores: Dict[str, float],
                                      regulatory_disruption_score: float) -> float:
        """
        Calculate overall GOF score based on regulatory mechanisms
        """
        # Focus on the core GOF mechanisms
        core_mechanisms = {
            'constitutive_activation': 0.35,
            'increased_binding_affinity': 0.25,
            'autoinhibition_loss': 0.35,
            'degradation_resistance': 0.05
        }

        weighted_score = 0.0
        for mechanism, weight in core_mechanisms.items():
            if mechanism in gof_scores:
                weighted_score += gof_scores[mechanism] * weight

        # Boost by overall regulatory disruption potential
        regulatory_boost = 1.0 + (regulatory_disruption_score * 0.5)  # Up to 50% boost

        final_score = weighted_score * regulatory_boost

        return min(final_score, 1.0)  # Cap at 1.0

    def _get_structural_confidence_with_fallbacks(self, uniprot_id: str, position: int) -> Dict:
        """
        🛡️ BULLETPROOF structural confidence analysis with smart fallbacks

        Hierarchy:
        1. AlphaFold pLDDT (best) - When position is covered
        2. Domain annotations (good) - When we know the domain type
        3. Conservation scores (okay) - When Ensembl works
        4. Protein type defaults (fallback) - Always available
        5. Emergency fallback (never fails) - 1.0x multiplier

        Returns:
            Dict with score, source, confidence, multiplier, rationale
        """

        # LEVEL 1: Try AlphaFold pLDDT (best quality)
        try:
            plddt_score = self._get_alphafold_plddt(uniprot_id, position)
            if plddt_score is not None:
                return {
                    'score': plddt_score,
                    'source': 'alphafold_plddt',
                    'confidence': 'high',
                    'multiplier': self._calculate_plddt_gof_multiplier(plddt_score),
                    'rationale': f'AlphaFold pLDDT confidence score for structured region'
                }
        except Exception as e:
            logger.debug(f"AlphaFold pLDDT failed: {e}")

        # LEVEL 2: Try domain-based confidence estimation
        try:
            domain_result = self._estimate_domain_confidence(uniprot_id, position)
            if domain_result:
                return domain_result
        except Exception as e:
            logger.debug(f"Domain analysis failed: {e}")

        # LEVEL 3: Try conservation if available and not offline
        if not self.offline_mode and self.conservation_db:
            try:
                conservation_data = self.conservation_db.get_variant_conservation(uniprot_id, position)
                if conservation_data and 'conservation_scores' in conservation_data and conservation_data['conservation_scores']:
                    phylop = conservation_data['conservation_scores']['phyloP']
                    return {
                        'score': phylop * 10,  # Convert to 0-100 scale
                        'source': 'phylop_conservation',
                        'confidence': 'medium',
                        'multiplier': self._calculate_conservation_gof_multiplier(phylop),
                        'rationale': f'PhyloP evolutionary conservation score'
                    }
            except Exception as e:
                logger.debug(f"Conservation analysis failed: {e}")

        # LEVEL 4: Protein type-based defaults
        try:
            protein_default = self._default_confidence_by_protein_type(uniprot_id)
            if protein_default:
                return protein_default
        except Exception as e:
            logger.debug(f"Protein type analysis failed: {e}")

        # LEVEL 5: Emergency fallback (never fails!)
        return {
            'score': 50.0,
            'source': 'emergency_fallback',
            'confidence': 'none',
            'multiplier': 1.0,
            'rationale': 'No structural or evolutionary information available - using neutral multiplier'
        }

    def _get_alphafold_plddt(self, uniprot_id: str, position: int) -> float:
        """
        Extract pLDDT confidence score from AlphaFold structure

        Args:
            uniprot_id: UniProt ID (e.g., 'P04637')
            position: Protein position (1-based)

        Returns:
            pLDDT score (0-100) or None if not available
        """
        try:
            import gzip
            from pathlib import Path

            # AlphaFold structure path
            alphafold_path = Path("/mnt/Arcana/alphafold_human/structures")
            pdb_file = alphafold_path / f"AF-{uniprot_id}-F1-model_v4.pdb.gz"

            if not pdb_file.exists():
                logger.debug(f"AlphaFold structure not found: {pdb_file}")
                return None

            # Parse PDB file to extract pLDDT for specific position
            with gzip.open(pdb_file, 'rt') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        # Parse PDB ATOM line
                        atom_name = line[12:16].strip()
                        residue_num = int(line[22:26].strip())
                        plddt_score = float(line[60:66].strip())  # B-factor column = pLDDT

                        # Get pLDDT for CA atom of target residue
                        if atom_name == 'CA' and residue_num == position:
                            logger.debug(f"Found pLDDT for {uniprot_id}:{position} = {plddt_score:.1f}")
                            return plddt_score

            logger.debug(f"Position {position} not found in {uniprot_id} structure")
            return None

        except Exception as e:
            logger.warning(f"Failed to extract pLDDT for {uniprot_id}:{position}: {e}")
            return None

    def _calculate_plddt_gof_multiplier(self, plddt: float) -> float:
        """
        Calculate GOF multiplier based on AlphaFold pLDDT confidence

        Logic: High confidence regions (blue) = mutations more disruptive
               Low confidence regions (orange) = mutations less impactful

        Args:
            plddt: pLDDT score (0-100)

        Returns:
            Multiplier for GOF score (0.5-2.0)
        """
        if plddt >= 90:
            # Very high confidence (blue) - mutations very disruptive
            multiplier = 2.0
            logger.debug(f"VERY HIGH CONFIDENCE: pLDDT {plddt:.1f} → {multiplier:.1f}x GOF boost!")
        elif plddt >= 70:
            # Confident (green) - mutations moderately disruptive
            multiplier = 1.5
            logger.debug(f"HIGH CONFIDENCE: pLDDT {plddt:.1f} → {multiplier:.1f}x GOF boost!")
        elif plddt >= 50:
            # Low confidence (yellow) - mutations less impactful
            multiplier = 1.0
            logger.debug(f"LOW CONFIDENCE: pLDDT {plddt:.1f} → {multiplier:.1f}x (no boost)")
        else:
            # Very low confidence (orange) - mutations minimal impact
            multiplier = 0.8
            logger.debug(f"VERY LOW CONFIDENCE: pLDDT {plddt:.1f} → {multiplier:.1f}x (reduced)")

        return multiplier

    def _estimate_domain_confidence(self, uniprot_id: str, position: int) -> Dict:
        """
        Estimate structural confidence based on known domain types and annotations
        """
        try:
            # Try to get UniProt annotations (we have this in our system!)
            from nova_dn.universal_context import UniversalContext
            context = UniversalContext()

            # Get protein annotations
            annotations = context.get_context_for_protein("", uniprot_id)  # Gene name not needed

            if "error" in annotations:
                return None

            function_text = annotations.get("function", "").lower()
            go_terms = [term.lower() for term in annotations.get("go_terms", [])]

            # Check for low-confidence regions (disordered/flexible)
            low_confidence_keywords = [
                'transmembrane', 'signal peptide', 'coiled coil', 'linker',
                'disordered', 'flexible', 'loop', 'unstructured'
            ]

            for keyword in low_confidence_keywords:
                if keyword in function_text or any(keyword in term for term in go_terms):
                    return {
                        'score': 35.0,  # Orange zone
                        'source': 'domain_annotation_low',
                        'confidence': 'medium',
                        'multiplier': 0.8,
                        'rationale': f'Position likely in disordered region ({keyword} detected)'
                    }

            # Check for high-confidence regions (structured/functional)
            high_confidence_keywords = [
                'catalytic', 'binding', 'active site', 'zinc finger', 'domain',
                'kinase', 'phosphatase', 'enzyme', 'structural'
            ]

            for keyword in high_confidence_keywords:
                if keyword in function_text or any(keyword in term for term in go_terms):
                    return {
                        'score': 80.0,  # Green zone
                        'source': 'domain_annotation_high',
                        'confidence': 'medium',
                        'multiplier': 1.4,
                        'rationale': f'Position likely in structured region ({keyword} detected)'
                    }

            # Medium confidence default for annotated proteins
            return {
                'score': 60.0,  # Yellow zone
                'source': 'domain_annotation_medium',
                'confidence': 'low',
                'multiplier': 1.1,
                'rationale': 'Protein has functional annotations but no specific domain info'
            }

        except Exception as e:
            logger.debug(f"🔍 Domain confidence estimation failed: {e}")
            return None

    def _default_confidence_by_protein_type(self, uniprot_id: str) -> Dict:
        """
        Default confidence based on protein family/type - ALWAYS returns something!
        """
        try:
            from nova_dn.universal_context import UniversalContext
            context = UniversalContext()

            # Get protein annotations
            annotations = context.get_context_for_protein("", uniprot_id)

            if "error" not in annotations:
                function = annotations.get('function', '').lower()

                # Membrane proteins - mixed confidence
                if any(keyword in function for keyword in ['channel', 'receptor', 'transporter', 'membrane']):
                    return {
                        'score': 55.0,  # Yellow zone
                        'source': 'protein_type_membrane',
                        'confidence': 'low',
                        'multiplier': 1.0,
                        'rationale': 'Membrane protein (mixed structured/disordered regions)'
                    }

                # Enzymes - generally well-structured
                elif any(keyword in function for keyword in ['enzyme', 'kinase', 'phosphatase', 'catalytic']):
                    return {
                        'score': 70.0,  # Green zone
                        'source': 'protein_type_enzyme',
                        'confidence': 'low',
                        'multiplier': 1.2,
                        'rationale': 'Enzyme (typically well-structured active sites)'
                    }

                # Structural proteins - high confidence
                elif any(keyword in function for keyword in ['collagen', 'structural', 'fibrillar', 'matrix']):
                    return {
                        'score': 75.0,  # Green zone
                        'source': 'protein_type_structural',
                        'confidence': 'low',
                        'multiplier': 1.3,
                        'rationale': 'Structural protein (typically well-organized)'
                    }

                # Transcription factors - mixed (structured domains + flexible linkers)
                elif any(keyword in function for keyword in ['transcription', 'dna binding', 'regulatory']):
                    return {
                        'score': 65.0,  # Green-yellow zone
                        'source': 'protein_type_transcription',
                        'confidence': 'low',
                        'multiplier': 1.1,
                        'rationale': 'Transcription factor (structured domains + flexible regions)'
                    }

        except Exception as e:
            logger.debug(f"🔍 Protein type analysis failed: {e}")

        # Ultimate fallback - NEVER fails!
        return {
            'score': 50.0,  # Neutral
            'source': 'global_default',
            'confidence': 'none',
            'multiplier': 1.0,
            'rationale': 'No protein information available - using neutral confidence'
        }

    def _calculate_conservation_gof_multiplier(self, phylop: float) -> float:
        """
        REVOLUTIONARY CONSERVATION-BASED GOF MULTIPLIER!

        Ultra-conserved positions get MASSIVE GOF boosts because ANY change is devastating!
        This catches variants like V1363L that are pathogenic despite low Grantham distance.
        """
        if phylop > 5.0:
            # EXTREMELY conserved - ANY change is devastating!
            multiplier = 3.0  # MASSIVE boost for ultra-conserved positions!
            logger.info(f"ULTRA-CONSERVED POSITION: PhyloP {phylop:.3f} → {multiplier:.1f}x GOF boost!")
        elif phylop > 2.0:
            # Highly conserved - significant GOF potential
            multiplier = 2.0  # Major boost
            logger.info(f"HIGHLY CONSERVED: PhyloP {phylop:.3f} → {multiplier:.1f}x GOF boost!")
        elif phylop > 1.0:
            # Moderately conserved - moderate GOF boost
            multiplier = 1.5  # Moderate boost
            logger.info(f"MODERATELY CONSERVED: PhyloP {phylop:.3f} → {multiplier:.1f}x GOF boost!")
        elif phylop > 0.5:
            # Somewhat conserved - minor boost
            multiplier = 1.2  # Minor boost
        else:
            # Not conserved - no boost
            multiplier = 1.0  # No boost

        return multiplier

    def _run_enhanced_regulatory_analysis(self, original_aa: str, mutant_aa: str, position: int,
                                        conservation_enhanced_disruption: float,
                                        regulatory_gof_scores: Dict[str, float],
                                        conservation_multiplier: float,
                                        conservation_data: Dict) -> Dict[str, Any]:
        """
        Enhanced regulatory analysis for sequence mismatch cases
        """
        # Apply additional regulatory enhancement
        enhanced_scores = regulatory_gof_scores.copy()

        # Boost scores based on conservation-enhanced regulatory disruption potential
        regulatory_multiplier = 1.0 + conservation_enhanced_disruption  # Up to 2x boost

        core_mechanisms = ['constitutive_activation', 'increased_binding_affinity',
                          'autoinhibition_loss', 'degradation_resistance']

        for mechanism in core_mechanisms:
            if mechanism in enhanced_scores:
                enhanced_scores[mechanism] *= regulatory_multiplier
                enhanced_scores[mechanism] = min(enhanced_scores[mechanism], 1.0)

        # Calculate enhanced overall score
        overall_score = self._calculate_regulatory_gof_score(enhanced_scores, conservation_enhanced_disruption)

        return {
            'mutation': f"{original_aa}{position}{mutant_aa}",
            'grantham_distance': self.get_grantham_distance(original_aa, mutant_aa),
            'conservation_enhanced_disruption': conservation_enhanced_disruption,
            'conservation_multiplier': conservation_multiplier,
            'conservation_data': conservation_data,
            'gof_mechanisms': enhanced_scores,
            'gof_score': overall_score,
            'prediction': 'GOF_LIKELY' if overall_score > 0.6 else 'GOF_POSSIBLE' if overall_score > 0.3 else 'GOF_UNLIKELY',
            'confidence': 0.95,  # Highest confidence with conservation analysis
            'analysis_level': 'GATE_3_CONSERVATION_ENHANCED_REGULATORY_SEQUENCE_MISMATCH',
            'sequence_mismatch': True,
            'regulatory_multiplier': regulatory_multiplier
        }

    def _run_structural_regulatory_analysis(self, original_aa: str, mutant_aa: str, position: int,
                                          sequence: str, uniprot_id: str,
                                          regulatory_disruption_score: float,
                                          regulatory_gof_scores: Dict[str, float],
                                          conservation_multiplier: float = 1.0,
                                          conservation_data: Dict = None) -> Dict[str, Any]:
        """
        Combined structural and regulatory analysis for perfect sequence matches
        """
        try:
            # Start with regulatory scores as baseline
            final_scores = regulatory_gof_scores.copy()

            # TODO: Add structural analysis enhancement here
            # This would analyze AlphaFold structures to refine regulatory predictions

            # For now, apply modest structural context enhancement
            structural_enhancement = self._estimate_structural_context(position, sequence, uniprot_id)

            core_mechanisms = ['constitutive_activation', 'increased_binding_affinity',
                              'autoinhibition_loss', 'degradation_resistance']

            for mechanism in core_mechanisms:
                if mechanism in final_scores:
                    final_scores[mechanism] *= structural_enhancement
                    final_scores[mechanism] = min(final_scores[mechanism], 1.0)

            # Calculate final score
            overall_score = self._calculate_regulatory_gof_score(final_scores, regulatory_disruption_score)

            return {
                'mutation': f"{original_aa}{position}{mutant_aa}",
                'grantham_distance': self.get_grantham_distance(original_aa, mutant_aa),
                'regulatory_disruption_score': regulatory_disruption_score,
                'gof_mechanisms': final_scores,
                'gof_score': overall_score,
                'prediction': 'GOF_LIKELY' if overall_score > 0.6 else 'GOF_POSSIBLE' if overall_score > 0.3 else 'GOF_UNLIKELY',
                'confidence': 0.95,  # Highest confidence with both regulatory and structural
                'analysis_level': 'GATE_3_STRUCTURAL_REGULATORY_ANALYSIS',
                'sequence_mismatch': False,
                'structural_enhancement': structural_enhancement,
                'regulatory_scores': regulatory_gof_scores  # Include original regulatory scores
            }

        except Exception as e:
            # Fallback to regulatory analysis
            overall_score = self._calculate_regulatory_gof_score(regulatory_gof_scores, regulatory_disruption_score)
            return {
                'mutation': f"{original_aa}{position}{mutant_aa}",
                'grantham_distance': self.get_grantham_distance(original_aa, mutant_aa),
                'regulatory_disruption_score': regulatory_disruption_score,
                'gof_mechanisms': regulatory_gof_scores,
                'gof_score': overall_score,
                'prediction': 'GOF_LIKELY' if overall_score > 0.6 else 'GOF_POSSIBLE' if overall_score > 0.3 else 'GOF_UNLIKELY',
                'confidence': 0.9,
                'analysis_level': 'GATE_3_REGULATORY_FALLBACK',
                'sequence_mismatch': False,
                'error': f'Structural analysis failed: {str(e)}'
            }

    def _get_domain_context(self, uniprot_id: str, gene_symbol: str) -> Dict[str, Any]:
        """Get domain context from UniProt for GOF analysis"""
        if not self.protein_annotator:
            return {"domains": [], "signal_peptide": [], "active_sites": [], "binding_sites": []}

        try:
            domain_data = self.protein_annotator.annotate_protein(gene_symbol, uniprot_id)
            return domain_data
        except Exception as e:
            logger.warning(f"Domain annotation error: {e}")
            return {"domains": [], "signal_peptide": [], "active_sites": [], "binding_sites": []}

    def _get_gof_domain_multiplier(self, position: int, domain_context: Dict[str, Any]) -> float:
        """Calculate domain-aware multiplier for GOF scoring

        GOF mechanisms are most dangerous in:
        - Regulatory regions (autoinhibition domains)
        - Active sites (constitutive activation)
        - Binding sites (increased affinity)
        - Less dangerous in cleavable regions (propeptides)
        """
        multiplier = 1.0

        # UNIVERSAL PROPEPTIDE LOGIC - Downweight cleavable regions
        for propeptide in domain_context.get("propeptides", []):
            if position in range(propeptide["start"], propeptide["end"]+1):
                if "n-terminal" in propeptide["description"].lower():
                    multiplier *= 0.5  # N-terminal propeptide - gets cleaved
                elif "c-terminal" in propeptide["description"].lower():
                    multiplier *= 0.3  # C-terminal propeptide - less critical for GOF

        # Signal peptides get cleaved off - minimal GOF impact
        for signal in domain_context.get("signal_peptide", []):
            if position in range(signal["start"], signal["end"]+1):
                multiplier *= 0.3

        # UPWEIGHT CRITICAL REGIONS FOR GOF MECHANISMS

        # Active sites - critical for constitutive activation
        for site in domain_context.get("active_sites", []):
            if position == site:  # Active sites are single positions
                multiplier *= 2.0  # Strong upweight for GOF activation

        # Binding sites - important for increased affinity
        for site in domain_context.get("binding_sites", []):
            if position == site:  # Binding sites are single positions
                multiplier *= 1.7  # Strong upweight for binding disruption

        # Functional regions - context-dependent
        for region in domain_context.get("regions", []):
            if position in range(region["start"], region["end"]+1):
                desc = region["description"].lower()
                if any(term in desc for term in ["regulatory", "autoinhibition", "inhibitory"]):
                    multiplier *= 1.8  # Regulatory regions - prime GOF targets
                elif any(term in desc for term in ["catalytic", "enzyme", "kinase"]):
                    multiplier *= 1.5  # Catalytic regions - good GOF targets
                elif "disordered" in desc:
                    multiplier *= 0.7  # Disordered regions - less GOF impact

        # Domains - generally important for GOF mechanisms
        for domain in domain_context.get("domains", []):
            if position in range(domain["start"], domain["end"]+1):
                desc = domain["description"].lower()
                if any(term in desc for term in ["kinase", "catalytic", "enzyme"]):
                    multiplier *= 1.4  # Catalytic domains - prime GOF targets
                elif any(term in desc for term in ["regulatory", "inhibitor"]):
                    multiplier *= 1.6  # Regulatory domains - excellent GOF targets
                elif any(term in desc for term in ["immunoglobulin", "fibronectin", "repeat"]):
                    multiplier *= 1.1  # Structural domains - mild GOF risk

        return max(multiplier, 0.1)  # Don't go below 0.1

    def _apply_domain_awareness_to_gof_score(self, gof_score: float, position: int,
                                           uniprot_id: str, gene_symbol: str) -> Tuple[float, float]:
        """Apply domain awareness to final GOF score"""
        if not uniprot_id or not self.protein_annotator:
            return gof_score, 1.0

        try:
            domain_context = self._get_domain_context(uniprot_id, gene_symbol or "")
            domain_multiplier = self._get_gof_domain_multiplier(position, domain_context)

            # Apply domain multiplier
            domain_aware_score = min(gof_score * domain_multiplier, 1.0)

            logger.info(f"GOF Domain multiplier for {gene_symbol or 'unknown'} position {position}: {domain_multiplier:.3f}")
            logger.info(f"   GOF score: {gof_score:.3f} → {domain_aware_score:.3f}")

            return domain_aware_score, domain_multiplier

        except Exception as e:
            logger.warning(f"Domain awareness failed: {e}")
            return gof_score, 1.0
