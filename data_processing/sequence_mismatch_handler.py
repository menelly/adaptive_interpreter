#!/usr/bin/env python3
"""
Universal Sequence Mismatch Handler
Handles transcript/isoform differences gracefully across all analyzers
Following Ren's Fundamental Laziness Principle: Do it right ONCE!
"""

import logging
from typing import Dict, Any, Optional, Tuple
import requests
import time

logger = logging.getLogger(__name__)

class SequenceMismatchHandler:
    """
    Universal handler for sequence mismatches across all analyzers
    Provides fallback strategies when sequence doesn't match expected variant
    """
    
    def __init__(self):
        self.ensembl_cache = {}
        self.uniprot_cache = {}
    
    def check_sequence_match(self, sequence: str, position: int, expected_aa: str, 
                           gene_symbol: str, uniprot_id: str) -> Dict[str, Any]:
        """
        Check if sequence matches expected amino acid at position
        Returns comprehensive mismatch information and fallback strategies
        """
        try:
            # Basic bounds check
            if position < 1 or position > len(sequence):
                return {
                    'match': False,
                    'error': f'Position {position} out of range for sequence length {len(sequence)}',
                    'fallback_strategy': 'mathematical_only',
                    'confidence_penalty': 0.3
                }
            
            # Check actual sequence
            actual_aa = sequence[position - 1]
            
            if actual_aa == expected_aa:
                return {
                    'match': True,
                    'actual_aa': actual_aa,
                    'expected_aa': expected_aa,
                    'fallback_strategy': 'none',
                    'confidence_penalty': 0.0
                }
            
            # Sequence mismatch detected
            import os
            if os.getenv('SEQUENCE_MISMATCH_WARN', '0') == '1':
                logger.warning(f"⚠️ Sequence mismatch: Expected {expected_aa} at position {position}, found {actual_aa}")
                logger.warning(f"⚠️ Gene: {gene_symbol}, UniProt: {uniprot_id}")

            # Try to determine cause and best fallback strategy
            mismatch_info = self._analyze_mismatch(
                gene_symbol, uniprot_id, position, expected_aa, actual_aa
            )
            
            return {
                'match': False,
                'actual_aa': actual_aa,
                'expected_aa': expected_aa,
                'fallback_strategy': mismatch_info['strategy'],
                'confidence_penalty': mismatch_info['penalty'],
                'likely_cause': mismatch_info['cause'],
                'recommendation': mismatch_info['recommendation']
            }
            
        except Exception as e:
            logger.error(f"Error in sequence mismatch check: {e}")
            return {
                'match': False,
                'error': str(e),
                'fallback_strategy': 'mathematical_only',
                'confidence_penalty': 0.4
            }
    
    def _analyze_mismatch(self, gene_symbol: str, uniprot_id: str, position: int,
                         expected_aa: str, actual_aa: str) -> Dict[str, Any]:
        """
        Analyze the type of mismatch and recommend best fallback strategy
        """
        # Conservative amino acid changes (likely benign mismatches)
        conservative_changes = {
            ('I', 'L'), ('L', 'I'),  # Isoleucine <-> Leucine
            ('F', 'Y'), ('Y', 'F'),  # Phenylalanine <-> Tyrosine  
            ('K', 'R'), ('R', 'K'),  # Lysine <-> Arginine
            ('D', 'E'), ('E', 'D'),  # Aspartic acid <-> Glutamic acid
            ('N', 'Q'), ('Q', 'N'),  # Asparagine <-> Glutamine
            ('S', 'T'), ('T', 'S'),  # Serine <-> Threonine
        }
        
        change_pair = (expected_aa, actual_aa)
        
        if change_pair in conservative_changes:
            return {
                'strategy': 'mathematical_with_boost',
                'penalty': 0.1,  # Minimal penalty
                'cause': 'conservative_isoform_difference',
                'recommendation': 'Use mathematical analysis with confidence boost'
            }
        
        # Check if it's a common polymorphism position
        if self._is_common_polymorphism_site(gene_symbol, position):
            return {
                'strategy': 'mathematical_with_population_data',
                'penalty': 0.15,
                'cause': 'common_polymorphism_site',
                'recommendation': 'Use population frequency data if available'
            }
        
        # Radical amino acid changes (more concerning)
        if self._is_radical_change(expected_aa, actual_aa):
            return {
                'strategy': 'mathematical_conservative',
                'penalty': 0.25,
                'cause': 'radical_isoform_difference',
                'recommendation': 'Use conservative mathematical analysis'
            }
        
        # Default case
        return {
            'strategy': 'mathematical_only',
            'penalty': 0.2,
            'cause': 'transcript_isoform_difference',
            'recommendation': 'Use mathematical analysis only'
        }
    
    def _is_radical_change(self, aa1: str, aa2: str) -> bool:
        """Check if amino acid change is radical (different properties)"""
        # Amino acid property groups
        hydrophobic = set('AILMFWYV')
        polar = set('NQST')
        charged_pos = set('KRH')
        charged_neg = set('DE')
        special = set('CGP')
        
        groups = [hydrophobic, polar, charged_pos, charged_neg, special]
        
        # Find which groups each amino acid belongs to
        aa1_groups = [i for i, group in enumerate(groups) if aa1 in group]
        aa2_groups = [i for i, group in enumerate(groups) if aa2 in group]
        
        # Radical if they don't share any property groups
        return len(set(aa1_groups) & set(aa2_groups)) == 0
    
    def _is_common_polymorphism_site(self, gene_symbol: str, position: int) -> bool:
        """Check if position is known to have common polymorphisms"""
        # Known common polymorphism sites
        common_sites = {
            'TP53': [72, 47, 213],  # P72R, S47, etc.
            'APOE': [112, 158],     # E2/E3/E4 variants
            'MTHFR': [222, 429],    # C677T, A1298C
        }
        
        return gene_symbol in common_sites and position in common_sites[gene_symbol]
    
    def apply_fallback_strategy(self, strategy: str, base_score: float, 
                              confidence_penalty: float, analysis_context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Apply the recommended fallback strategy to adjust scores and confidence
        """
        adjusted_score = base_score
        adjusted_confidence = max(0.1, analysis_context.get('confidence', 0.8) - confidence_penalty)
        
        if strategy == 'mathematical_with_boost':
            # Conservative mismatch - slight boost to mathematical analysis
            adjusted_score = min(1.0, base_score * 1.1)
            adjusted_confidence = max(adjusted_confidence, 0.7)
            
        elif strategy == 'mathematical_conservative':
            # Radical mismatch - be more conservative
            adjusted_score = base_score * 0.9
            adjusted_confidence = min(adjusted_confidence, 0.6)
            
        elif strategy == 'mathematical_only':
            # Standard mathematical analysis
            adjusted_score = base_score
            
        elif strategy == 'mathematical_with_population_data':
            # Try to incorporate population data if available
            # For now, just use standard approach
            adjusted_score = base_score
        
        return {
            'adjusted_score': adjusted_score,
            'adjusted_confidence': adjusted_confidence,
            'strategy_applied': strategy,
            'sequence_mismatch_handled': True
        }
    
    def get_alternative_sequences(self, gene_symbol: str, uniprot_id: str) -> Dict[str, str]:
        """
        Try to fetch alternative transcript sequences
        This is a placeholder for future enhancement
        """
        # Future: Could query Ensembl for alternative transcripts
        # For now, return empty dict
        return {}

def create_mismatch_handler() -> SequenceMismatchHandler:
    """Factory function to create sequence mismatch handler"""
    return SequenceMismatchHandler()
