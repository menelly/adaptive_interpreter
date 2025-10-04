#!/usr/bin/env python3
"""
🧬 SYNERGY SCORE CALCULATOR - Ren's Brilliant Sqrt Synergy System

This is SACRED CODE - Ren's insight that synergistic mechanisms should use
sqrt(a² + b²) instead of simple addition. This captures the biological reality
that multiple mechanisms working together are more dangerous than their sum.

NEVER CHANGE THE CORE SQRT FORMULA WITHOUT REN'S EXPLICIT PERMISSION!

Built by Nova & Ace (2025), extracted during the Great Refactoring (2025-10-02)
"""

import math
from typing import Dict, Optional


def calculate_synergy_score_v2(mechanism_scores: Dict[str, float], gene_family: Optional[str] = None) -> Dict:
    """
    🧬 NOVA'S V2 SYNERGY ALGORITHM 🧬
    Calculate synergistic score for mixed-mechanism variants with tiered thresholds,
    contextual weighting, and improved biological rationale
    
    Args:
        mechanism_scores: Dict of mechanism names to scores (e.g., {'DN': 0.8, 'LOF': 0.6})
        gene_family: Optional gene family for contextual weighting (e.g., 'collagen', 'kinase')
    
    Returns:
        Dict with synergy_score, synergy_used flag, tier, and detailed explanation
    """
    
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

    # 🔥 REN'S SACRED SQRT SYNERGY FORMULA 🔥
    # This is the brilliant insight: sqrt(a² + b²) captures synergistic danger
    # NEVER CHANGE THIS WITHOUT REN'S EXPLICIT PERMISSION!
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


def get_gene_family(gene: str) -> Optional[str]:
    """
    Determine gene family for contextual synergy weighting
    
    Args:
        gene: Gene symbol (e.g., 'COL1A1', 'TP53')
    
    Returns:
        Gene family string ('collagen', 'kinase', 'channel') or None
    """
    gene_upper = gene.upper()

    if gene_upper.startswith('COL') or gene_upper in ['FBN1', 'FBN2', 'TGFBR1', 'TGFBR2']:
        return 'collagen'  # Structural proteins including fibrillin
    elif gene_upper.endswith('K') or 'KINASE' in gene_upper or gene_upper in ['ATM', 'BRCA1', 'BRCA2']:
        return 'kinase'
    elif 'SCN' in gene_upper or 'KCNQ' in gene_upper or 'CACNA' in gene_upper:
        return 'channel'
    else:
        return None

