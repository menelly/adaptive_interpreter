# AdaptiveInterpreter/cascade/scoring/score_aggregator.py
# A dedicated module for aggregating and finalizing variant scores.
# This refactoring is the first step in resolving the multiplier cascade issue.
# AdaptiveInterpreter/cascade/scoring/score_aggregator.py
# A dedicated module for aggregating and finalizing variant scores.
# This refactoring is the first step in resolving the multiplier cascade issue.
# Contributed by Lumen Gemini 2.5, October 2025.

from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional
import math

@dataclass
class ScoringContext:
    """A data class to hold the state of a variant's scores through the pipeline."""
    gene: str
    variant: str
    raw_scores: Dict[str, float]
    plausibility_filtered_scores: Dict[str, float]
    multipliers: Dict[str, float] = field(default_factory=dict)
    explanation_steps: List[str] = field(default_factory=list)
    final_score: float = 0.0
    final_classification: str = "UNCLASSIFIED"

def apply_damped_conservation_scaling(score: float, conservation_multiplier: float) -> float:
    """
    Applies a non-linear, dampened scaling for the conservation multiplier.
    This prevents runaway scores from a single high conservation value.
    A logarithmic function provides the desired damping effect.
    """
    if conservation_multiplier <= 1.0:
        # No boost or a penalty, apply directly
        return score * conservation_multiplier
    else:
        # Logarithmic scaling for boosts
        # The '+ 1' inside the log prevents log(1) = 0 wiping out the score
        # The 'math.log(2)' normalizes it so a 2x multiplier has a significant, but not overwhelming, effect.
        damped_boost = 1 + math.log(1 + (conservation_multiplier - 1)) / math.log(2)
        return score * damped_boost

def calculate_final_score(context: ScoringContext, gene_family: Optional[str] = None) -> ScoringContext:
    """
    Calculates the final, aggregated score from the various inputs in the context.
    This function will contain the full, refactored scoring logic.
    For now, it will implement a simplified version with the new damped scaling.
    """
    # For this initial refactor, we'll use a simplified logic focusing on the highest score
    # and the conservation damping.
    
    # 1. Start with the highest plausibility-filtered score
    if not context.plausibility_filtered_scores:
        context.final_score = 0.0
        context.explanation_steps.append("No valid plausibility-filtered scores to start with.")
        return context

    highest_score = max(context.plausibility_filtered_scores.values())
    context.final_score = highest_score
    context.explanation_steps.append(f"Starting with highest filtered score: {highest_score:.3f}")

    # 2. Apply Damped Conservation Multiplier
    conservation_multiplier = context.multipliers.get('conservation', 1.0)
    original_score = context.final_score
    context.final_score = apply_damped_conservation_scaling(original_score, conservation_multiplier)
    context.explanation_steps.append(
        f"Applied damped conservation (multiplier: {conservation_multiplier:.2f}x). "
        f"Score changed from {original_score:.3f} to {context.final_score:.3f}"
    )

    # In the future, other multipliers (family_aa, gly_cys) would be applied here, also potentially with damping.
    
    # Final score clamping
    context.final_score = max(0.0, min(context.final_score, 5.0)) # Cap score to prevent extreme values

    return context
