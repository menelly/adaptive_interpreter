# utils.py
"""
🔧 NOVA'S LATTICE DISRUPTION UTILITIES
Shared functions for Nova's Universal Lattice Disruption Framework.

Prevents circular imports between analyzers and router.
"""

import math
from typing import List, Dict


def aggregate_lattice_features(features: List[Dict]) -> float:
    """
    🔥 NOVA'S ROOT-SUM-SQUARE AGGREGATION
    
    Combines multiple lattice disruption features using root-sum-square
    to prevent score inflation while allowing synergy.
    
    Args:
        features: List of feature dictionaries with 'weight' keys
        
    Returns:
        Aggregated score (0.0 to 1.0)
    """
    if not features:
        return 0.0
    
    # Root-sum-square aggregation
    sum_of_squares = sum(f["weight"]**2 for f in features)
    score = math.sqrt(sum_of_squares)
    
    # Cap at 1.0
    return min(score, 1.0)


def format_lattice_explanation(mechanism_name: str, features: List[Dict]) -> str:
    """
    📝 FORMAT LATTICE DISRUPTION EXPLANATION
    
    Create human-readable explanation of lattice disruption features.
    
    Args:
        mechanism_name: Name of the specific mechanism
        features: List of feature dictionaries
        
    Returns:
        Formatted explanation string
    """
    if not features:
        return f"{mechanism_name}: no disruption detected"
    
    # Sort features by weight (highest first)
    sorted_features = sorted(features, key=lambda f: f["weight"], reverse=True)
    
    # Format top features
    explanations = []
    for feature in sorted_features[:3]:  # Top 3 features
        desc = feature.get("description", feature["feature"])
        weight = feature["weight"]
        explanations.append(f"{desc} ({weight:.2f})")
    
    return f"{mechanism_name}: " + "; ".join(explanations)
