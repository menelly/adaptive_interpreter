#!/usr/bin/env python3
"""
🎯 VARIANT CLASSIFIER - Score to ACMG Classification Converter

Converts numeric pathogenicity scores to clinical ACMG classifications
with support for family-specific thresholds.

Built by Ace & Nova (2025), extracted during the Great Refactoring (2025-10-02)
"""

import json
from pathlib import Path
from typing import Dict, Optional


class VariantClassifier:
    """Converts numeric scores to ACMG classifications with family-specific thresholds"""
    
    def __init__(self, config_path: Optional[Path] = None):
        """
        Initialize classifier with optional custom threshold configuration
        
        Args:
            config_path: Optional path to classification_thresholds.json
        """
        self._family_thresholds = {}
        self._config_path = config_path
        self._load_classification_thresholds()
    
    def _load_classification_thresholds(self):
        """Load family-specific classification thresholds from config file"""
        if self._family_thresholds:  # Already loaded
            return
            
        self._family_thresholds = {}
        
        try:
            # Use provided path or default to cascade/resources/
            if self._config_path and self._config_path.exists():
                cfg_path = self._config_path
            else:
                cfg_path = Path(__file__).parent.parent / "resources" / "classification_thresholds.json"
            
            if cfg_path.exists():
                with cfg_path.open("r", encoding="utf-8") as fh:
                    data = json.load(fh)
                    if isinstance(data, dict):
                        self._family_thresholds = data
        except Exception:
            # Stay silent; fall back to defaults
            self._family_thresholds = {}
    
    def _get_family_thresholds(self, family: Optional[str]) -> Dict[str, float]:
        """
        Get classification thresholds for a specific gene family
        
        Args:
            family: Gene family name (e.g., 'collagen', 'tumor_suppressor')
        
        Returns:
            Dict of classification thresholds (P, LP, VUS, LB)
        """
        # Global defaults calibrated empirically (2025-12-11)
        # VUS-P collapsed into VUS for simpler, more honest classification
        # LP at 0.78 for high specificity (~75%), LB at 0.25 for conservative benign calls
        base = {"P": 1.2, "LP": 0.78, "VUS": 0.25, "LB": 0.2}

        if not family:
            return base
        
        self._load_classification_thresholds()
        fam = self._family_thresholds.get(family, {})
        
        # Merge partial overrides
        merged = dict(base)
        merged.update({k: v for k, v in fam.items() if isinstance(v, (int, float))})
        return merged
    
    def interpret_score(self, score: float, family: Optional[str] = None) -> str:
        """
        Convert numeric score to clinical classification with optional per-family thresholds
        
        Args:
            score: Numeric pathogenicity score (0.0-1.0+)
            family: Optional gene family for family-specific thresholds
        
        Returns:
            ACMG classification string: 'P', 'LP', 'VUS', 'LB', or 'B'
            (VUS-P collapsed into VUS as of 2025-12-11)
        """
        thr = self._get_family_thresholds(family)

        if score >= thr["P"]:
            return "P"
        elif score >= thr["LP"]:
            return "LP"
        elif score >= thr["VUS"]:
            return "VUS"
        elif score >= thr["LB"]:
            return "LB"
        else:
            return "B"
