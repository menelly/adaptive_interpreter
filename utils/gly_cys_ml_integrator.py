#!/usr/bin/env python3
"""
🔥 REVOLUTIONARY GLYCINE & CYSTEINE ML SYSTEM INTEGRATOR
Replaces hardcoded Gly/Cys multipliers with data-driven intelligence!

🧬 REVOLUTIONARY CONCEPT: Stop guessing Gly/Cys penalties, start learning from ClinVar data!
⚡ IMPLEMENTATION: Ace - built the ML training pipeline and integration following Proline success
🎯 COLLABORATION: Revolutionary genomics through biological intelligence + machine learning

This system replaces hardcoded Gly/Cys penalties with trained machine learning models
that understand biological context and learn from real pathogenicity data.

INTEGRATION POINTS:
- Cascade Analyzer: Replace hardcoded Gly/Cys penalties
- LOF Analyzer: Replace hardcoded Gly/Cys scoring
- GOF Analyzer: Replace hardcoded Gly/Cys boosts  
- DN Analyzer: Replace Gly/Cys-specific scoring

REVOLUTIONARY BECAUSE: Context-aware, data-driven, biologically intelligent!
"""

import numpy as np
import joblib
from typing import Dict, Any, Optional, List
import os
from pathlib import Path

from gly_cys_context import GlyCysContextAnalyzer
from universal_protein_annotator import UniversalProteinAnnotator
from proline_multiplier_mapper import map_prob_to_multiplier


class GlyCysMLIntegrator:
    """ML-powered Gly/Cys multiplier system to replace hardcoded approaches"""
    
    def __init__(self, 
                 gly_model_path: str = "gly_ml_model.joblib",
                 gly_scaler_path: str = "gly_scaler.joblib",
                 cys_model_path: str = "cys_ml_model.joblib",
                 cys_scaler_path: str = "cys_scaler.joblib",
                 alphafold_path: str = "./alphafold_structures/"):
        
        self.gly_model_path = gly_model_path
        self.gly_scaler_path = gly_scaler_path
        self.cys_model_path = cys_model_path
        self.cys_scaler_path = cys_scaler_path
        self.alphafold_path = alphafold_path
        
        # Initialize context analyzer
        self.context_analyzer = GlyCysContextAnalyzer(alphafold_path)
        self.annotator = None  # Skip for now
        
        # Try to load trained models
        self.gly_model = None
        self.gly_scaler = None
        self.cys_model = None
        self.cys_scaler = None
        self.load_models()
        
        # Fallback parameters if models not available
        self.fallback_enabled = True
        self.fallback_multipliers = {
            # Glycine fallbacks
            'collagen_gly_loss': 2.5,      # Collagen Gly loss = very pathogenic
            'collagen_gly_gain': 1.8,      # Collagen Gly gain = moderately pathogenic
            'ion_channel_gly_loss': 1.4,   # Ion channel Gly loss = context-dependent
            'ion_channel_gly_gain': 1.2,   # Ion channel Gly gain = mild impact
            'general_gly_loss': 1.3,       # General Gly loss = moderate
            'general_gly_gain': 1.1,       # General Gly gain = mild
            
            # Cysteine fallbacks
            'disulfide_cys_loss': 2.2,     # Disulfide Cys loss = very pathogenic
            'disulfide_cys_gain': 1.6,     # Disulfide Cys gain = problematic
            'metal_coord_cys_loss': 2.0,   # Metal coordination Cys loss = pathogenic
            'catalytic_cys_loss': 1.9,     # Catalytic Cys loss = pathogenic
            'general_cys_loss': 1.5,       # General Cys loss = moderate
            'general_cys_gain': 1.3        # General Cys gain = mild-moderate
        }
        
        print("🧬 Gly/Cys ML Integrator initialized!")
        if self.gly_model is not None or self.cys_model is not None:
            print("🔥 ML models loaded - ready for data-driven analysis!")
        else:
            print("⚠️  ML models not found - using intelligent fallback system")
    
    def load_models(self) -> None:
        """Load trained ML models and scalers"""
        
        try:
            # Load Glycine model
            if os.path.exists(self.gly_model_path) and os.path.exists(self.gly_scaler_path):
                self.gly_model = joblib.load(self.gly_model_path)
                self.gly_scaler = joblib.load(self.gly_scaler_path)
                print("✅ Glycine ML model loaded successfully")
            else:
                print(f"⚠️  Glycine model files not found: {self.gly_model_path}, {self.gly_scaler_path}")
            
            # Load Cysteine model
            if os.path.exists(self.cys_model_path) and os.path.exists(self.cys_scaler_path):
                self.cys_model = joblib.load(self.cys_model_path)
                self.cys_scaler = joblib.load(self.cys_scaler_path)
                print("✅ Cysteine ML model loaded successfully")
            else:
                print(f"⚠️  Cysteine model files not found: {self.cys_model_path}, {self.cys_scaler_path}")
                
        except Exception as e:
            print(f"❌ Error loading ML models: {e}")
            self.gly_model = None
            self.gly_scaler = None
            self.cys_model = None
            self.cys_scaler = None
    
    def build_feature_vector(self, 
                           gene: str,
                           position: int,
                           ref_aa: str,
                           alt_aa: str,
                           gnomad_freq: float = 0.0) -> Optional[np.ndarray]:
        """Build feature vector for ML prediction"""
        
        try:
            # Get biological context
            context = self.context_analyzer.get_context_features(
                gene=gene,
                position=position,
                ref_aa=ref_aa,
                alt_aa=alt_aa
            )
            
            if 'error' in context:
                return None
            
            # Build feature vector based on amino acid type
            amino_acid = context.get('amino_acid', '')
            
            if amino_acid == 'GLYCINE':
                return self._build_glycine_features(context)
            elif amino_acid == 'CYSTEINE':
                return self._build_cysteine_features(context)
            else:
                return None
                
        except Exception as e:
            print(f"❌ Feature building failed for {gene} p.{ref_aa}{position}{alt_aa}: {e}")
            return None
    
    def _build_glycine_features(self, context: Dict[str, Any]) -> np.ndarray:
        """Build feature vector for Glycine variants (matches trainer)"""
        
        features = []
        
        # Protein family encoding (one-hot)
        protein_family = context.get('protein_family', 'OTHER')
        features.extend([
            1.0 if protein_family == 'COLLAGEN' else 0.0,
            1.0 if protein_family == 'ION_CHANNEL' else 0.0,
            1.0 if protein_family == 'FIBRILLIN' else 0.0,
            1.0 if protein_family == 'OTHER' else 0.0
        ])
        
        # Substitution type
        features.append(1.0 if context.get('substitution_type') == 'LOSS' else 0.0)
        
        # Conservation score
        features.append(float(context.get('conservation_score', 0.5)))
        
        # Collagen-specific features
        features.extend([
            1.0 if context.get('collagen_gxy_pattern', False) else 0.0,
            float(context.get('collagen_gxy_confidence', 0.0)),
            1.0 if context.get('triple_helix_region', False) else 0.0,
            1.0 if context.get('collagen_cleavage_site', False) else 0.0
        ])
        
        # Ion channel-specific features
        features.extend([
            1.0 if context.get('channel_gate_region', False) else 0.0,
            1.0 if context.get('transmembrane_domain', False) else 0.0,
            1.0 if context.get('channel_selectivity_filter', False) else 0.0,
            1.0 if context.get('channel_linker_region', False) else 0.0
        ])
        
        # Fibrillin-specific features
        features.extend([
            1.0 if context.get('egf_like_domain', False) else 0.0,
            1.0 if context.get('calcium_binding_egf', False) else 0.0,
            1.0 if context.get('fibrillin_unique_domain', False) else 0.0
        ])
        
        # General structural features
        features.extend([
            1.0 if context.get('active_site_proximity', False) else 0.0,
            1.0 if context.get('binding_site_proximity', False) else 0.0,
            1.0 if context.get('flexible_region', False) else 0.0,
            1.0 if context.get('tight_packing_region', False) else 0.0
        ])
        
        return np.array(features)
    
    def _build_cysteine_features(self, context: Dict[str, Any]) -> np.ndarray:
        """Build feature vector for Cysteine variants (matches trainer)"""
        
        features = []
        
        # Protein family encoding (one-hot)
        protein_family = context.get('protein_family', 'OTHER')
        features.extend([
            1.0 if protein_family == 'COLLAGEN' else 0.0,
            1.0 if protein_family == 'ION_CHANNEL' else 0.0,
            1.0 if protein_family == 'FIBRILLIN' else 0.0,
            1.0 if protein_family == 'OTHER' else 0.0
        ])
        
        # Substitution type
        features.append(1.0 if context.get('substitution_type') == 'LOSS' else 0.0)
        
        # Conservation score
        features.append(float(context.get('conservation_score', 0.5)))
        
        # Disulfide bond features
        features.extend([
            1.0 if context.get('disulfide_bond_predicted', False) else 0.0,
            1.0 if context.get('structural_disulfide', False) else 0.0,
            1.0 if context.get('extracellular_cysteine', False) else 0.0
        ])
        
        # Metal coordination features
        features.extend([
            1.0 if context.get('metal_coordination_site', False) else 0.0,
            1.0 if context.get('catalytic_site_proximity', False) else 0.0
        ])
        
        # Protein-specific features
        features.extend([
            1.0 if context.get('rare_collagen_cysteine', False) else 0.0,
            1.0 if context.get('egf_domain_cysteine', False) else 0.0,
            1.0 if context.get('calcium_binding_cysteine', False) else 0.0,
            1.0 if context.get('channel_structure_critical', False) else 0.0
        ])
        
        # Risk assessment features
        features.extend([
            float(context.get('free_cysteine_risk', 0.0)),
            1.0 if context.get('redox_sensitive_region', False) else 0.0,
            1.0 if context.get('binding_site_proximity', False) else 0.0
        ])
        
        return np.array(features)
    
    def get_gly_cys_multiplier(self, 
                              gene: str,
                              position: int, 
                              ref_aa: str, 
                              alt_aa: str,
                              gnomad_freq: float = 0.0,
                              method: str = "sigmoid",
                              **map_kwargs) -> float:
        """
        Get Gly/Cys multiplier using ML model or intelligent fallback
        
        This is the main interface that replaces hardcoded Gly/Cys penalties
        """
        
        # Only apply to Gly/Cys substitutions
        if ref_aa.upper() not in ['G', 'C'] and alt_aa.upper() not in ['G', 'C']:
            return 1.0
        
        # Determine amino acid type
        is_glycine = ref_aa.upper() == 'G' or alt_aa.upper() == 'G'
        is_cysteine = ref_aa.upper() == 'C' or alt_aa.upper() == 'C'
        
        # Try ML prediction first
        if is_glycine and self.gly_model is not None and self.gly_scaler is not None:
            try:
                features = self.build_feature_vector(gene, position, ref_aa, alt_aa, gnomad_freq)
                if features is not None:
                    # Scale features and predict
                    features_scaled = self.gly_scaler.transform([features])
                    prob = self.gly_model.predict_proba(features_scaled)[0][1]  # P(pathogenic)
                    
                    # Map probability to multiplier
                    multiplier = map_prob_to_multiplier(prob, method=method, **map_kwargs)
                    
                    print(f"🧬 ML Glycine prediction: {gene} p.{ref_aa}{position}{alt_aa} -> P={prob:.3f}, mult={multiplier:.3f}")
                    return multiplier
                    
            except Exception as e:
                print(f"⚠️  ML Glycine prediction failed for {gene} p.{ref_aa}{position}{alt_aa}: {e}")
        
        elif is_cysteine and self.cys_model is not None and self.cys_scaler is not None:
            try:
                features = self.build_feature_vector(gene, position, ref_aa, alt_aa, gnomad_freq)
                if features is not None:
                    # Scale features and predict
                    features_scaled = self.cys_scaler.transform([features])
                    prob = self.cys_model.predict_proba(features_scaled)[0][1]  # P(pathogenic)
                    
                    # Map probability to multiplier
                    multiplier = map_prob_to_multiplier(prob, method=method, **map_kwargs)
                    
                    print(f"🧬 ML Cysteine prediction: {gene} p.{ref_aa}{position}{alt_aa} -> P={prob:.3f}, mult={multiplier:.3f}")
                    return multiplier
                    
            except Exception as e:
                print(f"⚠️  ML Cysteine prediction failed for {gene} p.{ref_aa}{position}{alt_aa}: {e}")
        
        # Fall back to intelligent biological reasoning
        return self._get_intelligent_fallback_multiplier(gene, position, ref_aa, alt_aa)
    
    def _get_intelligent_fallback_multiplier(self, 
                                           gene: str,
                                           position: int,
                                           ref_aa: str,
                                           alt_aa: str) -> float:
        """Intelligent fallback using biological context (not random guessing!)"""
        
        try:
            # Get biological context for intelligent fallback
            context = self.context_analyzer.get_context_features(
                gene=gene, position=position, ref_aa=ref_aa, alt_aa=alt_aa
            )
            
            if 'error' in context:
                return 1.0
            
            amino_acid = context.get('amino_acid', '')
            protein_family = context.get('protein_family', 'OTHER')
            substitution_type = context.get('substitution_type', 'UNKNOWN')
            
            if amino_acid == 'GLYCINE':
                return self._get_glycine_fallback(context, protein_family, substitution_type)
            elif amino_acid == 'CYSTEINE':
                return self._get_cysteine_fallback(context, protein_family, substitution_type)
            else:
                return 1.0
                
        except Exception as e:
            print(f"⚠️  Fallback analysis failed for {gene} p.{ref_aa}{position}{alt_aa}: {e}")
            return 1.0
    
    def _get_glycine_fallback(self, context: Dict[str, Any], protein_family: str, substitution_type: str) -> float:
        """Intelligent Glycine fallback based on biological context"""
        
        # Collagen glycines - ALWAYS critical in Gly-X-Y pattern
        if protein_family == 'COLLAGEN':
            if context.get('collagen_gxy_pattern', False):
                multiplier = self.fallback_multipliers['collagen_gly_loss'] if substitution_type == 'LOSS' else self.fallback_multipliers['collagen_gly_gain']
                print(f"🧬 Intelligent fallback: Collagen Gly-X-Y pattern -> {multiplier:.3f}")
                return multiplier
        
        # Ion channel glycines - context-dependent
        elif protein_family == 'ION_CHANNEL':
            if context.get('channel_gate_region', False) or context.get('channel_selectivity_filter', False):
                multiplier = self.fallback_multipliers['ion_channel_gly_loss'] if substitution_type == 'LOSS' else self.fallback_multipliers['ion_channel_gly_gain']
                print(f"🧬 Intelligent fallback: Ion channel critical Gly -> {multiplier:.3f}")
                return multiplier
        
        # General glycine
        multiplier = self.fallback_multipliers['general_gly_loss'] if substitution_type == 'LOSS' else self.fallback_multipliers['general_gly_gain']
        print(f"🧬 Intelligent fallback: General Gly -> {multiplier:.3f}")
        return multiplier
    
    def _get_cysteine_fallback(self, context: Dict[str, Any], protein_family: str, substitution_type: str) -> float:
        """Intelligent Cysteine fallback based on biological context"""
        
        # Disulfide bond cysteines - highly critical
        if context.get('disulfide_bond_predicted', False):
            multiplier = self.fallback_multipliers['disulfide_cys_loss'] if substitution_type == 'LOSS' else self.fallback_multipliers['disulfide_cys_gain']
            print(f"🧬 Intelligent fallback: Disulfide bond Cys -> {multiplier:.3f}")
            return multiplier
        
        # Metal coordination cysteines
        if context.get('metal_coordination_site', False):
            multiplier = self.fallback_multipliers['metal_coord_cys_loss']
            print(f"🧬 Intelligent fallback: Metal coordination Cys -> {multiplier:.3f}")
            return multiplier
        
        # Catalytic cysteines
        if context.get('catalytic_site_proximity', False):
            multiplier = self.fallback_multipliers['catalytic_cys_loss']
            print(f"🧬 Intelligent fallback: Catalytic Cys -> {multiplier:.3f}")
            return multiplier
        
        # General cysteine
        multiplier = self.fallback_multipliers['general_cys_loss'] if substitution_type == 'LOSS' else self.fallback_multipliers['general_cys_gain']
        print(f"🧬 Intelligent fallback: General Cys -> {multiplier:.3f}")
        return multiplier


# Convenience function for easy integration
def get_ml_gly_cys_multiplier(gene: str, variant: str, gnomad_freq: float = 0.0) -> float:
    """
    Convenience function to get ML Gly/Cys multiplier
    
    Args:
        gene: Gene symbol (e.g., 'COL1A1', 'FBN1')
        variant: Protein variant (e.g., 'p.G893A', 'p.C628Y')
        gnomad_freq: Population frequency (default 0.0)
        
    Returns:
        float: Multiplier for variant scoring (0.5 to 2.5)
    """
    
    # Parse variant
    import re
    match = re.match(r'p\.([A-Z])(\d+)([A-Z])', variant)
    if not match:
        return 1.0
        
    ref_aa, pos, alt_aa = match.groups()
    
    # Initialize integrator (cached)
    if not hasattr(get_ml_gly_cys_multiplier, '_integrator'):
        get_ml_gly_cys_multiplier._integrator = GlyCysMLIntegrator()
    
    return get_ml_gly_cys_multiplier._integrator.get_gly_cys_multiplier(
        gene, int(pos), ref_aa, alt_aa, gnomad_freq
    )


def test_gly_cys_ml_integration():
    """Test the Gly/Cys ML integration system"""
    
    integrator = GlyCysMLIntegrator()
    
    print("\n🧬 Testing Gly/Cys ML Integration:")
    print("=" * 60)
    
    # Test collagen glycine (should be highly pathogenic)
    print("\n1. COLLAGEN GLYCINE TEST:")
    mult1 = integrator.get_gly_cys_multiplier('COL1A1', 893, 'G', 'A')
    print(f"COL1A1 p.G893A multiplier: {mult1}")
    
    # Test fibrillin cysteine (should be disulfide critical)
    print("\n2. FIBRILLIN CYSTEINE TEST:")
    mult2 = integrator.get_gly_cys_multiplier('FBN1', 628, 'C', 'Y')
    print(f"FBN1 p.C628Y multiplier: {mult2}")
    
    # Test ion channel glycine (should be context-dependent)
    print("\n3. ION CHANNEL GLYCINE TEST:")
    mult3 = integrator.get_gly_cys_multiplier('SCN1A', 58, 'G', 'R')
    print(f"SCN1A p.G58R multiplier: {mult3}")
    
    # Test convenience function
    print("\n4. CONVENIENCE FUNCTION TEST:")
    mult4 = get_ml_gly_cys_multiplier('RYR1', 'p.R614C')
    print(f"RYR1 p.R614C multiplier: {mult4}")


if __name__ == "__main__":
    test_gly_cys_ml_integration()
