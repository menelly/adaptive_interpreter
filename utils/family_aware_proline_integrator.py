#!/usr/bin/env python3
"""
🔥 FAMILY-AWARE PROLINE ML INTEGRATOR
Use family-specific proline models instead of one-size-fits-all!

Built by Ace (2025) for the Family-Aware ML Revolution
Replaces hardcoded proline logic with family-specific learned patterns
"""

import json
import numpy as np
import joblib
from typing import Dict, Optional, Tuple, List
from pathlib import Path
import re

class FamilyAwareProlineIntegrator:
    """Family-specific proline multiplier system"""
    
    def __init__(self, models_dir: str = "resources"):
        self.models_dir = Path(models_dir)
        self.models = {}  # Family -> trained model
        self.scalers = {}  # Family -> scaler
        self.family_mappings = {}
        
        self.load_models()
        self.load_family_mappings()
        
        print("🔥 Family-Aware Proline Integrator initialized!")
        print(f"📊 Loaded {len(self.models)} family-specific models")
    
    def load_models(self):
        """Load all family-specific proline models"""
        if not self.models_dir.exists():
            print(f"⚠️ Models directory not found: {self.models_dir}")
            return
        
        # Find all proline model files
        model_files = list(self.models_dir.glob("proline_model_*.joblib"))
        
        for model_file in model_files:
            # Extract family name from filename
            family = model_file.stem.replace("proline_model_", "").upper()
            scaler_file = self.models_dir / f"proline_scaler_{family.lower()}.joblib"
            
            try:
                self.models[family] = joblib.load(model_file)
                if scaler_file.exists():
                    self.scalers[family] = joblib.load(scaler_file)
                print(f"✅ Loaded {family} proline model")
            except Exception as e:
                print(f"⚠️ Error loading {family} model: {e}")
    
    def load_family_mappings(self):
        """Load gene -> family mappings"""
        mapping_file = self.models_dir / "proline_family_mappings.json"
        
        if mapping_file.exists():
            try:
                with open(mapping_file, 'r') as f:
                    self.family_mappings = json.load(f)
                print(f"✅ Loaded {len(self.family_mappings)} gene family mappings")
            except Exception as e:
                print(f"⚠️ Error loading family mappings: {e}")
    
    def guess_gene_family(self, gene: str) -> str:
        """Guess gene family with fallback logic"""
        gene = gene.upper()
        
        # Try exact mapping first
        if gene in self.family_mappings:
            return self.family_mappings[gene]
        
        # Fallback heuristics
        if gene.startswith('COL') and any(x in gene for x in ['1A1', '1A2', '3A1']):
            return 'COLLAGEN_FIBRILLAR'
        elif gene.startswith('SCN') or gene.startswith('KCNQ') or gene in ['CACNA1A', 'CACNA1C']:
            return 'ION_CHANNEL'
        elif gene in ['TP53', 'APC', 'BRCA1', 'BRCA2', 'TGFBR2']:
            return 'TUMOR_SUPPRESSOR'
        elif gene in ['KIT', 'EGFR', 'MET']:
            return 'ONCOGENE'
        elif gene in ['FBN1']:
            return 'FIBRILLIN'
        else:
            return 'GENERAL'
    
    def get_proline_multiplier(self, gene: str, variant: str, position: int = None) -> float:
        """
        Get family-aware proline multiplier
        
        Args:
            gene: Gene symbol
            variant: Variant string (e.g., "p.Pro123Leu")
            position: Amino acid position (optional, extracted from variant if not provided)
            
        Returns:
            Proline multiplier (1.0 = neutral, >1.0 = pathogenic boost, <1.0 = benign)
        """
        
        # Extract amino acid change - handle 3-letter codes too
        aa_match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', variant)
        if not aa_match:
            # Try single letter codes
            aa_match = re.search(r'p\.([A-Z])(\d+)([A-Z])', variant)
            if not aa_match:
                return 1.0  # No amino acid change detected
            ref_aa, pos_str, alt_aa = aa_match.groups()
        else:
            ref_aa_3, pos_str, alt_aa_3 = aa_match.groups()
            # Convert 3-letter to 1-letter codes
            aa_map = {'Pro': 'P', 'Leu': 'L', 'Ser': 'S', 'Arg': 'R', 'Ala': 'A', 'Gly': 'G'}
            ref_aa = aa_map.get(ref_aa_3, ref_aa_3[0])
            alt_aa = aa_map.get(alt_aa_3, alt_aa_3[0])
        position = position or int(pos_str)
        
        # Only process proline variants
        if ref_aa != 'P' and alt_aa != 'P':
            return 1.0  # No proline change
        
        # Determine gene family
        gene_family = self.guess_gene_family(gene)
        
        # Check if we have a trained model for this family
        if gene_family not in self.models:
            print(f"⚠️ No proline model for {gene_family}, using fallback")
            return self._fallback_proline_multiplier(ref_aa, alt_aa, gene_family)
        
        # Prepare features
        is_proline_loss = 1 if ref_aa == 'P' else 0
        is_proline_gain = 1 if alt_aa == 'P' else 0
        
        features = np.array([[position, is_proline_loss, is_proline_gain]])
        
        # Scale features if scaler available
        if gene_family in self.scalers:
            features = self.scalers[gene_family].transform(features)
        
        # Predict pathogenicity score
        try:
            pathogenicity_score = self.models[gene_family].predict(features)[0]
            
            # Convert pathogenicity score to multiplier
            # High pathogenicity = high multiplier
            multiplier = 1.0 + (pathogenicity_score * 2.0)  # Scale 0-1 to 1-3
            
            # Clamp to reasonable range
            multiplier = max(0.5, min(3.0, multiplier))
            
            print(f"🔥 {gene_family} proline ML: {variant} -> {multiplier:.2f}x (pathogenicity: {pathogenicity_score:.3f})")
            
            return multiplier
            
        except Exception as e:
            print(f"⚠️ Error predicting proline multiplier: {e}")
            return self._fallback_proline_multiplier(ref_aa, alt_aa, gene_family)
    
    def _fallback_proline_multiplier(self, ref_aa: str, alt_aa: str, gene_family: str) -> float:
        """Fallback proline multipliers when ML models unavailable"""
        
        # Family-specific fallback logic
        if gene_family == 'COLLAGEN_FIBRILLAR':
            if ref_aa == 'P':  # Proline loss in collagen
                return 2.5  # Very pathogenic
            else:  # Proline gain in collagen
                return 1.8  # Moderately pathogenic
        
        elif gene_family == 'ION_CHANNEL':
            if ref_aa == 'P':  # Proline loss in ion channel
                return 1.4  # Mildly pathogenic
            else:  # Proline gain in ion channel
                return 1.6  # Moderately pathogenic
        
        elif gene_family == 'FIBRILLIN':
            if ref_aa == 'P':  # Proline loss in fibrillin
                return 2.0  # Pathogenic
            else:  # Proline gain in fibrillin
                return 1.5  # Moderately pathogenic
        
        else:  # General fallback
            if ref_aa == 'P':  # Proline loss
                return 1.6  # Moderately pathogenic
            else:  # Proline gain
                return 1.3  # Mildly pathogenic
    
    def test_family_multipliers(self, test_variants: Dict[str, List[str]] = None):
        """Test proline multipliers across gene families"""
        
        if test_variants is None:
            test_variants = {
                'COL1A1': ['p.Pro123Leu', 'p.Leu456Pro'],
                'SCN5A': ['p.Pro789Ser', 'p.Ala321Pro'],
                'FBN1': ['p.Pro555Arg', 'p.Gly777Pro'],
                'TP53': ['p.Pro72Arg', 'p.Arg175Pro']
            }
        
        print("\n🔥 Testing Family-Aware Proline Multipliers")
        print("=" * 60)
        
        for gene, variants in test_variants.items():
            gene_family = self.guess_gene_family(gene)
            print(f"\n🧬 {gene} ({gene_family})")
            print("-" * 40)
            
            for variant in variants:
                multiplier = self.get_proline_multiplier(gene, variant)
                
                if multiplier > 1.5:
                    effect = "HIGH IMPACT"
                elif multiplier > 1.2:
                    effect = "MODERATE"
                elif multiplier > 0.9:
                    effect = "MILD"
                else:
                    effect = "PROTECTIVE"
                
                print(f"  {variant:<15} -> {multiplier:.2f}x ({effect})")

def main():
    """Test the family-aware proline integrator"""
    integrator = FamilyAwareProlineIntegrator()
    
    # Test with example variants
    integrator.test_family_multipliers()

if __name__ == "__main__":
    main()
