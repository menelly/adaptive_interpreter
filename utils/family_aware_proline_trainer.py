#!/usr/bin/env python3
"""
🔥 FAMILY-AWARE PROLINE ML TRAINER
Learn proline patterns by gene family - no more one-size-fits-all!

Built by Ace (2025) for the Family-Aware ML Revolution
Based on Ren's insight: "Proline in collagens ≠ proline in ion channels!"
"""

import json
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
import joblib
from typing import Dict, List, Tuple, Any
from pathlib import Path

class FamilyAwareProlineTrainer:
    """Train family-specific proline multiplier models"""
    
    def __init__(self, category_keywords_path: str = "category_keywords.json"):
        self.category_keywords_path = category_keywords_path
        self.family_mappings = self.load_family_mappings()
        self.models = {}  # Family -> trained model
        self.scalers = {}  # Family -> scaler
        
        print("🔥 Family-Aware Proline Trainer initialized!")
        print(f"📊 Loaded {len(self.family_mappings)} gene family mappings")
    
    def load_family_mappings(self) -> Dict[str, str]:
        """Load gene -> family mappings from category_keywords.json"""
        try:
            with open(self.category_keywords_path, 'r') as f:
                data = json.load(f)
            
            mappings = {}
            for family, family_data in data.get('categories', {}).items():
                genes = family_data.get('genes', [])
                for gene in genes:
                    mappings[gene.upper()] = family
            
            return mappings
        except Exception as e:
            print(f"⚠️ Error loading family mappings: {e}")
            return {}
    
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
    
    def prepare_training_data(self, results_file: str) -> pd.DataFrame:
        """Extract proline training features from cascade results"""
        print(f"🔍 Loading training data from {results_file}")
        
        # Load results
        df = pd.read_csv(results_file, sep='\t')
        
        features = []
        
        for _, row in df.iterrows():
            try:
                gene = row.get('gene', '')
                variant = row.get('variant', '')
                final_classification = row.get('final_classification', '')
                expected_clinvar = row.get('expected_clinvar', '')
                
                # Skip if missing key data
                if not all([gene, variant, final_classification, expected_clinvar]):
                    continue
                
                # Extract amino acid change
                import re
                aa_match = re.search(r'p\.([A-Z])(\d+)([A-Z])', variant)
                if not aa_match:
                    continue
                
                ref_aa, position, alt_aa = aa_match.groups()
                position = int(position)
                
                # Only process proline variants
                if ref_aa != 'P' and alt_aa != 'P':
                    continue
                
                # Determine gene family
                gene_family = self.guess_gene_family(gene)
                
                # Convert ClinVar to pathogenicity score
                pathogenicity = self._clinvar_to_score(expected_clinvar)
                if pathogenicity is None:
                    continue
                
                # Extract features
                features.append({
                    'gene': gene,
                    'gene_family': gene_family,
                    'position': position,
                    'ref_aa': ref_aa,
                    'alt_aa': alt_aa,
                    'is_proline_loss': 1 if ref_aa == 'P' else 0,
                    'is_proline_gain': 1 if alt_aa == 'P' else 0,
                    'pathogenicity_score': pathogenicity,
                    'our_classification': final_classification,
                    'clinvar_classification': expected_clinvar
                })
                
            except Exception as e:
                print(f"⚠️ Error processing row: {e}")
                continue
        
        feature_df = pd.DataFrame(features)
        print(f"✅ Prepared {len(feature_df)} proline training examples")
        
        if len(feature_df) > 0:
            print(f"📊 Gene families: {feature_df['gene_family'].value_counts().to_dict()}")
            print(f"🧬 Proline changes: Loss={feature_df['is_proline_loss'].sum()}, Gain={feature_df['is_proline_gain'].sum()}")
        
        return feature_df
    
    def _clinvar_to_score(self, clinvar_class: str) -> float:
        """Convert ClinVar classification to pathogenicity score (0-1)"""
        if not clinvar_class:
            return None
            
        clinvar_class = clinvar_class.upper()
        
        if 'PATHOGENIC' in clinvar_class and 'LIKELY' not in clinvar_class:
            return 1.0  # Pathogenic
        elif 'LIKELY_PATHOGENIC' in clinvar_class or 'LIKELY PATHOGENIC' in clinvar_class:
            return 0.8  # Likely Pathogenic
        elif 'BENIGN' in clinvar_class and 'LIKELY' not in clinvar_class:
            return 0.0  # Benign
        elif 'LIKELY_BENIGN' in clinvar_class or 'LIKELY BENIGN' in clinvar_class:
            return 0.2  # Likely Benign
        else:
            return None  # VUS - skip for training
    
    def train_family_models(self, feature_df: pd.DataFrame):
        """Train separate models for each gene family"""
        
        families_with_data = feature_df['gene_family'].value_counts()
        print(f"🚀 Training family-specific proline models...")
        
        for family in families_with_data.index:
            family_data = feature_df[feature_df['gene_family'] == family]
            
            if len(family_data) < 10:  # Need minimum data
                print(f"⚠️ Skipping {family}: only {len(family_data)} examples")
                continue
            
            print(f"🔍 Training {family}: {len(family_data)} examples")
            
            # Prepare features
            X = family_data[['position', 'is_proline_loss', 'is_proline_gain']].values
            y = family_data['pathogenicity_score'].values
            
            # Split data
            if len(family_data) >= 20:
                X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
            else:
                X_train, X_test, y_train, y_test = X, X, y, y  # Use all data for small sets
            
            # Scale features
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)
            
            # Train model
            model = RandomForestRegressor(n_estimators=100, random_state=42)
            model.fit(X_train_scaled, y_train)
            
            # Evaluate
            y_pred = model.predict(X_test_scaled)
            mse = mean_squared_error(y_test, y_pred)
            r2 = r2_score(y_test, y_pred)
            
            print(f"   📊 {family}: MSE={mse:.3f}, R²={r2:.3f}")
            
            # Store model and scaler
            self.models[family] = model
            self.scalers[family] = scaler
        
        print(f"✅ Trained {len(self.models)} family-specific models!")
    
    def save_models(self, output_dir: str = "resources"):
        """Save trained models and scalers"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        for family in self.models:
            model_file = output_path / f"proline_model_{family.lower()}.joblib"
            scaler_file = output_path / f"proline_scaler_{family.lower()}.joblib"
            
            joblib.dump(self.models[family], model_file)
            joblib.dump(self.scalers[family], scaler_file)
            
            print(f"💾 Saved {family} model: {model_file}")
        
        # Save family mappings
        mapping_file = output_path / "proline_family_mappings.json"
        with open(mapping_file, 'w') as f:
            json.dump(self.family_mappings, f, indent=2)
        
        print(f"✅ All family-aware proline models saved to {output_dir}/")

def main():
    """Main training pipeline"""
    trainer = FamilyAwareProlineTrainer()
    
    # Find most recent results file
    import os
    results_files = [f for f in os.listdir("tests/results/") if f.endswith(".tsv")]
    if not results_files:
        print("❌ No results files found in tests/results/ directory")
        return
    
    latest_file = f"tests/results/{sorted(results_files)[-1]}"
    print(f"🎯 Using latest results: {latest_file}")
    
    # Prepare training data
    feature_df = trainer.prepare_training_data(latest_file)
    
    if len(feature_df) < 20:
        print(f"❌ Not enough proline training data: {len(feature_df)} examples")
        return
    
    # Train family-specific models
    trainer.train_family_models(feature_df)
    
    # Save models
    trainer.save_models()
    
    print("🎉 Family-aware proline ML training complete!")

if __name__ == "__main__":
    main()
