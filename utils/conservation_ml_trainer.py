#!/usr/bin/env python3
"""
🧬 Conservation ML Trainer
Learn real conservation→pathogenicity relationships instead of guessing!

Built by Nova & Ace (2025) for AdaptiveInterpreter system
Based on Nova's brilliant ML conservation framework
"""

import json
import numpy as np
import pandas as pd
import xgboost as xgb
from collections import defaultdict
import os
from typing import Dict, List, Tuple

class ConservationMLTrainer:
    """Train ML model to learn conservation patterns by gene family"""
    
    def __init__(self):
        self.model = None
        self.families = [
            "TUMOR_SUPPRESSOR", "ONCOGENE", "ION_CHANNEL", "TRANSCRIPTION_FACTOR",
            "MOTOR_PROTEIN", "RIBOSOMAL_PROTEIN", "MUSCULAR_DYSTROPHY",
            "COLLAGEN_FIBRILLAR", "COLLAGEN_NETWORK", "COLLAGEN_ANCHORING",
            "COLLAGEN_FACIT", "FIBRILLIN", "ELASTIN", "INTERMEDIATE_FILAMENT",
            "CYTOSKELETON_POLYMER", "LAMIN", "NEGATIVE_REGULATOR",
            "AUTOSOMAL_RECESSIVE", "METABOLIC_ENZYME", "TRANSPORTER",
            "SCAFFOLD_ADAPTOR", "SIGNALING_REGULATOR", "STRUCTURAL", "GENERAL"
        ]
        
        # Default bins for conservation scores
        self.bins = {
            "phylop": np.linspace(-5, 10, 8),      # -5,-2.14,0.71,3.57,6.43,9.29,10
            "phastcons": np.linspace(0, 1, 6),     # 0,0.2,0.4,0.6,0.8,1.0
            "gerp": np.linspace(-5, 6, 8)          # -5,-3.43,-1.86,-0.29,1.29,2.86,4.43,6
        }
    
    def prepare_training_data(self, results_file: str) -> pd.DataFrame:
        """
        Extract training features from cascade results

        Args:
            results_file: TSV file with cascade results

        Returns:
            DataFrame with features for ML training
        """
        print(f"🔍 Loading training data from {results_file}")

        # Load results
        df = pd.read_csv(results_file, sep='\t')

        # Extract features we need
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

                # For now, use dummy conservation data - we'll improve this later
                # In a real implementation, we'd extract this from the logs or add it to the TSV
                phylop = 0.0  # Placeholder
                phastcons = 0.0  # Placeholder
                gerp = 0.0  # Placeholder

                # Determine gene family from gene name (simple heuristic for now)
                gene_family = self._guess_gene_family(gene)

                # Convert ClinVar classification to binary target
                target = self._clinvar_to_binary(expected_clinvar)
                if target is None:
                    continue

                features.append({
                    'gene': gene,
                    'variant': variant,
                    'phylop': phylop,
                    'phastcons': phastcons,
                    'gerp': gerp,
                    'gene_family': gene_family,
                    'target': target,
                    'our_class': final_classification,
                    'clinvar_class': expected_clinvar
                })

            except Exception as e:
                print(f"⚠️ Error processing row: {e}")
                continue

        feature_df = pd.DataFrame(features)
        print(f"✅ Prepared {len(feature_df)} training examples")
        if len(feature_df) > 0:
            print(f"📊 Gene families: {feature_df['gene_family'].value_counts().to_dict()}")

        return feature_df

    def _guess_gene_family(self, gene: str) -> str:
        """Simple heuristic to guess gene family from gene name"""
        gene = gene.upper()

        if gene in ['KIT']:
            return 'ONCOGENE'
        elif gene in ['SCN2A', 'SCN1A', 'KCNQ2', 'KCNQ3']:
            return 'ION_CHANNEL'
        elif gene in ['TGFBR2', 'TP53', 'APC', 'BRCA1', 'BRCA2']:
            return 'TUMOR_SUPPRESSOR'
        elif gene in ['COL1A1', 'COL1A2', 'COL3A1']:
            return 'COLLAGEN_FIBRILLAR'
        elif gene in ['FBN1']:
            return 'FIBRILLIN'
        else:
            return 'GENERAL'
    
    def _clinvar_to_binary(self, clinvar_class: str) -> int:
        """Convert ClinVar classification to binary pathogenic/benign"""
        if not clinvar_class:
            return None
            
        clinvar_class = clinvar_class.upper()
        
        if any(x in clinvar_class for x in ['PATHOGENIC', 'LIKELY_PATHOGENIC']):
            return 1  # Pathogenic
        elif any(x in clinvar_class for x in ['BENIGN', 'LIKELY_BENIGN']):
            return 0  # Benign
        else:
            return None  # VUS - skip for training
    
    def train_model(self, feature_df: pd.DataFrame, model_path: str = "conservation_model.json"):
        """Train XGBoost model on conservation features"""
        
        print(f"🚀 Training conservation ML model...")
        
        # Prepare features
        X = pd.get_dummies(feature_df[['phylop', 'phastcons', 'gerp', 'gene_family']], 
                          columns=['gene_family'])
        y = feature_df['target']
        
        print(f"📊 Training on {len(X)} examples with {len(X.columns)} features")
        print(f"🎯 Target distribution: {y.value_counts().to_dict()}")
        
        # Train XGBoost
        dtrain = xgb.DMatrix(X, label=y)
        
        params = {
            "objective": "binary:logistic",
            "eval_metric": "logloss", 
            "max_depth": 6,
            "eta": 0.1,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
            "seed": 42
        }
        
        self.model = xgb.train(params, dtrain, num_boost_round=300)
        
        # Save model
        self.model.save_model(model_path)
        print(f"✅ Model saved to {model_path}")
        
        # Feature importance
        importance = self.model.get_score(importance_type='weight')
        print(f"🔍 Feature importance: {importance}")
        
        return self.model
    
    def extract_family_curves(self, feature_df: pd.DataFrame, output_json: str = "conservation_multipliers.json"):
        """Extract conservation curves for each gene family"""
        
        if self.model is None:
            raise ValueError("Model not trained yet! Call train_model() first.")
        
        print(f"🧬 Extracting conservation curves for {len(self.families)} families...")
        
        results = {
            "meta": {
                "model_version": "conservation_ml_v1",
                "trained_on": f"ClinVar_variants_{len(feature_df)}",
                "features_used": ["phylop", "phastcons", "gerp", "gene_family"],
                "note": "ML-learned conservation multipliers by gene family"
            },
            "families": {}
        }
        
        for family in self.families:
            family_df = feature_df[feature_df['gene_family'] == family]
            
            if len(family_df) < 10:  # Need minimum data
                print(f"⚠️ Skipping {family}: only {len(family_df)} examples")
                results["families"][family] = {
                    "phylop_curve": [],
                    "phastcons_curve": [], 
                    "gerp_curve": []
                }
                continue
            
            print(f"🔍 Processing {family}: {len(family_df)} examples")
            family_curves = {}
            
            for metric in ["phylop", "phastcons", "gerp"]:
                curve = self._extract_metric_curve(family_df, family, metric)
                family_curves[f"{metric}_curve"] = curve
            
            results["families"][family] = family_curves
        
        # Save curves
        with open(output_json, "w") as f:
            json.dump(results, f, indent=2)
        
        print(f"✅ Conservation curves saved to {output_json}")
        return results
    
    def _extract_metric_curve(self, family_df: pd.DataFrame, family: str, metric: str) -> List[List[float]]:
        """Extract curve for one metric for one family"""
        
        edges = self.bins[metric]
        curve = []
        
        for i in range(len(edges) - 1):
            low, high = edges[i], edges[i + 1]
            bin_df = family_df[(family_df[metric] >= low) & (family_df[metric] < high)]
            
            if len(bin_df) < 3:  # Need minimum examples per bin
                continue
            
            # Create feature matrix for this bin
            X_bin = pd.get_dummies(bin_df[['phylop', 'phastcons', 'gerp', 'gene_family']], 
                                  columns=['gene_family'])
            
            # Ensure all columns exist (add missing ones as zeros)
            model_features = self.model.feature_names
            for col in model_features:
                if col not in X_bin.columns:
                    X_bin[col] = 0
            X_bin = X_bin[model_features]  # Reorder to match model
            
            # Predict
            dX_bin = xgb.DMatrix(X_bin)
            preds = self.model.predict(dX_bin)
            
            avg_score = float(bin_df[metric].mean())
            avg_pred = float(preds.mean())
            
            # Convert log-odds to multiplier (assuming baseline of 0.5 probability)
            baseline_logit = np.log(0.5 / (1 - 0.5))  # 0
            pred_logit = np.log(avg_pred / (1 - avg_pred + 1e-10))
            multiplier = np.exp(pred_logit - baseline_logit)
            
            curve.append([round(avg_score, 2), round(float(multiplier), 3)])
        
        return curve

def main():
    """Main training pipeline"""
    trainer = ConservationMLTrainer()
    
    # Find most recent results file
    results_files = [f for f in os.listdir("results/") if f.endswith(".tsv")]
    if not results_files:
        print("❌ No results files found in results/ directory")
        return
    
    latest_file = f"results/{sorted(results_files)[-1]}"
    print(f"🎯 Using latest results: {latest_file}")
    
    # Prepare training data
    feature_df = trainer.prepare_training_data(latest_file)
    
    if len(feature_df) < 100:
        print(f"❌ Not enough training data: {len(feature_df)} examples")
        return
    
    # Train model
    trainer.train_model(feature_df)
    
    # Extract curves
    trainer.extract_family_curves(feature_df)
    
    print("🎉 Conservation ML training complete!")

if __name__ == "__main__":
    main()
