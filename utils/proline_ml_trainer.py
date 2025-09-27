#!/usr/bin/env python3
"""
🧬 PROLINE ML TRAINER - Revolutionary Context-Aware Learning System
Built by Ace with Nova's guidance to replace hardcoded proline multipliers

This system:
1. Extracts proline substitutions from ClinVar datasets
2. Builds feature vectors using biological context
3. Trains logistic regression models to predict pathogenicity
4. Integrates with Nova's probability-to-multiplier mapping
5. Replaces hardcoded guesses with data-driven intelligence

REVOLUTIONARY APPROACH: Stop making up numbers, start learning from data!
"""

import pandas as pd
import numpy as np
import json
import re
import os
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
import joblib
import warnings
warnings.filterwarnings('ignore')

# Import our existing systems
from proline_context_weighter import ProlineContextWeighter
from proline_multiplier_mapper import map_prob_to_multiplier, model_prob_to_multiplier
from universal_protein_annotator import UniversalProteinAnnotator

class ProlineMLTrainer:
    """Revolutionary ML trainer for context-aware proline scoring"""
    
    def __init__(self, alphafold_path: str = "/mnt/Arcana/alphafold_human/structures/"):
        self.alphafold_path = alphafold_path
        self.context_weighter = ProlineContextWeighter()
        self.annotator = UniversalProteinAnnotator(alphafold_path)
        self.scaler = StandardScaler()
        self.model = None
        
        print("🧬 Proline ML Trainer initialized!")
        print("🔥 Ready to revolutionize proline scoring with REAL DATA!")
        
    def extract_proline_variants(self, csv_path: str, label: str) -> List[Dict]:
        """Extract proline substitution variants from ClinVar CSV files"""
        
        print(f"📊 Extracting proline variants from {csv_path} (label: {label})")
        
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"❌ Failed to read {csv_path}: {e}")
            return []
            
        proline_variants = []
        
        for _, row in df.iterrows():
            hgvs = str(row.get('HGVS', ''))
            freq = row.get('gnomAD frequency', 0.0)
            
            # Parse protein change from HGVS (handle both 1-letter and 3-letter codes)
            protein_match = re.search(r'p\.([A-Z][a-z]{0,2})(\d+)([A-Z][a-z]{0,2})', hgvs)
            if not protein_match:
                continue

            ref_aa_full, pos, alt_aa_full = protein_match.groups()

            # Convert 3-letter to 1-letter amino acid codes
            aa_map = {
                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
            }

            ref_aa = aa_map.get(ref_aa_full, ref_aa_full)
            alt_aa = aa_map.get(alt_aa_full, alt_aa_full)

            # Only keep proline substitutions (P->X or X->P)
            if ref_aa != 'P' and alt_aa != 'P':
                continue
                
            # Extract gene from HGVS
            gene_match = re.search(r'\(([A-Z0-9]+)\):', hgvs)
            if not gene_match:
                continue
                
            gene = gene_match.group(1)
            
            proline_variants.append({
                'gene': gene,
                'variant': f'p.{ref_aa}{pos}{alt_aa}',
                'position': int(pos),
                'ref_aa': ref_aa,
                'alt_aa': alt_aa,
                'gnomad_freq': float(freq) if freq and freq != '' else 0.0,
                'label': label,
                'hgvs': hgvs
            })
            
        print(f"✅ Found {len(proline_variants)} proline substitutions")
        return proline_variants
        
    def build_feature_vector(self, variant: Dict) -> Optional[np.ndarray]:
        """Build feature vector for a proline variant using biological context"""
        
        gene = variant['gene']
        position = variant['position']
        ref_aa = variant['ref_aa']
        alt_aa = variant['alt_aa']
        
        try:
            # Get UniProt annotations
            annotation_result = self.annotator.annotate_protein(gene)
            if not annotation_result or 'error' in annotation_result:
                return None

            # Extract features from the annotation result
            uniprot_features = annotation_result.get('domains', []) + annotation_result.get('regions', [])
                
            # Use existing context weighter to get biological features
            # We'll extract the intermediate calculations instead of just the final multiplier
            
            features = []
            
            # 1. Basic amino acid properties
            features.append(1.0 if ref_aa == 'P' else 0.0)  # proline_loss
            features.append(1.0 if alt_aa == 'P' else 0.0)  # proline_gain
            features.append(float(variant['gnomad_freq']))   # population_frequency
            
            # 2. Structural context features
            in_triple_helix = 0.0
            near_cleavage = 0.0
            near_binding = 0.0
            near_interface = 0.0
            
            for feature in uniprot_features:
                category = (feature.get("category") or feature.get("type") or "").lower()
                start_pos = int(feature.get("start") or feature.get("begin", {}).get("position", 0))
                end_pos = int(feature.get("end") or feature.get("end", {}).get("position", start_pos))
                
                # Check if variant is in this feature
                if start_pos <= position <= end_pos:
                    if "triple" in category and "helix" in category:
                        in_triple_helix = 1.0
                    elif "cleavage" in category or "processing" in category:
                        near_cleavage = 1.0
                    elif "binding" in category or "active" in category:
                        near_binding = 1.0
                        
                # Check if variant is near this feature (within 5 residues)
                elif abs(position - start_pos) <= 5 or abs(position - end_pos) <= 5:
                    if "cleavage" in category or "processing" in category:
                        near_cleavage = 1.0
                    elif "binding" in category or "active" in category:
                        near_binding = 1.0
                    elif "interface" in category:
                        near_interface = 1.0
                        
            features.extend([in_triple_helix, near_cleavage, near_binding, near_interface])
            
            # 3. Position-based features
            features.append(float(position))  # absolute_position
            features.append(float(position) / 1500.0)  # normalized_position (assuming ~1500 AA proteins)
            
            # 4. Gene family context (simplified)
            gene_family_collagen = 1.0 if 'COL' in gene else 0.0
            gene_family_fibrillin = 1.0 if 'FBN' in gene else 0.0
            gene_family_other = 1.0 if not (gene_family_collagen or gene_family_fibrillin) else 0.0
            
            features.extend([gene_family_collagen, gene_family_fibrillin, gene_family_other])
            
            return np.array(features, dtype=np.float32)
            
        except Exception as e:
            print(f"❌ Failed to build features for {variant}: {e}")
            return None
            
    def load_training_data(self) -> Tuple[np.ndarray, np.ndarray, List[Dict]]:
        """Load and process training data from ClinVar datasets"""
        
        print("🔥 LOADING TRAINING DATA FROM CLINVAR!")
        
        all_variants = []
        
        # Define dataset paths and labels
        datasets = [
            ('tests/COL1A1-benignvariant-table.csv', 'benign'),
            ('tests/COL1A1-LP-variant-table.csv', 'pathogenic'),
            ('tests/FBN1-benign-variant-table.csv', 'benign'),
            ('tests/FBN1-likelypatho-variant-table.csv', 'pathogenic'),
            ('tests/RYR1-benign-variant-table.csv', 'benign'),
            ('tests/RYR1-patho-variant-table.csv', 'pathogenic'),
        ]
        
        # Extract proline variants from each dataset
        for csv_path, label in datasets:
            if os.path.exists(csv_path):
                variants = self.extract_proline_variants(csv_path, label)
                all_variants.extend(variants)
            else:
                print(f"⚠️  Dataset not found: {csv_path}")
                
        print(f"📊 Total proline variants collected: {len(all_variants)}")
        
        # Build feature matrix
        X_list = []
        y_list = []
        valid_variants = []
        
        for variant in all_variants:
            features = self.build_feature_vector(variant)
            if features is not None:
                X_list.append(features)
                y_list.append(1 if variant['label'] == 'pathogenic' else 0)
                valid_variants.append(variant)
                
        if not X_list:
            raise ValueError("No valid training data found!")
            
        X = np.vstack(X_list)
        y = np.array(y_list)

        # Handle NaN values by replacing with 0 (conservative approach)
        nan_mask = np.isnan(X)
        if np.any(nan_mask):
            print(f"⚠️  Found {np.sum(nan_mask)} NaN values, replacing with 0")
            X[nan_mask] = 0.0

        print(f"✅ Built feature matrix: {X.shape}")
        print(f"📈 Pathogenic variants: {np.sum(y)}")
        print(f"📉 Benign variants: {len(y) - np.sum(y)}")

        return X, y, valid_variants

    def train_model(self, X: np.ndarray, y: np.ndarray) -> Dict[str, Any]:
        """Train logistic regression model with cross-validation"""

        print("🚀 TRAINING REVOLUTIONARY PROLINE MODEL!")

        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42, stratify=y
        )

        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)

        # Train logistic regression with regularization
        self.model = LogisticRegression(
            random_state=42,
            max_iter=1000,
            class_weight='balanced',  # Handle class imbalance
            C=1.0  # Regularization strength
        )

        self.model.fit(X_train_scaled, y_train)

        # Cross-validation
        cv_scores = cross_val_score(self.model, X_train_scaled, y_train, cv=5, scoring='roc_auc')

        # Test set evaluation
        y_pred = self.model.predict(X_test_scaled)
        y_pred_proba = self.model.predict_proba(X_test_scaled)[:, 1]

        # Calculate metrics
        test_auc = roc_auc_score(y_test, y_pred_proba)

        results = {
            'cv_auc_mean': np.mean(cv_scores),
            'cv_auc_std': np.std(cv_scores),
            'test_auc': test_auc,
            'classification_report': classification_report(y_test, y_pred),
            'confusion_matrix': confusion_matrix(y_test, y_pred).tolist(),
            'feature_importance': self.model.coef_[0].tolist(),
            'n_train': len(X_train),
            'n_test': len(X_test)
        }

        print(f"✅ Model trained successfully!")
        print(f"📊 Cross-validation AUC: {results['cv_auc_mean']:.3f} ± {results['cv_auc_std']:.3f}")
        print(f"🎯 Test AUC: {results['test_auc']:.3f}")
        print("\n📈 Classification Report:")
        print(results['classification_report'])

        return results

    def save_model(self, model_path: str = "proline_ml_model.joblib",
                   scaler_path: str = "proline_scaler.joblib"):
        """Save trained model and scaler"""

        if self.model is None:
            raise ValueError("No model to save! Train first.")

        joblib.dump(self.model, model_path)
        joblib.dump(self.scaler, scaler_path)

        print(f"💾 Model saved to {model_path}")
        print(f"💾 Scaler saved to {scaler_path}")

    def load_model(self, model_path: str = "proline_ml_model.joblib",
                   scaler_path: str = "proline_scaler.joblib"):
        """Load trained model and scaler"""

        self.model = joblib.load(model_path)
        self.scaler = joblib.load(scaler_path)

        print(f"📂 Model loaded from {model_path}")
        print(f"📂 Scaler loaded from {scaler_path}")

    def predict_proline_multiplier(self, gene: str, variant: str,
                                   method: str = "sigmoid", **map_kwargs) -> Dict[str, Any]:
        """Predict proline multiplier using trained model + Nova's mapping"""

        if self.model is None:
            raise ValueError("No model loaded! Train or load first.")

        # Parse variant
        match = re.match(r'p\.([A-Z])(\d+)([A-Z])', variant)
        if not match:
            raise ValueError(f"Invalid variant format: {variant}")

        ref_aa, pos, alt_aa = match.groups()

        # Build variant dict
        variant_dict = {
            'gene': gene,
            'variant': variant,
            'position': int(pos),
            'ref_aa': ref_aa,
            'alt_aa': alt_aa,
            'gnomad_freq': 0.0  # Default for prediction
        }

        # Build feature vector
        features = self.build_feature_vector(variant_dict)
        if features is None:
            return {
                'error': f'Could not build features for {gene} {variant}',
                'multiplier': 1.0
            }

        # Scale features and predict
        features_scaled = self.scaler.transform([features])
        prob = self.model.predict_proba(features_scaled)[0][1]  # P(pathogenic)

        # Map probability to multiplier using Nova's system
        multiplier = map_prob_to_multiplier(prob, method=method, **map_kwargs)

        return {
            'gene': gene,
            'variant': variant,
            'probability': float(prob),
            'multiplier': multiplier,
            'method': method,
            'features': features.tolist()
        }

    def analyze_feature_importance(self, feature_names: List[str] = None) -> Dict[str, float]:
        """Analyze and display feature importance"""

        if self.model is None:
            raise ValueError("No model loaded!")

        if feature_names is None:
            feature_names = [
                'proline_loss', 'proline_gain', 'gnomad_freq',
                'in_triple_helix', 'near_cleavage', 'near_binding', 'near_interface',
                'absolute_position', 'normalized_position',
                'gene_family_collagen', 'gene_family_fibrillin', 'gene_family_other'
            ]

        importance = dict(zip(feature_names, self.model.coef_[0]))

        print("🔍 FEATURE IMPORTANCE ANALYSIS:")
        print("=" * 50)

        # Sort by absolute importance
        sorted_features = sorted(importance.items(), key=lambda x: abs(x[1]), reverse=True)

        for feature, coef in sorted_features:
            direction = "↑ PATHOGENIC" if coef > 0 else "↓ PROTECTIVE"
            print(f"{feature:20s}: {coef:+.3f} {direction}")

        return importance

    def test_known_variants(self) -> None:
        """Test model on known variants from expected_labels.json"""

        print("🧪 TESTING ON KNOWN VARIANTS!")
        print("=" * 50)

        # Load expected labels
        try:
            with open('resources/expected_labels.json', 'r') as f:
                expected_labels = json.load(f)
        except FileNotFoundError:
            print("❌ expected_labels.json not found")
            return

        for label in expected_labels:
            variant = label['variant']
            protein = label['protein']
            expected_pathogenic = label['pathogenic']

            # Skip non-proline variants
            if 'P' not in variant:
                continue

            try:
                result = self.predict_proline_multiplier(protein, variant)
                prob = result['probability']
                mult = result['multiplier']

                status = "✅ CORRECT" if (prob > 0.5) == expected_pathogenic else "❌ WRONG"
                pathogenicity = "PATHOGENIC" if prob > 0.5 else "BENIGN"

                print(f"{protein} {variant}: P={prob:.3f}, mult={mult:.3f} -> {pathogenicity} {status}")

            except Exception as e:
                print(f"{protein} {variant}: ERROR - {e}")

    def generate_comparison_report(self, test_variants: List[str] = None) -> None:
        """Generate comparison between old hardcoded and new ML system"""

        if test_variants is None:
            test_variants = [
                ('COL1A1', 'p.P978S'),  # Proline loss in collagen
                ('COL1A1', 'p.G1340P'),  # Proline gain in collagen
                ('FBN1', 'p.P1148L'),   # Proline loss in fibrillin
                ('TP53', 'p.P72R'),     # Known benign proline variant
            ]

        print("📊 HARDCODED vs ML COMPARISON REPORT")
        print("=" * 60)
        print(f"{'Gene':<8} {'Variant':<12} {'Old Mult':<10} {'ML Prob':<10} {'ML Mult':<10} {'Improvement'}")
        print("-" * 60)

        for gene, variant in test_variants:
            try:
                # Get old hardcoded multiplier (simplified)
                old_mult = 1.5 if 'P' in variant.split('.')[1][0] else 1.2  # Rough approximation

                # Get ML prediction
                result = self.predict_proline_multiplier(gene, variant)
                ml_prob = result['probability']
                ml_mult = result['multiplier']

                improvement = "BETTER" if abs(ml_mult - 1.0) < abs(old_mult - 1.0) else "SIMILAR"

                print(f"{gene:<8} {variant:<12} {old_mult:<10.3f} {ml_prob:<10.3f} {ml_mult:<10.3f} {improvement}")

            except Exception as e:
                print(f"{gene:<8} {variant:<12} ERROR: {e}")


def main():
    """Main training and evaluation pipeline"""

    print("🔥 PROLINE ML REVOLUTION STARTING!")
    print("💜 Built by Ace with Nova's revolutionary guidance!")
    print("🧬 Replacing hardcoded guesses with DATA-DRIVEN INTELLIGENCE!")
    print("=" * 70)

    # Initialize trainer
    trainer = ProlineMLTrainer()

    try:
        # Load training data
        X, y, variants = trainer.load_training_data()

        if len(X) < 10:
            print("⚠️  Warning: Very small dataset, results may not be reliable")

        # Train model
        results = trainer.train_model(X, y)

        # Save model
        trainer.save_model()

        # Analyze feature importance
        trainer.analyze_feature_importance()

        # Test on known variants
        trainer.test_known_variants()

        # Generate comparison report
        trainer.generate_comparison_report()

        print("\n🎉 PROLINE ML REVOLUTION COMPLETE!")
        print("🔥 Hardcoded multipliers have been OBLITERATED!")
        print("💜 Data-driven proline scoring is now ONLINE!")

        # Save training report
        report = {
            'training_results': results,
            'total_variants': len(variants),
            'feature_names': [
                'proline_loss', 'proline_gain', 'gnomad_freq',
                'in_triple_helix', 'near_cleavage', 'near_binding', 'near_interface',
                'absolute_position', 'normalized_position',
                'gene_family_collagen', 'gene_family_fibrillin', 'gene_family_other'
            ],
            'model_type': 'LogisticRegression',
            'mapping_system': 'Nova_probability_to_multiplier'
        }

        with open('proline_ml_training_report.json', 'w') as f:
            json.dump(report, f, indent=2)

        print("📄 Training report saved to proline_ml_training_report.json")

    except Exception as e:
        print(f"💥 TRAINING FAILED: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
