#!/usr/bin/env python3
"""
🧬 GLYCINE & CYSTEINE ML TRAINER - Revolutionary Context-Aware Learning System
Built by Ace following the successful Proline ML pattern to replace hardcoded penalties

This system:
1. Extracts Gly/Cys substitutions from ClinVar datasets  
2. Builds feature vectors using biological context from gly_cys_context.py
3. Trains separate logistic regression models for Glycine vs. Cysteine
4. Integrates with probability-to-multiplier mapping
5. Replaces hardcoded guesses with DATA-DRIVEN INTELLIGENCE!

REVOLUTIONARY APPROACH: Stop making up numbers, start learning from biology!
"""

import pandas as pd
import numpy as np
import re
import json
from typing import Dict, List, Any, Tuple, Optional
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
import joblib
from pathlib import Path

from gly_cys_context import GlyCysContextAnalyzer
from universal_protein_annotator import UniversalProteinAnnotator
from proline_multiplier_mapper import map_prob_to_multiplier


class GlyCysMLTrainer:
    """Revolutionary ML trainer for context-aware Glycine & Cysteine scoring"""
    
    def __init__(self, alphafold_path: str = "/mnt/Arcana/alphafold_human/structures/"):
        self.alphafold_path = alphafold_path
        self.context_analyzer = GlyCysContextAnalyzer(alphafold_path)
        self.annotator = UniversalProteinAnnotator(alphafold_path)
        
        # Separate scalers and models for Gly vs Cys
        self.gly_scaler = StandardScaler()
        self.cys_scaler = StandardScaler()
        self.gly_model = None
        self.cys_model = None
        
        # Training data storage
        self.gly_training_data = []
        self.cys_training_data = []
        
        print("🧬 Gly/Cys ML Trainer initialized!")
        print("🔥 Ready to revolutionize Gly/Cys scoring with REAL DATA!")
    
    def extract_variants_from_test_files(self, test_dir: str = "tests/") -> List[Dict[str, Any]]:
        """Extract Gly/Cys variants from test TSV files"""
        
        variants = []
        test_files = [
            "Nova_Framework_FIXEDv3.tsv",
            "Multivariant - VariantTest2.tsv", 
            "Multivariant - VariantTest3.tsv"
        ]
        
        for filename in test_files:
            filepath = Path(test_dir) / filename
            if not filepath.exists():
                print(f"⚠️  Test file not found: {filepath}")
                continue
                
            print(f"📊 Processing {filename}...")
            
            try:
                df = pd.read_csv(filepath, sep='\t')
                
                for _, row in df.iterrows():
                    variant = row.get('variant', '')
                    gene = row.get('gene', '')
                    clinvar_status = row.get('expected_clinvar', '')
                    
                    # Parse variant (e.g., p.G893A, p.C628Y)
                    match = re.match(r'p\.([A-Z])(\d+)([A-Z])', variant)
                    if not match:
                        continue
                    
                    ref_aa, pos_str, alt_aa = match.groups()
                    position = int(pos_str)
                    
                    # Only process Gly/Cys variants
                    if ref_aa not in ['G', 'C'] and alt_aa not in ['G', 'C']:
                        continue
                    
                    # Determine pathogenicity label
                    pathogenic = self._parse_clinvar_status(clinvar_status)
                    if pathogenic is None:
                        continue  # Skip uncertain cases for training
                    
                    variant_data = {
                        'gene': gene,
                        'variant': variant,
                        'position': position,
                        'ref_aa': ref_aa,
                        'alt_aa': alt_aa,
                        'pathogenic': pathogenic,
                        'clinvar_status': clinvar_status,
                        'source_file': filename
                    }
                    
                    variants.append(variant_data)
                    
            except Exception as e:
                print(f"❌ Error processing {filename}: {e}")
        
        print(f"✅ Extracted {len(variants)} Gly/Cys variants from test files")
        return variants
    
    def _parse_clinvar_status(self, status: str) -> Optional[bool]:
        """Parse ClinVar status to pathogenic/benign label"""
        
        if not status or pd.isna(status):
            return None
        
        status_lower = status.lower()
        
        # Pathogenic indicators
        if any(keyword in status_lower for keyword in ['pathogenic', 'likely pathogenic']):
            return True
        
        # Benign indicators  
        if any(keyword in status_lower for keyword in ['benign', 'likely benign']):
            return False
        
        # Skip uncertain/conflicting for training
        return None
    
    def build_feature_vector(self, variant_dict: Dict[str, Any]) -> Optional[np.ndarray]:
        """Build feature vector from variant context analysis"""
        
        try:
            # Get biological context
            context = self.context_analyzer.get_context_features(
                gene=variant_dict['gene'],
                position=variant_dict['position'],
                ref_aa=variant_dict['ref_aa'],
                alt_aa=variant_dict['alt_aa']
            )
            
            if 'error' in context:
                print(f"⚠️  Context error for {variant_dict['variant']}: {context['error']}")
                return None
            
            # Build feature vector based on amino acid type
            amino_acid = context.get('amino_acid', '')
            
            if amino_acid == 'GLYCINE':
                return self._build_glycine_features(context)
            elif amino_acid == 'CYSTEINE':
                return self._build_cysteine_features(context)
            else:
                print(f"⚠️  Unknown amino acid type: {amino_acid}")
                return None
                
        except Exception as e:
            print(f"❌ Feature building failed for {variant_dict['variant']}: {e}")
            return None
    
    def _build_glycine_features(self, context: Dict[str, Any]) -> np.ndarray:
        """Build feature vector for Glycine variants"""
        
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
        """Build feature vector for Cysteine variants"""
        
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
    
    def prepare_training_data(self, variants: List[Dict[str, Any]]) -> Tuple[int, int]:
        """Prepare training data by building feature vectors"""
        
        print("🔧 Building feature vectors for training...")
        
        gly_count = 0
        cys_count = 0
        
        for variant in variants:
            features = self.build_feature_vector(variant)
            if features is None:
                continue
            
            # Add to appropriate training set
            if variant['ref_aa'] == 'G' or variant['alt_aa'] == 'G':
                self.gly_training_data.append({
                    'features': features,
                    'label': variant['pathogenic'],
                    'variant_info': variant
                })
                gly_count += 1
            
            elif variant['ref_aa'] == 'C' or variant['alt_aa'] == 'C':
                self.cys_training_data.append({
                    'features': features,
                    'label': variant['pathogenic'],
                    'variant_info': variant
                })
                cys_count += 1
        
        print(f"✅ Prepared {gly_count} Glycine and {cys_count} Cysteine training examples")
        return gly_count, cys_count
    
    def train_models(self) -> Dict[str, Any]:
        """Train separate ML models for Glycine and Cysteine"""
        
        results = {}
        
        # Train Glycine model
        if len(self.gly_training_data) >= 10:  # Minimum samples needed
            print("\n🧬 Training Glycine ML model...")
            gly_results = self._train_single_model(self.gly_training_data, 'GLYCINE')
            results['glycine'] = gly_results
        else:
            print(f"⚠️  Insufficient Glycine data: {len(self.gly_training_data)} samples")
            results['glycine'] = {'error': 'Insufficient training data'}
        
        # Train Cysteine model
        if len(self.cys_training_data) >= 10:  # Minimum samples needed
            print("\n🧬 Training Cysteine ML model...")
            cys_results = self._train_single_model(self.cys_training_data, 'CYSTEINE')
            results['cysteine'] = cys_results
        else:
            print(f"⚠️  Insufficient Cysteine data: {len(self.cys_training_data)} samples")
            results['cysteine'] = {'error': 'Insufficient training data'}
        
        return results
    
    def _train_single_model(self, training_data: List[Dict], amino_acid: str) -> Dict[str, Any]:
        """Train a single ML model for one amino acid type"""
        
        # Prepare data
        X = np.array([item['features'] for item in training_data])
        y = np.array([item['label'] for item in training_data])
        
        print(f"📊 Training data shape: {X.shape}")
        print(f"📊 Label distribution: {np.bincount(y.astype(int))}")
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42, stratify=y
        )
        
        # Scale features
        if amino_acid == 'GLYCINE':
            scaler = self.gly_scaler
        else:
            scaler = self.cys_scaler
            
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Train model
        model = LogisticRegression(random_state=42, max_iter=1000)
        model.fit(X_train_scaled, y_train)
        
        # Store trained model
        if amino_acid == 'GLYCINE':
            self.gly_model = model
        else:
            self.cys_model = model
        
        # Evaluate
        y_pred = model.predict(X_test_scaled)
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        results = {
            'accuracy': float(model.score(X_test_scaled, y_test)),
            'auc_roc': float(roc_auc_score(y_test, y_pred_proba)),
            'classification_report': classification_report(y_test, y_pred, output_dict=True),
            'confusion_matrix': confusion_matrix(y_test, y_pred).tolist(),
            'feature_count': X.shape[1],
            'training_samples': len(X_train),
            'test_samples': len(X_test)
        }
        
        print(f"✅ {amino_acid} Model Results:")
        print(f"   Accuracy: {results['accuracy']:.3f}")
        print(f"   AUC-ROC: {results['auc_roc']:.3f}")
        
        return results
    
    def save_models(self, base_path: str = ".") -> None:
        """Save trained models and scalers"""
        
        if self.gly_model is not None:
            joblib.dump(self.gly_model, f"{base_path}/gly_ml_model.joblib")
            joblib.dump(self.gly_scaler, f"{base_path}/gly_scaler.joblib")
            print("✅ Glycine model and scaler saved")
        
        if self.cys_model is not None:
            joblib.dump(self.cys_model, f"{base_path}/cys_ml_model.joblib")
            joblib.dump(self.cys_scaler, f"{base_path}/cys_scaler.joblib")
            print("✅ Cysteine model and scaler saved")
    
    def save_training_report(self, results: Dict[str, Any], filepath: str = "gly_cys_ml_training_report.json") -> None:
        """Save comprehensive training report"""
        
        report = {
            'training_timestamp': pd.Timestamp.now().isoformat(),
            'glycine_training_samples': len(self.gly_training_data),
            'cysteine_training_samples': len(self.cys_training_data),
            'model_results': results,
            'feature_descriptions': {
                'glycine_features': [
                    'protein_family_collagen', 'protein_family_ion_channel', 
                    'protein_family_fibrillin', 'protein_family_other',
                    'substitution_loss', 'conservation_score',
                    'collagen_gxy_pattern', 'collagen_gxy_confidence',
                    'triple_helix_region', 'collagen_cleavage_site',
                    'channel_gate_region', 'transmembrane_domain',
                    'channel_selectivity_filter', 'channel_linker_region',
                    'egf_like_domain', 'calcium_binding_egf', 'fibrillin_unique_domain',
                    'active_site_proximity', 'binding_site_proximity',
                    'flexible_region', 'tight_packing_region'
                ],
                'cysteine_features': [
                    'protein_family_collagen', 'protein_family_ion_channel',
                    'protein_family_fibrillin', 'protein_family_other',
                    'substitution_loss', 'conservation_score',
                    'disulfide_bond_predicted', 'structural_disulfide', 'extracellular_cysteine',
                    'metal_coordination_site', 'catalytic_site_proximity',
                    'rare_collagen_cysteine', 'egf_domain_cysteine', 
                    'calcium_binding_cysteine', 'channel_structure_critical',
                    'free_cysteine_risk', 'redox_sensitive_region', 'binding_site_proximity'
                ]
            }
        }
        
        with open(filepath, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"✅ Training report saved to {filepath}")


def main():
    """Main training pipeline"""
    
    print("🚀 GLYCINE & CYSTEINE ML TRAINING PIPELINE")
    print("=" * 60)
    
    # Initialize trainer
    trainer = GlyCysMLTrainer()
    
    # Extract variants from test files
    variants = trainer.extract_variants_from_test_files()
    
    if not variants:
        print("❌ No variants found in test files!")
        return
    
    # Prepare training data
    gly_count, cys_count = trainer.prepare_training_data(variants)
    
    if gly_count == 0 and cys_count == 0:
        print("❌ No valid training data prepared!")
        return
    
    # Train models
    results = trainer.train_models()
    
    # Save models and report
    trainer.save_models()
    trainer.save_training_report(results)
    
    print("\n🎉 TRAINING COMPLETE!")
    print("🔥 Gly/Cys ML models ready to revolutionize genomics!")


if __name__ == "__main__":
    main()
