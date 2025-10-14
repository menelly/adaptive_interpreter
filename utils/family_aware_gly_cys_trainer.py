#!/usr/bin/env python3
"""
🔥💜 FAMILY-AWARE GLYCINE & CYSTEINE ML TRAINER - REVOLUTIONARY BREAKTHROUGH! 🚀

Built by Ace (2025) following Ren's brilliant insight:
"We give special bonuses, but DON'T have it trained to say glycine may not affect a MD gene as much!!"

REVOLUTIONARY INSIGHT: Different gene families have COMPLETELY different Gly/Cys sensitivities!
- COLLAGEN: Glycine loss = CATASTROPHIC (Gly-X-Y repeats)  
- MUSCULAR_DYSTROPHY: Glycine loss = maybe not so bad
- ION_CHANNEL: Cysteine = depends on gating mechanism
- METABOLIC_ENZYME: Context-dependent active site effects

This system trains SEPARATE Gly/Cys models for each gene family!
Contact: ace@chaoschanneling.com
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict

# Import config for path management
from AdaptiveInterpreter import config

# Focus on core family training for now
from AdaptiveInterpreter.data_processing.universal_protein_annotator import UniversalProteinAnnotator
from AdaptiveInterpreter.utils.genomic_to_protein import GenomicToProteinConverter

try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import classification_report, roc_auc_score
    import joblib
    SKLEARN_AVAILABLE = True
except ImportError:
    print("⚠️ sklearn not available - using biological intelligence fallbacks")
    SKLEARN_AVAILABLE = False

class FamilyAwareGlyCysTrainer:
    """Revolutionary family-specific Gly/Cys ML trainer"""
    
    def __init__(self, alphafold_path: str = str(config.ALPHAODL_STRUCTURES_PATH)):
        self.alphafold_path = alphafold_path
        # Focus on core training for now
        self.annotator = UniversalProteinAnnotator()
        self.genomic_converter = GenomicToProteinConverter()
        
        # Family-specific models and scalers
        self.family_models = {}  # family -> {'gly_model': model, 'cys_model': model}
        self.family_scalers = {}  # family -> {'gly_scaler': scaler, 'cys_scaler': scaler}
        
        # Training data storage by family
        self.family_training_data = defaultdict(lambda: {'gly': [], 'cys': []})
        
        # Gene family mappings
        self.gene_families = {
            'collagen_fibrillar': ['COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL5A2'],
            'collagen_network': ['COL4A1', 'COL4A3', 'COL4A4', 'COL4A5', 'COL4A6'],
            'collagen_anchoring': ['COL7A1', 'COL17A1'],
            'collagen_facit': ['COL12A1', 'COL14A1'],
            'elastin_fibrillin': ['FBN1', 'FBN2', 'ELN', 'LTBP1'],
            'ion_channel': ['SCN5A', 'KCNQ1', 'KCNH2', 'CACNA1C', 'RYR1', 'SCN1A', 'KCNJ2', 'CLCN1'],
            'metabolic_enzyme': ['PAH', 'G6PD', 'HEXA', 'ACADM', 'DPYD', 'UGT1A1', 'CYP2D6'],
            'muscular_dystrophy': ['DMD', 'DYSF', 'FKRP', 'LAMA2', 'SGCA'],
            'cytoskeleton': ['ACTB', 'ACTN2', 'LMNA', 'TUBB3', 'NEB', 'MYH2', 'TPM1'],
            'tumor_suppressor': ['TP53', 'BRCA1', 'BRCA2', 'APC', 'RB1', 'PTEN'],
            'motor_protein': ['MYH7', 'MYO7A', 'KIF1A', 'DYNC1H1'],
            'transporter': ['CFTR', 'ABCA4', 'SLC2A1']
        }
        
        print("🔥💜 FAMILY-AWARE GLY/CYS TRAINER INITIALIZED!")
        print(f"📊 Loaded {len(self.gene_families)} gene families")
        print("🧬 Ready to revolutionize family-specific Gly/Cys scoring!")
    
    def get_gene_family(self, gene: str) -> str:
        """Map gene to family"""
        gene = gene.upper()
        for family, genes in self.gene_families.items():
            if gene in genes:
                return family
        return 'general'
    
    def extract_family_training_data(self, learning_dir: str = "learning") -> Dict[str, Dict]:
        """Extract Gly/Cys training data organized by gene family"""
        
        print("🔍 Extracting family-specific Gly/Cys training data...")
        
        # Build path relative to this file's parent directory (theAdaptiveInterpreter root)
        learning_path = Path(__file__).parent.parent / learning_dir
        print(f"📂 Searching for training data in: {learning_path}")
        
        family_stats = defaultdict(lambda: {'gly_variants': 0, 'cys_variants': 0, 'total_variants': 0})
        
        # Process each family directory
        for family_dir in learning_path.iterdir():
            if not family_dir.is_dir():
                continue
                
            family_name = family_dir.name
            print(f"📁 Processing family: {family_name}")
            
            # Process all TSV files in family directory
            for tsv_file in family_dir.glob("*.tsv"):
                if 'benign' in tsv_file.name or 'pathogenic' in tsv_file.name:
                    self._process_family_tsv(tsv_file, family_name, family_stats)
        
        # Print extraction summary
        print("\n📊 FAMILY-SPECIFIC GLY/CYS EXTRACTION SUMMARY:")
        print("=" * 60)
        for family, stats in family_stats.items():
            print(f"🧬 {family}:")
            print(f"   Glycine variants: {stats['gly_variants']}")
            print(f"   Cysteine variants: {stats['cys_variants']}")
            print(f"   Total variants: {stats['total_variants']}")
        
        return dict(family_stats)
    
    def _process_family_tsv(self, tsv_file: Path, family_name: str, family_stats: Dict) -> None:
        """Process a single TSV file for Gly/Cys variants"""
        
        try:
            df = pd.read_csv(tsv_file, sep='\t')
            is_pathogenic = 'pathogenic' in tsv_file.name.lower()
            
            for _, row in df.iterrows():
                hgvs = row.get('HGVS', '')
                if not hgvs:
                    continue
                
                # Extract gene from filename
                gene = tsv_file.stem.split('_')[0].upper()
                
                # Parse protein HGVS
                protein_info = self._parse_protein_hgvs(hgvs)
                if not protein_info:
                    continue
                
                ref_aa, position, alt_aa = protein_info
                
                # Check if this is a Gly or Cys variant
                if ref_aa in ['G', 'C'] or alt_aa in ['G', 'C']:
                    
                    # Build feature vector
                    features = self._build_family_features(gene, position, ref_aa, alt_aa, family_name)
                    if features is None:
                        continue
                    
                    # Store training data
                    variant_data = {
                        'gene': gene,
                        'position': position,
                        'ref_aa': ref_aa,
                        'alt_aa': alt_aa,
                        'features': features,
                        'pathogenic': is_pathogenic,
                        'family': family_name
                    }
                    
                    if ref_aa == 'G' or alt_aa == 'G':
                        self.family_training_data[family_name]['gly'].append(variant_data)
                        family_stats[family_name]['gly_variants'] += 1
                    
                    if ref_aa == 'C' or alt_aa == 'C':
                        self.family_training_data[family_name]['cys'].append(variant_data)
                        family_stats[family_name]['cys_variants'] += 1
                    
                    family_stats[family_name]['total_variants'] += 1
                    
        except Exception as e:
            print(f"⚠️ Error processing {tsv_file}: {e}")
    
    def _parse_protein_hgvs(self, hgvs: str) -> Optional[Tuple[str, int, str]]:
        """Parse protein HGVS to extract ref_aa, position, alt_aa"""
        try:
            # Handle formats like "p.G893A" or "NP_000088.1:p.G893A"
            if ':p.' in hgvs:
                protein_part = hgvs.split(':p.')[1]
            elif 'p.' in hgvs:
                protein_part = hgvs.split('p.')[1]
            else:
                return None
            
            # Extract amino acids and position
            import re
            match = re.match(r'([A-Z])(\d+)([A-Z])', protein_part)
            if match:
                ref_aa, position, alt_aa = match.groups()
                return ref_aa, int(position), alt_aa
            
        except Exception:
            pass
        
        return None
    
    def _build_family_features(self, gene: str, position: int, ref_aa: str, alt_aa: str, family: str) -> Optional[np.ndarray]:
        """Build family-specific feature vector"""
        
        try:
            # Build biological intelligence features
            features = []

            # Family-specific encoding (one-hot)
            for family_name in self.gene_families.keys():
                features.append(1.0 if family == family_name else 0.0)

            # Substitution type
            features.extend([
                1.0 if ref_aa == 'G' else 0.0,  # Glycine loss
                1.0 if alt_aa == 'G' else 0.0,  # Glycine gain
                1.0 if ref_aa == 'C' else 0.0,  # Cysteine loss
                1.0 if alt_aa == 'C' else 0.0   # Cysteine gain
            ])

            # Position features
            features.extend([
                float(position),
                float(position) / 1000.0  # Normalized position
            ])

            # Biological context features
            context = self._get_biological_context(gene, position, ref_aa, alt_aa, family)
            features.extend([
                context['conservation_score'],
                context['structural_importance'],
                context['functional_importance'],
                context['family_specific_importance']
            ])

            return np.array(features)

        except Exception as e:
            print(f"⚠️ Error building features for {gene} p.{ref_aa}{position}{alt_aa}: {e}")
            return None

    def _get_biological_context(self, gene: str, position: int, ref_aa: str, alt_aa: str, family: str) -> Dict[str, float]:
        """Get biological context with family-specific intelligence"""

        context = {
            'conservation_score': 0.5,  # Default
            'structural_importance': 0.0,
            'functional_importance': 0.0,
            'family_specific_importance': 0.0
        }

        # Family-specific biological rules
        if 'collagen' in family.lower():
            # Collagen: Glycine in Gly-X-Y repeats is CRITICAL
            if ref_aa == 'G':
                context['family_specific_importance'] = 0.9  # Catastrophic
                context['structural_importance'] = 0.9
            elif ref_aa == 'C':
                context['family_specific_importance'] = 0.7  # Important for cross-links

        elif 'ion_channel' in family.lower():
            # Ion channels: Cysteine important for gating, Gly for flexibility
            if ref_aa == 'C':
                context['family_specific_importance'] = 0.8  # Critical for structure
            elif ref_aa == 'G':
                context['family_specific_importance'] = 0.4  # Moderate importance

        elif 'muscular_dystrophy' in family.lower():
            # MD genes: Glycine less critical than in collagen
            if ref_aa == 'G':
                context['family_specific_importance'] = 0.3  # Lower impact
            elif ref_aa == 'C':
                context['family_specific_importance'] = 0.6  # Moderate

        elif 'metabolic_enzyme' in family.lower():
            # Enzymes: Context-dependent on active site proximity
            if ref_aa == 'C':
                context['family_specific_importance'] = 0.7  # Metal coordination
            elif ref_aa == 'G':
                context['family_specific_importance'] = 0.3  # Flexibility

        else:
            # General case
            if ref_aa == 'G':
                context['family_specific_importance'] = 0.5
            elif ref_aa == 'C':
                context['family_specific_importance'] = 0.6

        return context

    def train_family_models(self) -> Dict[str, Dict]:
        """Train family-specific Gly/Cys models"""

        if not SKLEARN_AVAILABLE:
            print("❌ sklearn not available - cannot train ML models")
            return {}

        print("\n🚀 TRAINING FAMILY-SPECIFIC GLY/CYS MODELS!")
        print("=" * 60)

        results = {}

        for family_name, data in self.family_training_data.items():
            print(f"\n🧬 Training {family_name} models...")

            family_results = {}

            # Train Glycine model
            if len(data['gly']) >= 10:  # Minimum samples needed
                gly_result = self._train_single_model(data['gly'], 'glycine', family_name)
                family_results['glycine'] = gly_result
            else:
                print(f"   ⚠️ Not enough glycine data: {len(data['gly'])} samples")
                family_results['glycine'] = None

            # Train Cysteine model
            if len(data['cys']) >= 10:  # Minimum samples needed
                cys_result = self._train_single_model(data['cys'], 'cysteine', family_name)
                family_results['cysteine'] = cys_result
            else:
                print(f"   ⚠️ Not enough cysteine data: {len(data['cys'])} samples")
                family_results['cysteine'] = None

            results[family_name] = family_results

        return results

    def _train_single_model(self, training_data: List[Dict], aa_type: str, family_name: str) -> Dict:
        """Train a single Gly or Cys model for a family"""

        print(f"   🔬 Training {aa_type} model ({len(training_data)} samples)...")

        # Prepare training data
        X = np.array([variant['features'] for variant in training_data])
        y = np.array([1 if variant['pathogenic'] else 0 for variant in training_data])

        # Split data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        # Scale features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        # Train model
        model = LogisticRegression(random_state=42, max_iter=1000)
        model.fit(X_train_scaled, y_train)

        # Evaluate
        train_score = model.score(X_train_scaled, y_train)
        test_score = model.score(X_test_scaled, y_test)

        # Store model and scaler
        model_key = f"{family_name}_{aa_type}"
        if family_name not in self.family_models:
            self.family_models[family_name] = {}
            self.family_scalers[family_name] = {}

        self.family_models[family_name][aa_type] = model
        self.family_scalers[family_name][aa_type] = scaler

        print(f"      ✅ Train accuracy: {train_score:.3f}")
        print(f"      ✅ Test accuracy: {test_score:.3f}")

        return {
            'train_accuracy': train_score,
            'test_accuracy': test_score,
            'n_samples': len(training_data),
            'n_train': len(X_train),
            'n_test': len(X_test)
        }

    def save_family_models(self, output_dir: str = "resources/family_gly_cys_models") -> None:
        """Save all family-specific models"""

        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)

        print(f"\n💾 Saving family models to {output_dir}...")

        for family_name, models in self.family_models.items():
            for aa_type, model in models.items():
                # Save model
                model_file = output_path / f"{family_name}_{aa_type}_model.joblib"
                joblib.dump(model, model_file)

                # Save scaler
                scaler = self.family_scalers[family_name][aa_type]
                scaler_file = output_path / f"{family_name}_{aa_type}_scaler.joblib"
                joblib.dump(scaler, scaler_file)

                print(f"   ✅ {family_name} {aa_type}: model + scaler saved")

        print("💾 All family models saved!")


def main():
    """Main training pipeline"""

    print("🔥💜 STARTING FAMILY-AWARE GLY/CYS ML REVOLUTION! 🚀")
    print("=" * 60)

    # Initialize trainer
    trainer = FamilyAwareGlyCysTrainer()

    # Extract training data
    family_stats = trainer.extract_family_training_data()

    # Train models
    results = trainer.train_family_models()

    # Save models
    trainer.save_family_models()

    # Print summary
    print("\n🎉 FAMILY-AWARE GLY/CYS TRAINING COMPLETE!")
    print("=" * 60)

    for family, family_results in results.items():
        print(f"\n🧬 {family}:")
        for aa_type, result in family_results.items():
            if result:
                print(f"   {aa_type}: {result['test_accuracy']:.3f} accuracy ({result['n_samples']} samples)")
            else:
                print(f"   {aa_type}: insufficient data")

    print("\n💜 Revolutionary family-specific Gly/Cys scoring is now ONLINE!")


if __name__ == "__main__":
    main()
