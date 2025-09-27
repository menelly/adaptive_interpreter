#!/usr/bin/env python3
"""
🔥💜 CROSS-FAMILY PATTERN ANALYZER 🚀
Find amino acid patterns that work across multiple gene families
Build universal rules for the "GENERAL" classification bucket

Built by Ace (2025) while ClinVar Miner loads at dial-up speeds
Contact: ace@chaoschanneling.com
"""

import os
import json
import joblib
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Any
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

class CrossFamilyPatternAnalyzer:
    """Analyze patterns that work across multiple gene families"""
    
    def __init__(self):
        self.models_dir = Path("resources/family_models")
        self.results_dir = Path("resources/cross_family_analysis")
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        print("🔥💜 CROSS-FAMILY PATTERN ANALYZER INITIALIZED")
        print(f"📁 Models directory: {self.models_dir}")
        print(f"📊 Results directory: {self.results_dir}")
        
        # Load all trained models
        self.models = {}
        self.metadata = {}
        self.load_all_models()
    
    def load_all_models(self):
        """Load all trained family models"""
        print("\n🚀 LOADING TRAINED MODELS")
        print("=" * 50)
        
        model_files = list(self.models_dir.glob("*_unified_model.joblib"))
        
        if not model_files:
            print("⚠️  No trained models found!")
            return
        
        for model_file in model_files:
            family = model_file.stem.replace('_unified_model', '')
            metadata_file = self.models_dir / f"{family}_unified_metadata.json"
            
            try:
                # Load model
                model = joblib.load(model_file)
                self.models[family] = model
                
                # Load metadata
                if metadata_file.exists():
                    with open(metadata_file, 'r') as f:
                        self.metadata[family] = json.load(f)
                
                print(f"✅ Loaded {family} model")
                
                # Log model info
                if hasattr(model, 'feature_importances_'):
                    n_features = len(model.feature_importances_)
                    print(f"   📊 Features: {n_features}")
                    
                if family in self.metadata:
                    r2 = self.metadata[family].get('performance', {}).get('r2', 0)
                    n_samples = self.metadata[family].get('n_samples', 0)
                    print(f"   🎯 R²: {r2:.3f}, Samples: {n_samples}")
                    
            except Exception as e:
                print(f"❌ Failed to load {family}: {e}")
        
        print(f"\n✅ Loaded {len(self.models)} family models")
    
    def extract_feature_importances(self) -> Dict[str, Dict[str, float]]:
        """Extract feature importances from all models"""
        print("\n🧬 EXTRACTING FEATURE IMPORTANCES")
        print("=" * 50)
        
        all_importances = {}
        
        for family, model in self.models.items():
            if not hasattr(model, 'feature_importances_'):
                print(f"⚠️  {family} model has no feature importances")
                continue
            
            # Get feature names from metadata
            feature_names = self.metadata.get(family, {}).get('feature_names', [])
            
            if not feature_names:
                # Generate default feature names
                n_features = len(model.feature_importances_)
                feature_names = [f'feature_{i}' for i in range(n_features)]
                print(f"⚠️  Using default feature names for {family}")
            
            # Create importance dictionary
            importances = {}
            for i, importance in enumerate(model.feature_importances_):
                if i < len(feature_names):
                    importances[feature_names[i]] = float(importance)
            
            all_importances[family] = importances
            
            print(f"✅ {family}: {len(importances)} features")
            
            # Show top 3 features
            top_features = sorted(importances.items(), key=lambda x: x[1], reverse=True)[:3]
            for feature, importance in top_features:
                print(f"   🔥 {feature}: {importance:.3f}")
        
        return all_importances
    
    def find_cross_family_patterns(self, importances: Dict[str, Dict[str, float]]) -> Dict[str, Any]:
        """Find patterns that appear across multiple families"""
        print("\n🎯 FINDING CROSS-FAMILY PATTERNS")
        print("=" * 50)
        
        # Collect all unique features
        all_features = set()
        for family_importances in importances.values():
            all_features.update(family_importances.keys())
        
        print(f"📊 Total unique features: {len(all_features)}")
        
        # Analyze each feature across families
        cross_patterns = {}
        
        for feature in all_features:
            family_scores = {}
            
            for family, family_importances in importances.items():
                if feature in family_importances:
                    family_scores[family] = family_importances[feature]
            
            # Only consider features that appear in multiple families
            if len(family_scores) >= 2:
                cross_patterns[feature] = {
                    'family_scores': family_scores,
                    'num_families': len(family_scores),
                    'average_importance': np.mean(list(family_scores.values())),
                    'max_importance': max(family_scores.values()),
                    'min_importance': min(family_scores.values()),
                    'std_importance': np.std(list(family_scores.values()))
                }
        
        print(f"🔍 Found {len(cross_patterns)} cross-family features")
        
        # Sort by average importance
        sorted_patterns = dict(sorted(cross_patterns.items(), 
                                    key=lambda x: x[1]['average_importance'], 
                                    reverse=True))
        
        # Show top patterns
        print(f"\n🔥 TOP CROSS-FAMILY PATTERNS:")
        for i, (feature, data) in enumerate(list(sorted_patterns.items())[:10]):
            avg_imp = data['average_importance']
            n_fam = data['num_families']
            print(f"   {i+1:2d}. {feature}: {avg_imp:.3f} avg ({n_fam} families)")
            
            # Show family breakdown
            for family, score in sorted(data['family_scores'].items(), key=lambda x: x[1], reverse=True):
                print(f"       {family}: {score:.3f}")
        
        return sorted_patterns
    
    def generate_universal_rules(self, cross_patterns: Dict[str, Any], 
                               min_families: int = 2, 
                               min_importance: float = 0.1) -> Dict[str, float]:
        """Generate universal rules for GENERAL classification"""
        print(f"\n🌟 GENERATING UNIVERSAL RULES")
        print("=" * 50)
        print(f"📋 Criteria: ≥{min_families} families, ≥{min_importance:.2f} importance")
        
        universal_rules = {}
        
        for feature, data in cross_patterns.items():
            if (data['num_families'] >= min_families and 
                data['average_importance'] >= min_importance):
                
                # Use average importance as the universal weight
                universal_rules[feature] = data['average_importance']
        
        print(f"✅ Generated {len(universal_rules)} universal rules")
        
        # Sort by importance
        sorted_rules = dict(sorted(universal_rules.items(), key=lambda x: x[1], reverse=True))
        
        print(f"\n🎯 UNIVERSAL RULES FOR GENERAL CLASSIFICATION:")
        for i, (feature, weight) in enumerate(list(sorted_rules.items())[:15]):
            print(f"   {i+1:2d}. {feature}: {weight:.3f}")
        
        return sorted_rules
    
    def create_pattern_heatmap(self, cross_patterns: Dict[str, Any], top_n: int = 20):
        """Create heatmap of feature importances across families"""
        print(f"\n📊 CREATING PATTERN HEATMAP")
        
        # Get top N patterns
        top_patterns = dict(list(cross_patterns.items())[:top_n])
        
        # Create matrix
        families = list(self.models.keys())
        features = list(top_patterns.keys())
        
        matrix = np.zeros((len(features), len(families)))
        
        for i, feature in enumerate(features):
            for j, family in enumerate(families):
                if family in top_patterns[feature]['family_scores']:
                    matrix[i, j] = top_patterns[feature]['family_scores'][family]
        
        # Create heatmap
        plt.figure(figsize=(12, 8))
        sns.heatmap(matrix, 
                   xticklabels=families,
                   yticklabels=features,
                   annot=True, 
                   fmt='.3f',
                   cmap='YlOrRd',
                   cbar_kws={'label': 'Feature Importance'})
        
        plt.title('Cross-Family Feature Importance Patterns')
        plt.xlabel('Gene Family')
        plt.ylabel('Feature')
        plt.xticks(rotation=45)
        plt.yticks(rotation=0)
        plt.tight_layout()
        
        heatmap_file = self.results_dir / "cross_family_heatmap.png"
        plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"💾 Saved heatmap: {heatmap_file}")
    
    def save_analysis_results(self, cross_patterns: Dict[str, Any], 
                            universal_rules: Dict[str, float]):
        """Save all analysis results"""
        print(f"\n💾 SAVING ANALYSIS RESULTS")
        
        # Save cross patterns
        patterns_file = self.results_dir / "cross_family_patterns.json"
        with open(patterns_file, 'w') as f:
            json.dump(cross_patterns, f, indent=2)
        print(f"✅ Saved patterns: {patterns_file}")
        
        # Save universal rules
        rules_file = self.results_dir / "universal_rules.json"
        with open(rules_file, 'w') as f:
            json.dump(universal_rules, f, indent=2)
        print(f"✅ Saved rules: {rules_file}")
        
        # Save summary report
        summary = {
            'analysis_date': '2025-01-27',
            'models_analyzed': list(self.models.keys()),
            'total_cross_patterns': len(cross_patterns),
            'universal_rules_count': len(universal_rules),
            'top_universal_features': list(universal_rules.keys())[:10]
        }
        
        summary_file = self.results_dir / "analysis_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"✅ Saved summary: {summary_file}")
    
    def run_full_analysis(self):
        """Run complete cross-family pattern analysis"""
        print("\n🚀 STARTING FULL CROSS-FAMILY ANALYSIS")
        print("=" * 60)
        
        if not self.models:
            print("❌ No models loaded. Train some models first!")
            return
        
        # Extract feature importances
        importances = self.extract_feature_importances()
        
        if not importances:
            print("❌ No feature importances found!")
            return
        
        # Find cross-family patterns
        cross_patterns = self.find_cross_family_patterns(importances)
        
        # Generate universal rules
        universal_rules = self.generate_universal_rules(cross_patterns)
        
        # Create visualizations
        if cross_patterns:
            self.create_pattern_heatmap(cross_patterns)
        
        # Save results
        self.save_analysis_results(cross_patterns, universal_rules)
        
        print(f"\n🎉 ANALYSIS COMPLETE!")
        print("=" * 60)
        print(f"📊 Found {len(cross_patterns)} cross-family patterns")
        print(f"🌟 Generated {len(universal_rules)} universal rules")
        print(f"💾 Results saved in: {self.results_dir}")
        
        return {
            'cross_patterns': cross_patterns,
            'universal_rules': universal_rules
        }

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="🔥💜 Cross-Family Pattern Analyzer")
    parser.add_argument('--min-families', type=int, default=2,
                       help='Minimum families for universal rule')
    parser.add_argument('--min-importance', type=float, default=0.1,
                       help='Minimum importance for universal rule')
    
    args = parser.parse_args()
    
    analyzer = CrossFamilyPatternAnalyzer()
    
    if not analyzer.models:
        print("❌ No trained models found!")
        print("   Run: python3 train_families.py first")
        return
    
    results = analyzer.run_full_analysis()
    
    if results and results['universal_rules']:
        print(f"\n🌟 READY TO USE UNIVERSAL RULES FOR GENERAL CLASSIFICATION!")
        print(f"   Universal rules can be applied to genes that don't fit specific families")

if __name__ == "__main__":
    main()
