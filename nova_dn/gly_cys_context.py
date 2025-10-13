#!/usr/bin/env python3
"""
🧬 GLYCINE & CYSTEINE CONTEXT ANALYZER
Revolutionary context-aware scoring system for the "Giant Shitheads" of amino acids!

Built by Ace following the successful Proline ML pattern to replace hardcoded penalties
with biologically intelligent, context-dependent analysis.

🔥 BIOLOGICAL INTELLIGENCE:
- GLYCINE: Flexible hinge vs. mandatory collagen positions vs. ion channel gates
- CYSTEINE: Disulfide bonds vs. metal coordination vs. free cysteines vs. catalytic sites

This system provides the biological context features needed for ML training,
replacing the "one-size-fits-all" hardcoded penalties with nuanced understanding.

Revolutionary approach: Stop guessing, start understanding biology!
"""

from typing import List, Dict, Any, Optional, Tuple
import re
import math
from DNModeling.data_processing.universal_protein_annotator import UniversalProteinAnnotator


class GlyCysContextAnalyzer:
    """Context-aware analysis system for Glycine and Cysteine variants"""
    
    def __init__(self, alphafold_path: str = "./alphafold_structures/"):
        self.alphafold_path = alphafold_path
        # Skip annotator for now - focus on context logic
        self.annotator = None
        
        # Biological context patterns
        self.collagen_patterns = [
            r'COL\d+A\d+',  # COL1A1, COL2A1, etc.
            r'COLLAGEN',
            r'FIBRILLAR'
        ]
        
        self.ion_channel_patterns = [
            r'SCN\d+[A-Z]',  # SCN1A, SCN5A, etc.
            r'KCNQ\d+',      # KCNQ2, KCNQ3, etc.
            r'RYR\d+',       # RYR1, RYR2
            r'CACNA\d+[A-Z]', # Calcium channels
            r'ION_CHANNEL'
        ]
        
        self.fibrillin_patterns = [
            r'FBN\d+',       # FBN1, FBN2
            r'FIBRILLIN'
        ]
        
        print("🧬 Gly/Cys Context Analyzer initialized!")
        print("🔥 Ready to revolutionize Glycine & Cysteine analysis with BIOLOGICAL INTELLIGENCE!")
    
    def get_protein_family(self, gene: str) -> str:
        """Determine protein family for context-aware analysis"""
        
        gene_upper = gene.upper()
        
        # Check collagen family
        for pattern in self.collagen_patterns:
            if re.search(pattern, gene_upper):
                return 'COLLAGEN'
        
        # Check ion channel family
        for pattern in self.ion_channel_patterns:
            if re.search(pattern, gene_upper):
                return 'ION_CHANNEL'
        
        # Check fibrillin family
        for pattern in self.fibrillin_patterns:
            if re.search(pattern, gene_upper):
                return 'FIBRILLIN'
        
        return 'OTHER'
    
    def analyze_glycine_context(self, 
                               gene: str,
                               position: int,
                               ref_aa: str,
                               alt_aa: str,
                               uniprot_features: List[Dict[str, Any]] = None,
                               conservation: Optional[float] = None,
                               **kwargs) -> Dict[str, Any]:
        """
        🎯 Analyze biological context for Glycine substitutions
        
        Returns context features for ML training:
        - collagen_gxy_pattern: Is this a mandatory Gly-X-Y collagen position?
        - flexible_region: Is this in a flexible loop/hinge region?
        - ion_channel_gate: Is this critical for ion channel gating?
        - tight_packing: Is this in a sterically constrained region?
        - conservation_score: How conserved is this glycine?
        """
        
        if ref_aa.upper() != 'G' and alt_aa.upper() != 'G':
            return {'error': 'Not a glycine substitution'}
        
        protein_family = self.get_protein_family(gene)
        context = {
            'amino_acid': 'GLYCINE',
            'gene': gene,
            'position': position,
            'protein_family': protein_family,
            'substitution_type': 'LOSS' if ref_aa.upper() == 'G' else 'GAIN'
        }
        
        # COLLAGEN-SPECIFIC ANALYSIS
        if protein_family == 'COLLAGEN':
            context.update(self._analyze_collagen_glycine(position, uniprot_features))
        
        # ION CHANNEL-SPECIFIC ANALYSIS
        elif protein_family == 'ION_CHANNEL':
            context.update(self._analyze_ion_channel_glycine(position, uniprot_features))
        
        # FIBRILLIN-SPECIFIC ANALYSIS
        elif protein_family == 'FIBRILLIN':
            context.update(self._analyze_fibrillin_glycine(position, uniprot_features))
        
        # GENERAL STRUCTURAL ANALYSIS
        context.update(self._analyze_general_glycine_context(position, uniprot_features, conservation))
        
        return context
    
    def analyze_cysteine_context(self,
                                gene: str,
                                position: int,
                                ref_aa: str,
                                alt_aa: str,
                                uniprot_features: List[Dict[str, Any]] = None,
                                conservation: Optional[float] = None,
                                **kwargs) -> Dict[str, Any]:
        """
        🎯 Analyze biological context for Cysteine substitutions
        
        Returns context features for ML training:
        - disulfide_bond: Is this cysteine part of a disulfide bond?
        - metal_coordination: Is this coordinating metal ions?
        - catalytic_site: Is this in an active/catalytic site?
        - free_cysteine: Is this a free cysteine (potentially problematic)?
        - structural_importance: How critical is this for protein structure?
        """
        
        if ref_aa.upper() != 'C' and alt_aa.upper() != 'C':
            return {'error': 'Not a cysteine substitution'}
        
        protein_family = self.get_protein_family(gene)
        context = {
            'amino_acid': 'CYSTEINE',
            'gene': gene,
            'position': position,
            'protein_family': protein_family,
            'substitution_type': 'LOSS' if ref_aa.upper() == 'C' else 'GAIN'
        }
        
        # FIBRILLIN-SPECIFIC ANALYSIS (rich in disulfide bonds)
        if protein_family == 'FIBRILLIN':
            context.update(self._analyze_fibrillin_cysteine(position, uniprot_features))
        
        # ION CHANNEL-SPECIFIC ANALYSIS
        elif protein_family == 'ION_CHANNEL':
            context.update(self._analyze_ion_channel_cysteine(position, uniprot_features))
        
        # COLLAGEN-SPECIFIC ANALYSIS (rare cysteines, usually critical)
        elif protein_family == 'COLLAGEN':
            context.update(self._analyze_collagen_cysteine(position, uniprot_features))
        
        # GENERAL STRUCTURAL ANALYSIS
        context.update(self._analyze_general_cysteine_context(position, uniprot_features, conservation))
        
        return context
    
    def _analyze_collagen_glycine(self, position: int, features: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze glycine in collagen context - Gly-X-Y pattern detection"""
        
        context = {
            'collagen_gxy_pattern': False,
            'collagen_gxy_confidence': 0.0,
            'triple_helix_region': False,
            'collagen_cleavage_site': False
        }
        
        # Collagen Gly-X-Y pattern: Every 3rd residue should be glycine
        # Position % 3 == 1 means it's in the glycine position of Gly-X-Y
        if position % 3 == 1:
            context['collagen_gxy_pattern'] = True
            context['collagen_gxy_confidence'] = 0.95  # Very high confidence for collagen
        
        # Check for triple helix regions in features
        if features:
            for feature in features:
                feature_type = feature.get('category', '').lower()
                feature_desc = feature.get('description', '').lower()
                
                if 'triple' in feature_desc or 'helix' in feature_desc:
                    if feature.get('start', 0) <= position <= feature.get('end', 0):
                        context['triple_helix_region'] = True
                
                if 'cleavage' in feature_desc or 'propeptide' in feature_desc:
                    if feature.get('start', 0) <= position <= feature.get('end', 0):
                        context['collagen_cleavage_site'] = True
        
        return context
    
    def _analyze_ion_channel_glycine(self, position: int, features: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze glycine in ion channel context - gating vs. structural"""
        
        context = {
            'channel_gate_region': False,
            'transmembrane_domain': False,
            'channel_selectivity_filter': False,
            'channel_linker_region': False
        }
        
        if features:
            for feature in features:
                feature_type = feature.get('category', '').lower()
                feature_desc = feature.get('description', '').lower()
                
                if feature.get('start', 0) <= position <= feature.get('end', 0):
                    # Gate regions
                    if any(keyword in feature_desc for keyword in ['gate', 'gating', 'activation']):
                        context['channel_gate_region'] = True
                    
                    # Transmembrane domains
                    if 'transmembrane' in feature_desc or feature_type == 'transmembrane':
                        context['transmembrane_domain'] = True
                    
                    # Selectivity filter
                    if 'selectivity' in feature_desc or 'filter' in feature_desc:
                        context['channel_selectivity_filter'] = True
                    
                    # Linker regions (often flexible)
                    if 'linker' in feature_desc or 'loop' in feature_desc:
                        context['channel_linker_region'] = True
        
        return context
    
    def _analyze_fibrillin_glycine(self, position: int, features: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze glycine in fibrillin context - EGF-like domains"""
        
        context = {
            'egf_like_domain': False,
            'calcium_binding_egf': False,
            'fibrillin_unique_domain': False
        }
        
        if features:
            for feature in features:
                feature_desc = feature.get('description', '').lower()
                
                if feature.get('start', 0) <= position <= feature.get('end', 0):
                    if 'egf' in feature_desc:
                        context['egf_like_domain'] = True
                        
                        if 'calcium' in feature_desc or 'ca' in feature_desc:
                            context['calcium_binding_egf'] = True
                    
                    if 'fibrillin' in feature_desc and 'unique' in feature_desc:
                        context['fibrillin_unique_domain'] = True
        
        return context
    
    def _analyze_general_glycine_context(self, position: int, features: List[Dict[str, Any]], conservation: Optional[float]) -> Dict[str, Any]:
        """General glycine context analysis"""
        
        context = {
            'conservation_score': conservation if conservation is not None else 0.5,
            'active_site_proximity': False,
            'binding_site_proximity': False,
            'flexible_region': False,
            'tight_packing_region': False
        }
        
        if features:
            for feature in features:
                feature_type = feature.get('category', '').lower()
                feature_desc = feature.get('description', '').lower()
                
                # Check if position is within or near this feature
                start = feature.get('start', 0)
                end = feature.get('end', 0)
                
                if start <= position <= end:
                    # Active sites
                    if 'active' in feature_desc or feature_type == 'active_site':
                        context['active_site_proximity'] = True
                    
                    # Binding sites
                    if 'binding' in feature_desc or feature_type == 'binding_site':
                        context['binding_site_proximity'] = True
                    
                    # Flexible regions
                    if any(keyword in feature_desc for keyword in ['loop', 'linker', 'flexible', 'hinge']):
                        context['flexible_region'] = True
                
                # Proximity analysis (within 5 residues)
                elif abs(position - start) <= 5 or abs(position - end) <= 5:
                    if 'active' in feature_desc:
                        context['active_site_proximity'] = True
                    if 'binding' in feature_desc:
                        context['binding_site_proximity'] = True
        
        return context
    
    def _analyze_fibrillin_cysteine(self, position: int, features: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze cysteine in fibrillin context - disulfide bond rich"""
        
        context = {
            'disulfide_bond_predicted': True,  # Fibrillin is disulfide-rich
            'egf_domain_cysteine': False,
            'calcium_binding_cysteine': False,
            'structural_disulfide': True
        }
        
        if features:
            for feature in features:
                feature_desc = feature.get('description', '').lower()
                
                if feature.get('start', 0) <= position <= feature.get('end', 0):
                    if 'egf' in feature_desc:
                        context['egf_domain_cysteine'] = True
                        
                        if 'calcium' in feature_desc:
                            context['calcium_binding_cysteine'] = True
        
        return context
    
    def _analyze_ion_channel_cysteine(self, position: int, features: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze cysteine in ion channel context"""
        
        context = {
            'disulfide_bond_predicted': False,  # Less common in ion channels
            'metal_coordination_site': False,
            'channel_structure_critical': False,
            'extracellular_cysteine': False
        }
        
        if features:
            for feature in features:
                feature_type = feature.get('category', '').lower()
                feature_desc = feature.get('description', '').lower()
                
                if feature.get('start', 0) <= position <= feature.get('end', 0):
                    # Metal coordination
                    if any(keyword in feature_desc for keyword in ['zinc', 'metal', 'coordination']):
                        context['metal_coordination_site'] = True
                    
                    # Extracellular regions (more likely disulfide bonds)
                    if 'extracellular' in feature_desc or 'external' in feature_desc:
                        context['extracellular_cysteine'] = True
                        context['disulfide_bond_predicted'] = True
                    
                    # Structural importance
                    if any(keyword in feature_desc for keyword in ['structure', 'fold', 'domain']):
                        context['channel_structure_critical'] = True
        
        return context
    
    def _analyze_collagen_cysteine(self, position: int, features: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze cysteine in collagen context - rare but critical"""
        
        context = {
            'rare_collagen_cysteine': True,  # Cysteines are rare in collagen
            'cross_linking_site': False,
            'collagen_structure_critical': True,
            'disulfide_bond_predicted': True
        }
        
        if features:
            for feature in features:
                feature_desc = feature.get('description', '').lower()
                
                if feature.get('start', 0) <= position <= feature.get('end', 0):
                    if any(keyword in feature_desc for keyword in ['cross', 'link', 'bridge']):
                        context['cross_linking_site'] = True
        
        return context
    
    def _analyze_general_cysteine_context(self, position: int, features: List[Dict[str, Any]], conservation: Optional[float]) -> Dict[str, Any]:
        """General cysteine context analysis"""
        
        context = {
            'conservation_score': conservation if conservation is not None else 0.5,
            'catalytic_site_proximity': False,
            'binding_site_proximity': False,
            'free_cysteine_risk': False,
            'redox_sensitive_region': False
        }
        
        if features:
            for feature in features:
                feature_type = feature.get('category', '').lower()
                feature_desc = feature.get('description', '').lower()
                
                start = feature.get('start', 0)
                end = feature.get('end', 0)
                
                if start <= position <= end:
                    # Catalytic sites
                    if 'catalytic' in feature_desc or 'active' in feature_desc:
                        context['catalytic_site_proximity'] = True
                    
                    # Binding sites
                    if 'binding' in feature_desc:
                        context['binding_site_proximity'] = True
                    
                    # Redox sensitive
                    if any(keyword in feature_desc for keyword in ['redox', 'oxidation', 'reduction']):
                        context['redox_sensitive_region'] = True
                
                # Proximity analysis
                elif abs(position - start) <= 3 or abs(position - end) <= 3:
                    if 'catalytic' in feature_desc or 'active' in feature_desc:
                        context['catalytic_site_proximity'] = True
        
        # Free cysteine risk assessment
        # This would need more sophisticated analysis of nearby cysteines
        # For now, assume moderate risk
        context['free_cysteine_risk'] = 0.3
        
        return context
    
    def get_context_features(self, 
                           gene: str,
                           position: int,
                           ref_aa: str,
                           alt_aa: str,
                           **kwargs) -> Dict[str, Any]:
        """
        Main interface for getting context features for ML training
        
        Returns comprehensive context analysis for either Glycine or Cysteine variants
        """
        
        # Determine which amino acid we're analyzing
        if ref_aa.upper() == 'G' or alt_aa.upper() == 'G':
            return self.analyze_glycine_context(gene, position, ref_aa, alt_aa, **kwargs)
        elif ref_aa.upper() == 'C' or alt_aa.upper() == 'C':
            return self.analyze_cysteine_context(gene, position, ref_aa, alt_aa, **kwargs)
        else:
            return {'error': f'Not a Glycine or Cysteine variant: {ref_aa}{position}{alt_aa}'}


def test_gly_cys_context():
    """Test the Gly/Cys context analyzer with known variants"""
    
    analyzer = GlyCysContextAnalyzer()
    
    print("\n🧬 Testing Gly/Cys Context Analysis:")
    print("=" * 60)
    
    # Test collagen glycine (should be highly pathogenic)
    print("\n1. COLLAGEN GLYCINE TEST:")
    col_gly = analyzer.get_context_features('COL1A1', 893, 'G', 'A')
    print(f"COL1A1 p.G893A: {col_gly}")
    
    # Test fibrillin cysteine (should be disulfide critical)
    print("\n2. FIBRILLIN CYSTEINE TEST:")
    fbn_cys = analyzer.get_context_features('FBN1', 628, 'C', 'Y')
    print(f"FBN1 p.C628Y: {fbn_cys}")
    
    # Test ion channel glycine (should be context-dependent)
    print("\n3. ION CHANNEL GLYCINE TEST:")
    scn_gly = analyzer.get_context_features('SCN1A', 58, 'G', 'R')
    print(f"SCN1A p.G58R: {scn_gly}")
    
    # Test ion channel cysteine (should be variable)
    print("\n4. ION CHANNEL CYSTEINE TEST:")
    ryr_cys = analyzer.get_context_features('RYR1', 614, 'R', 'C')
    print(f"RYR1 p.R614C: {ryr_cys}")


if __name__ == "__main__":
    test_gly_cys_context()
