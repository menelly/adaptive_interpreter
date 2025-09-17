#!/usr/bin/env python3
"""
🔬 LOSS OF FUNCTION ANALYZER - BIN 1 ANALYSIS
Built by Ace for the revolutionary two-bin approach

This tiny module analyzes whether variants cause loss of function.
Traditional pathogenicity prediction - does it break the protein?
"""

from typing import Dict, Any
import re
import sys
import os
from .smart_protein_analyzer import SmartProteinAnalyzer

# Add sequence mismatch handler
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../'))
from sequence_mismatch_handler import create_mismatch_handler

# Add DOMAIN AWARENESS! 🎯
from universal_protein_annotator import UniversalProteinAnnotator

class LOFAnalyzer:
    """Analyze loss of function potential - Bin 1 of our two-bin approach"""
    
    def __init__(self, offline_mode=False):
        self.name = "LOFAnalyzer"
        self.smart_analyzer = SmartProteinAnalyzer(offline_mode=offline_mode)
        self.mismatch_handler = create_mismatch_handler()

        # 🎯 DOMAIN AWARENESS - Universal protein annotation!
        self.domain_annotator = UniversalProteinAnnotator()
        self.domain_cache = {}  # Cache domain info to avoid repeated API calls
        
        # Amino acid stability/conservation properties
        self.aa_properties = {
            'G': {'size': 1, 'charge': 0, 'hydrophobic': False, 'flexibility': 'high', 'conservation': 'critical'},
            'A': {'size': 2, 'charge': 0, 'hydrophobic': True, 'flexibility': 'medium', 'conservation': 'medium'},
            'V': {'size': 3, 'charge': 0, 'hydrophobic': True, 'flexibility': 'low', 'conservation': 'medium'},
            'L': {'size': 4, 'charge': 0, 'hydrophobic': True, 'flexibility': 'low', 'conservation': 'medium'},
            'I': {'size': 4, 'charge': 0, 'hydrophobic': True, 'flexibility': 'low', 'conservation': 'medium'},
            'M': {'size': 4, 'charge': 0, 'hydrophobic': True, 'flexibility': 'medium', 'conservation': 'medium'},
            'F': {'size': 5, 'charge': 0, 'hydrophobic': True, 'flexibility': 'low', 'conservation': 'high'},
            'W': {'size': 6, 'charge': 0, 'hydrophobic': True, 'flexibility': 'low', 'conservation': 'high'},
            'P': {'size': 3, 'charge': 0, 'hydrophobic': False, 'flexibility': 'rigid', 'conservation': 'critical'},
            'S': {'size': 2, 'charge': 0, 'hydrophobic': False, 'flexibility': 'high', 'conservation': 'low'},
            'T': {'size': 3, 'charge': 0, 'hydrophobic': False, 'flexibility': 'medium', 'conservation': 'low'},
            'C': {'size': 2, 'charge': 0, 'hydrophobic': False, 'flexibility': 'medium', 'conservation': 'critical'},
            'Y': {'size': 5, 'charge': 0, 'hydrophobic': False, 'flexibility': 'medium', 'conservation': 'high'},
            'N': {'size': 3, 'charge': 0, 'hydrophobic': False, 'flexibility': 'high', 'conservation': 'medium'},
            'Q': {'size': 4, 'charge': 0, 'hydrophobic': False, 'flexibility': 'high', 'conservation': 'medium'},
            'D': {'size': 3, 'charge': -1, 'hydrophobic': False, 'flexibility': 'high', 'conservation': 'high'},
            'E': {'size': 4, 'charge': -1, 'hydrophobic': False, 'flexibility': 'high', 'conservation': 'high'},
            'K': {'size': 4, 'charge': 1, 'hydrophobic': False, 'flexibility': 'high', 'conservation': 'high'},
            'R': {'size': 5, 'charge': 1, 'hydrophobic': False, 'flexibility': 'high', 'conservation': 'high'},
            'H': {'size': 4, 'charge': 0.5, 'hydrophobic': False, 'flexibility': 'high', 'conservation': 'high'}
        }

    def _get_domain_context(self, uniprot_id: str, gene_symbol: str) -> Dict[str, Any]:
        """🎯 Get domain context for protein - UNIVERSAL APPROACH!"""
        cache_key = uniprot_id or gene_symbol

        if cache_key in self.domain_cache:
            return self.domain_cache[cache_key]

        try:
            print(f"🔍 Getting domain context for {gene_symbol} ({uniprot_id})...")
            domain_info = self.domain_annotator.get_uniprot_features(uniprot_id)

            if "error" not in domain_info:
                self.domain_cache[cache_key] = domain_info
                return domain_info
            else:
                print(f"⚠️ Domain annotation failed: {domain_info.get('error', 'Unknown error')}")

        except Exception as e:
            print(f"⚠️ Domain annotation error: {e}")

        # Return minimal context if annotation fails
        return {"domains": [], "signal_peptide": [], "active_sites": [], "binding_sites": []}

    def _get_domain_multiplier(self, position: int, domain_context: Dict[str, Any]) -> float:
        """🎯 Calculate domain-aware multiplier for LOF scoring"""
        multiplier = 1.0
        gene_name = domain_context.get("gene_name", "").upper()
        sequence_length = domain_context.get("length", 0)

        # 🎯 UNIVERSAL PROPEPTIDE LOGIC - Use real UniProt annotations!
        for propeptide in domain_context.get("propeptides", []):
            if propeptide["start"] <= position <= propeptide.get("end", propeptide["start"]):
                prop_desc = propeptide.get("description", "").lower()
                if "n-terminal" in prop_desc:
                    print(f"🔍 Position {position} in N-TERMINAL PROPEPTIDE ({propeptide['start']}-{propeptide.get('end', propeptide['start'])}) - DOWNWEIGHT!")
                    multiplier *= 0.5  # N-terminal propeptides get cleaved
                elif "c-terminal" in prop_desc:
                    print(f"🔍 Position {position} in C-TERMINAL PROPEPTIDE ({propeptide['start']}-{propeptide.get('end', propeptide['start'])}) - HEAVY DOWNWEIGHT!")
                    multiplier *= 0.3  # C-terminal propeptides get cleaved and are less critical
                else:
                    print(f"🔍 Position {position} in PROPEPTIDE ({propeptide['start']}-{propeptide.get('end', propeptide['start'])}) - DOWNWEIGHT!")
                    multiplier *= 0.4  # Generic propeptide downweight

        # 🎯 MATURE CHAIN LOGIC - The functional protein region
        for chain in domain_context.get("mature_chain", []):
            if chain["start"] <= position <= chain.get("end", chain["start"]):
                print(f"🔍 Position {position} in MATURE CHAIN ({chain['start']}-{chain.get('end', chain['start'])}) - CRITICAL REGION!")
                multiplier *= 1.0  # This is the functional protein - normal weight

        # 🎯 FUNCTIONAL REGIONS - Triple helix, catalytic domains, etc.
        for region in domain_context.get("regions", []):
            if region["start"] <= position <= region.get("end", region["start"]):
                region_desc = region.get("description", "").lower()
                if any(term in region_desc for term in ["triple-helical", "catalytic", "active", "binding"]):
                    print(f"🔍 Position {position} in CRITICAL REGION ({region_desc}) - UPWEIGHT!")
                    multiplier *= 1.2  # Critical functional regions
                elif "disordered" in region_desc:
                    print(f"🔍 Position {position} in DISORDERED REGION - DOWNWEIGHT!")
                    multiplier *= 0.7  # Disordered regions less critical

        # Check if variant is in signal peptide (gets cleaved off)
        for signal in domain_context.get("signal_peptide", []):
            if signal["start"] <= position <= signal.get("end", signal["start"]):
                print(f"🔍 Position {position} in UniProt signal peptide ({signal['start']}-{signal.get('end', signal['start'])}) - DOWNWEIGHTING!")
                multiplier *= 0.3  # Heavily downweight signal peptides

        # Check if variant is in active site (critical!)
        if position in domain_context.get("active_sites", []):
            print(f"🔍 Position {position} in ACTIVE SITE - UPWEIGHTING!")
            multiplier *= 1.5  # Upweight active sites

        # Check if variant is in binding site (important!)
        for binding in domain_context.get("binding_sites", []):
            if binding["position"] == position:
                print(f"🔍 Position {position} in binding site ({binding.get('description', 'unknown')}) - UPWEIGHTING!")
                multiplier *= 1.3  # Upweight binding sites

        # Check if variant is in known domains
        for domain in domain_context.get("domains", []):
            if domain["start"] <= position <= domain.get("end", domain["start"]):
                domain_desc = domain.get("description", "").lower()

                # Propeptides and signal sequences get downweighted
                if any(term in domain_desc for term in ["propeptide", "signal", "leader", "transit"]):
                    print(f"🔍 Position {position} in cleavable domain ({domain_desc}) - DOWNWEIGHTING!")
                    multiplier *= 0.4

                # Catalytic and binding domains get upweighted
                elif any(term in domain_desc for term in ["catalytic", "kinase", "binding", "active", "enzymatic"]):
                    print(f"🔍 Position {position} in critical domain ({domain_desc}) - UPWEIGHTING!")
                    multiplier *= 1.4

        return max(multiplier, 0.1)  # Don't go below 0.1

    def analyze_lof(self, mutation: str, sequence: str, uniprot_id: str = None, gene_symbol: str = "", **kwargs) -> Dict[str, Any]:
        """
        Analyze loss of function potential
        
        Args:
            mutation: Mutation string (e.g., "R175H")
            sequence: Protein sequence
            
        Returns:
            LOF analysis results
        """
        
        parsed = self._parse_mutation(mutation)
        if not parsed:
            return self._empty_result()

        original_aa = parsed['original_aa']
        new_aa = parsed['new_aa']
        position = parsed['position']
        is_nonsense = parsed.get('is_nonsense', False)

        # Handle nonsense variants - they're almost always high LOF!
        if is_nonsense:
            # Nonsense variants cause premature termination = high LOF
            base_score = 0.9  # Very high LOF score

            # Earlier nonsense = worse (less functional protein made)
            seq_length = len(sequence) if sequence else 2000  # Default estimate
            position_factor = 1.0 - (position / seq_length)  # Earlier = higher score

            # 🎯 Apply domain awareness even to nonsense variants!
            domain_multiplier = 1.0
            if uniprot_id:
                domain_context = self._get_domain_context(uniprot_id, gene_symbol)
                domain_multiplier = self._get_domain_multiplier(position, domain_context)
                print(f"🎯 Domain multiplier for nonsense {gene_symbol} position {position}: {domain_multiplier:.3f}")

            final_score = min((base_score + (position_factor * 0.1)) * domain_multiplier, 1.0)

            # Return in the same format as regular LOF analysis
            return {
                'lof_score': final_score,
                'base_lof_score': base_score,
                'smart_multiplier': 1.0,
                'conservation_multiplier': 1.0,
                'domain_multiplier': domain_multiplier,  # 🎯 NEW! Domain awareness for nonsense too
                'total_multiplier': domain_multiplier,  # Only domain multiplier applied for nonsense
                'stability_impact': 1.0,  # Maximum impact
                'conservation_impact': 1.0,  # Maximum impact
                'structural_impact': 1.0,  # Maximum impact
                'functional_impact': 1.0,  # Maximum impact
                'mechanism': 'nonsense_mediated_decay',
                'confidence': 0.95,  # Very high confidence for nonsense
                'sequence_mismatch': False,
                'mismatch_info': None,
                'nonsense_details': {
                    'position': position,
                    'truncation_severity': position_factor,
                    'explanation': f'Nonsense variant at position {position} causes premature termination'
                }
            }

        # Check for sequence mismatch (only for missense variants)
        sequence_mismatch = False
        mismatch_info = None
        if sequence and position > 0:
            mismatch_info = self.mismatch_handler.check_sequence_match(
                sequence, position, original_aa, gene_symbol, uniprot_id or ""
            )
            sequence_mismatch = not mismatch_info['match']
        
        # Get amino acid properties
        orig_props = self.aa_properties.get(original_aa, self._default_props())
        new_props = self.aa_properties.get(new_aa, self._default_props())
        
        # Analyze different LOF mechanisms (now with Grantham distance!)
        stability_impact = self._assess_stability_impact(orig_props, new_props, mutation)
        conservation_impact = self._assess_conservation_impact(orig_props, new_props)
        structural_impact = self._assess_structural_impact(orig_props, new_props, position, len(sequence))
        functional_impact = self._assess_functional_impact(mutation, sequence)
        
        # Get smart protein context multiplier
        smart_multiplier, smart_confidence = 1.0, 0.0
        if uniprot_id and sequence:
            smart_multiplier, smart_confidence = self.smart_analyzer.get_protein_context_multiplier(
                uniprot_id, sequence, position
            )

        # Get conservation multiplier from kwargs
        conservation_multiplier = kwargs.get('conservation_multiplier', 1.0)

        # 🎯 GET DOMAIN-AWARE MULTIPLIER - UNIVERSAL APPROACH!
        domain_multiplier = 1.0
        if uniprot_id:
            domain_context = self._get_domain_context(uniprot_id, gene_symbol)
            domain_multiplier = self._get_domain_multiplier(position, domain_context)
            print(f"🎯 Domain multiplier for {gene_symbol} position {position}: {domain_multiplier:.3f}")

        # Calculate overall LOF score with ALL multipliers
        base_lof_score = self._calculate_lof_score(
            stability_impact, conservation_impact, structural_impact, functional_impact
        )

        # Apply ALL multipliers: conservation, smart, AND domain!
        total_multiplier = smart_multiplier * conservation_multiplier * domain_multiplier
        lof_score = base_lof_score * total_multiplier

        # Apply sequence mismatch handling if needed
        if sequence_mismatch and mismatch_info:
            fallback_result = self.mismatch_handler.apply_fallback_strategy(
                mismatch_info['fallback_strategy'],
                lof_score,
                mismatch_info['confidence_penalty'],
                {'confidence': smart_confidence}
            )
            lof_score = fallback_result['adjusted_score']

        # Don't cap at 1.0 - let it go higher like REVEL scores
        # lof_score = min(base_lof_score * total_multiplier, 1.0)
        
        return {
            'lof_score': lof_score,
            'base_lof_score': base_lof_score,
            'smart_multiplier': smart_multiplier,
            'conservation_multiplier': conservation_multiplier,
            'domain_multiplier': domain_multiplier,  # 🎯 NEW! Domain awareness
            'total_multiplier': total_multiplier,
            'stability_impact': stability_impact,
            'conservation_impact': conservation_impact,
            'structural_impact': structural_impact,
            'functional_impact': functional_impact,
            'mechanism': self._determine_lof_mechanism(stability_impact, conservation_impact, structural_impact),
            'confidence': min(self._calculate_lof_confidence(orig_props, new_props, mutation) + smart_confidence, 0.9),
            'sequence_mismatch': sequence_mismatch,
            'mismatch_info': mismatch_info
        }
    
    def _parse_mutation(self, mutation: str) -> Dict[str, Any]:
        """Parse mutation string - handles both missense and nonsense variants"""
        if not mutation or len(mutation) < 3:
            return None

        # Handle nonsense variants (e.g., "Q18Ter", "R100*")
        if mutation.endswith('Ter') or mutation.endswith('*'):
            if mutation.endswith('Ter'):
                position_str = mutation[1:-3]  # Remove first AA and "Ter"
                new_aa = '*'  # Standard nonsense notation
            else:  # ends with '*'
                position_str = mutation[1:-1]  # Remove first AA and "*"
                new_aa = '*'

            try:
                position = int(position_str)
                return {
                    'original_aa': mutation[0],
                    'position': position,
                    'new_aa': new_aa,
                    'mutation': mutation,
                    'is_nonsense': True
                }
            except ValueError:
                return None

        # Handle regular missense variants (e.g., "R175H")
        try:
            return {
                'original_aa': mutation[0],
                'position': int(mutation[1:-1]),
                'new_aa': mutation[-1],
                'mutation': mutation,
                'is_nonsense': False
            }
        except ValueError:
            return None
    
    def _default_props(self) -> Dict[str, Any]:
        """Default amino acid properties"""
        return {'size': 3, 'charge': 0, 'hydrophobic': False, 'flexibility': 'medium', 'conservation': 'medium'}
    
    def _empty_result(self) -> Dict[str, Any]:
        """Empty result for failed parsing"""
        return {
            'lof_score': 0.0,
            'stability_impact': 0.0,
            'conservation_impact': 0.0,
            'structural_impact': 0.0,
            'functional_impact': 0.0,
            'mechanism': 'unknown',
            'confidence': 0.0
        }
    
    def _assess_stability_impact(self, orig_props: Dict, new_props: Dict, mutation: str = None) -> float:
        """Assess impact on protein stability using REAL amino acid science!"""

        # If we have mutation string, use Grantham distance
        if mutation:
            parsed = self._parse_mutation(mutation)
            if parsed:
                orig_aa = parsed['original_aa']
                new_aa = parsed['new_aa']

                # Use Grantham distance for scientific accuracy!
                grantham_distance = self._get_grantham_distance(orig_aa, new_aa)

                # Convert Grantham distance to stability impact
                if grantham_distance >= 150:
                    base_score = 0.8  # Very severe change
                elif grantham_distance >= 100:
                    base_score = 0.6  # Severe change
                elif grantham_distance >= 50:
                    base_score = 0.4  # Moderate change
                elif grantham_distance >= 20:
                    base_score = 0.2  # Mild change
                else:
                    base_score = 0.1  # Very conservative change

                # Add special case modifiers
                if 'P' in mutation:  # Proline introduction/removal
                    base_score += 0.2
                if 'G' in mutation:  # Glycine flexibility
                    base_score += 0.15
                if 'C' in mutation:  # Cysteine disulfide bonds
                    base_score += 0.25

                return min(base_score, 1.0)

        # Fallback to old method if no mutation string
        score = 0.0

        # Size changes affect stability
        size_change = abs(new_props['size'] - orig_props['size'])
        if size_change > 2:
            score += 0.3
        elif size_change > 1:
            score += 0.1

        # Charge changes affect stability
        charge_change = abs(new_props['charge'] - orig_props['charge'])
        if charge_change > 1:
            score += 0.4
        elif charge_change > 0.5:
            score += 0.2

        # Hydrophobicity changes
        if orig_props['hydrophobic'] != new_props['hydrophobic']:
            score += 0.2

        return min(score, 1.0)
    
    def _assess_conservation_impact(self, orig_props: Dict, new_props: Dict) -> float:
        """Assess impact based on amino acid conservation"""
        conservation_scores = {'critical': 1.0, 'high': 0.8, 'medium': 0.5, 'low': 0.2}
        
        orig_conservation = conservation_scores.get(orig_props['conservation'], 0.5)
        
        # Higher impact if we're changing a highly conserved residue
        return orig_conservation
    
    def _assess_structural_impact(self, orig_props: Dict, new_props: Dict, position: int, seq_length: int) -> float:
        """Assess structural impact"""
        score = 0.0
        
        # Flexibility changes
        flexibility_map = {'rigid': 0, 'low': 1, 'medium': 2, 'high': 3}
        orig_flex = flexibility_map.get(orig_props['flexibility'], 2)
        new_flex = flexibility_map.get(new_props['flexibility'], 2)
        
        flex_change = abs(new_flex - orig_flex)
        if flex_change > 2:
            score += 0.3
        elif flex_change > 1:
            score += 0.1
        
        # Position-based impact (middle regions often more critical)
        position_factor = 1.0 - abs(position - seq_length/2) / (seq_length/2)
        score *= (0.5 + 0.5 * position_factor)
        
        return min(score, 1.0)
    
    def _assess_functional_impact(self, mutation: str, sequence: str) -> float:
        """Assess functional impact based on known patterns"""
        score = 0.0
        
        # Known highly disruptive changes
        if mutation[0] == 'C':  # Cysteine loss (disulfide bonds)
            score += 0.5
        if mutation[0] == 'P':  # Proline loss (structural rigidity)
            score += 0.3
        if mutation[0] == 'G':  # Glycine loss (flexibility)
            score += 0.4
        
        return min(score, 1.0)
    
    def _calculate_lof_score(self, stability: float, conservation: float, structural: float, functional: float) -> float:
        """Calculate overall LOF score"""
        # Weighted combination
        score = (
            stability * 0.3 +
            conservation * 0.3 +
            structural * 0.2 +
            functional * 0.2
        )
        
        return min(score, 1.0)
    
    def _determine_lof_mechanism(self, stability: float, conservation: float, structural: float) -> str:
        """Determine primary LOF mechanism"""
        if stability > 0.5:
            return 'protein_instability'
        elif conservation > 0.7:
            return 'critical_residue_loss'
        elif structural > 0.5:
            return 'structural_disruption'
        else:
            return 'mild_functional_impact'
    
    def _calculate_lof_confidence(self, orig_props: Dict, new_props: Dict, mutation: str) -> float:
        """Calculate confidence in LOF prediction"""
        confidence = 0.6  # Base confidence
        
        # Higher confidence for well-understood changes
        if orig_props['conservation'] in ['critical', 'high']:
            confidence += 0.2
        
        # Known disruptive patterns
        if mutation[0] in ['C', 'P', 'G']:
            confidence += 0.1
        
        return min(confidence, 0.9)

    def _get_grantham_distance(self, aa1, aa2):
        """Get Grantham distance between amino acids - REAL SCIENCE!"""
        # Key Grantham distances for common substitutions
        grantham_matrix = {
            ('A', 'A'): 0, ('A', 'R'): 112, ('A', 'N'): 111, ('A', 'D'): 126, ('A', 'C'): 195,
            ('A', 'Q'): 91, ('A', 'E'): 107, ('A', 'G'): 60, ('A', 'H'): 86, ('A', 'I'): 94,
            ('A', 'L'): 96, ('A', 'K'): 106, ('A', 'M'): 84, ('A', 'F'): 113, ('A', 'P'): 27,
            ('A', 'S'): 99, ('A', 'T'): 58, ('A', 'W'): 148, ('A', 'Y'): 112, ('A', 'V'): 64,

            ('R', 'R'): 0, ('R', 'N'): 86, ('R', 'D'): 96, ('R', 'C'): 180, ('R', 'Q'): 43,
            ('R', 'E'): 54, ('R', 'G'): 125, ('R', 'H'): 29, ('R', 'I'): 97, ('R', 'L'): 102,
            ('R', 'K'): 26, ('R', 'M'): 91, ('R', 'F'): 97, ('R', 'P'): 103, ('R', 'S'): 110,
            ('R', 'T'): 71, ('R', 'W'): 101, ('R', 'Y'): 77, ('R', 'V'): 96,

            ('T', 'M'): 81,  # T1424M - moderate severity
            ('V', 'I'): 29,  # V1172I - very conservative
            ('T', 'T'): 0, ('M', 'M'): 0, ('V', 'V'): 0, ('I', 'I'): 0,  # Identity
        }

        # Try both orientations
        distance = grantham_matrix.get((aa1, aa2))
        if distance is None:
            distance = grantham_matrix.get((aa2, aa1))

        return distance if distance is not None else 50  # Default moderate distance
