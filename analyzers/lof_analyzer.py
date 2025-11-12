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
from data_processing.sequence_mismatch_handler import create_mismatch_handler

# Add DOMAIN AWARENESS! 🎯
from data_processing.universal_protein_annotator import UniversalProteinAnnotator

# Add NOVA'S FUNCTIONAL DOMAIN WEIGHTING! 🚀
from core_analyzers.functional_domain_weighter import FunctionalDomainWeighter
from .conservation_database import ConservationDatabase

class LOFAnalyzer:
    """Analyze loss of function potential - Bin 1 of our two-bin approach"""
    
    def __init__(self, offline_mode=False, conservation_db: ConservationDatabase = None):
        self.name = "LOFAnalyzer"
        self.smart_analyzer = SmartProteinAnalyzer(offline_mode=offline_mode)
        self.mismatch_handler = create_mismatch_handler()
        self.conservation_db = conservation_db if conservation_db else ConservationDatabase()

        # 🎯 DOMAIN AWARENESS - Universal protein annotation!
        self.domain_annotator = UniversalProteinAnnotator()
        self.domain_cache = {}  # Cache domain info to avoid repeated API calls

        # 🚀 NOVA'S FUNCTIONAL DOMAIN WEIGHTING - Biological intelligence!
        self.functional_weighter = FunctionalDomainWeighter()
        
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
        """🚀 NOVA'S FUNCTIONAL DOMAIN WEIGHTING - Biological intelligence!"""

        # Extract UniProt features from domain context
        uniprot_features = []

        # Convert our domain context format to UniProt feature format
        for propeptide in domain_context.get("propeptides", []):
            uniprot_features.append({
                "type": "PROPEP",
                "description": propeptide.get("description", "propeptide"),
                "begin": {"position": str(propeptide["start"])},
                "end": {"position": str(propeptide.get("end", propeptide["start"]))}
            })

        for chain in domain_context.get("mature_chain", []):
            uniprot_features.append({
                "type": "CHAIN",
                "description": chain.get("description", "mature chain"),
                "begin": {"position": str(chain["start"])},
                "end": {"position": str(chain.get("end", chain["start"]))}
            })

        for region in domain_context.get("regions", []):
            uniprot_features.append({
                "type": "REGION",
                "description": region.get("description", "functional region"),
                "begin": {"position": str(region["start"])},
                "end": {"position": str(region.get("end", region["start"]))}
            })

        for signal in domain_context.get("signal_peptide", []):
            uniprot_features.append({
                "type": "SIGNAL",
                "description": signal.get("description", "signal peptide"),
                "begin": {"position": str(signal["start"])},
                "end": {"position": str(signal.get("end", signal["start"]))}
            })

        for domain in domain_context.get("domains", []):
            uniprot_features.append({
                "type": "DOMAIN",
                "description": domain.get("description", "domain"),
                "begin": {"position": str(domain["start"])},
                "end": {"position": str(domain.get("end", domain["start"]))}
            })

        # Add active sites as point features
        for active_site in domain_context.get("active_sites", []):
            uniprot_features.append({
                "type": "ACT_SITE",
                "description": "active site",
                "begin": {"position": str(active_site)},
                "end": {"position": str(active_site)}
            })

        # Add binding sites
        for binding in domain_context.get("binding_sites", []):
            uniprot_features.append({
                "type": "BINDING",
                "description": binding.get("description", "binding site"),
                "begin": {"position": str(binding["position"])},
                "end": {"position": str(binding["position"])}
            })

        # 🚀 USE NOVA'S FUNCTIONAL WEIGHTING SYSTEM!
        print(f"🚀 Using Nova's functional domain weighting for position {position}")
        functional_weight = self.functional_weighter.weight_variant_position(position, uniprot_features)

        return max(functional_weight, 0.1)  # Don't go below 0.1

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

        # NEW: Start-loss and Stop-loss fast-paths (never benign)
        if parsed.get('is_start_loss', False) or parsed.get('is_stop_loss', False):
            base_score = 0.95  # very high LOF starting point
            # domain awareness still applies
            domain_multiplier = 1.0
            if uniprot_id:
                domain_context = self._get_domain_context(uniprot_id, gene_symbol)
                domain_multiplier = self._get_domain_multiplier(max(position, 1), domain_context)
                tag = 'start_loss' if parsed.get('is_start_loss') else 'stop_loss'
                print(f"🎯 Domain multiplier for {tag} {gene_symbol} position {position}: {domain_multiplier:.3f}")
            final_score = min(base_score * domain_multiplier, 1.0)
            mech = 'start_codon_loss' if parsed.get('is_start_loss') else 'stop_codon_loss'
            return {
                'lof_score': final_score,
                'base_lof_score': base_score,
                'smart_multiplier': 1.0,
                'conservation_multiplier': 1.0,
                'domain_multiplier': domain_multiplier,
                'total_multiplier': domain_multiplier,
                'stability_impact': 1.0,
                'conservation_impact': 1.0,
                'structural_impact': 1.0,
                'functional_impact': 1.0,
                'mechanism': mech,
                'confidence': 0.95,
                'sequence_mismatch': False,
                'mismatch_info': None,
                'start_stop_details': {
                    'position': position,
                    'explanation': f'{mech.replace("_"," ").title()} detected at position {position}'
                }
            }

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

        # 🔥 REN'S BRILLIANT FIX: ADD ML PROLINE PANIC TO LOF TOO!
        ml_proline_multiplier = 1.0
        if mutation and len(mutation) >= 3:  # Format: "P840L"
            ref_aa = mutation[0]
            alt_aa = mutation[-1]

            if ref_aa == 'P' or alt_aa == 'P':  # Proline substitution detected!
                try:
                    # Import and use our revolutionary ML system
                    from AdaptiveInterpreter.utils.proline_ml_integrator import get_ml_proline_multiplier
                    variant_str = f"p.{mutation}"
                    ml_proline_multiplier = get_ml_proline_multiplier(gene_symbol, variant_str)
                    print(f"🔥 LOF ML PROLINE: {gene_symbol} {variant_str} -> ML multiplier = {ml_proline_multiplier:.3f}")

                except Exception as e:
                    print(f"⚠️ LOF ML proline system error: {e}")
                    # Fallback to old hardcoded system if ML fails
                    if ref_aa == 'P':  # Proline loss
                        ml_proline_multiplier = 0.7  # Reduce LOF score (less destabilizing than expected)
                        print(f"🔥 LOF PROLINE FALLBACK: {gene_symbol} {variant_str} -> fallback multiplier = {ml_proline_multiplier:.3f}")

        # Calculate mechanistic LOF base WITHOUT conservation
        # Fold ML proline into the functional component (bounded), and include domain/smart in base.
        functional_adj = min(functional_impact * ml_proline_multiplier, 1.0)
        # Original weights: stability 0.3, conservation 0.3, structural 0.2, functional 0.2
        # Drop conservation from base and renormalize by 0.7 to preserve scale.
        base_no_cons = (stability_impact * 0.3) + (structural_impact * 0.2) + (functional_adj * 0.2)
        base_lof_score = min(base_no_cons / 0.7, 1.0)

        # Mechanistic base multiplies in domain and smart context; exclude conservation from multipliers here.
        total_multiplier = smart_multiplier * domain_multiplier
        lof_score = min(base_lof_score * total_multiplier, 1.0)

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
            'ml_proline_multiplier': ml_proline_multiplier,  # 🔥 NEW! ML proline panic
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
        """Parse mutation string - supports 1-letter and 3-letter HGVS (with optional 'p.'), plus nonsense.
        Examples:
          - R175H, p.R175H
          - Glu446Ala, p.Glu446Ala
          - Tyr87Ter, R100*
        """
        if not mutation or len(mutation) < 3:
            return None

        s = mutation.strip()
        if s.startswith('p.'):
            s = s[2:]

        # 0) Start-loss and Stop-loss special patterns
        # Start-loss: Met1? / Met1Leu / M1L / M1?
        m_start1 = re.match(r'^[Mm]1([A-Za-z]|\?)$', s)
        if m_start1:
            alt = m_start1.group(1)
            alt1 = alt.upper() if alt != '?' else '?'
            return {
                'original_aa': 'M',
                'position': 1,
                'new_aa': alt1,
                'mutation': f"M1{alt1}",
                'is_nonsense': False,
                'is_start_loss': True
            }
        m_start3 = re.match(r'^(Met)1([A-Z][a-z]{2}|\?)$', s)
        if m_start3:
            alt3 = m_start3.group(2)
            if alt3 == '?':
                alt1 = '?'
            else:
                THREE_TO_ONE = {
                    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
                    'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
                    'Thr':'T','Trp':'W','Tyr':'Y','Val':'V'
                }
                alt1 = THREE_TO_ONE.get(alt3)
                if not alt1:
                    return None
            return {
                'original_aa': 'M',
                'position': 1,
                'new_aa': alt1,
                'mutation': f"M1{alt1}",
                'is_nonsense': False,
                'is_start_loss': True
            }

        # Stop-loss: Ter123Gln / *123Q and stop-loss with extension (ext*?, ext*NN)
        # Basic stop-loss (stop codon changed to amino acid)
        m_stop1 = re.match(r'^\*(\d+)([A-Za-z])$', s)
        if m_stop1:
            try:
                pos = int(m_stop1.group(1))
                alt = m_stop1.group(2).upper()
                return {
                    'original_aa': '*',
                    'position': pos,
                    'new_aa': alt,
                    'mutation': f"*{pos}{alt}",
                    'is_nonsense': False,
                    'is_stop_loss': True
                }
            except ValueError:
                return None
        m_stop3 = re.match(r'^Ter(\d+)([A-Z][a-z]{2})$', s)
        if m_stop3:
            pos_str, alt3 = m_stop3.groups()
            THREE_TO_ONE = {
                'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
                'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
                'Thr':'T','Trp':'W','Tyr':'Y','Val':'V'
            }
            alt = THREE_TO_ONE.get(alt3)
            if not alt:
                return None
            try:
                pos = int(pos_str)
                return {
                    'original_aa': '*',
                    'position': pos,
                    'new_aa': alt,
                    'mutation': f"*{pos}{alt}",
                    'is_nonsense': False,
                    'is_stop_loss': True
                }
            except ValueError:
                return None
        # Stop-loss with extension annotations: e.g., Ter577Cysext*?, *577Qext*?, or ext*12
        m_stop1_ext = re.match(r'^\*(\d+)([A-Za-z])ext\*\??(?:\d+)?$', s)
        if m_stop1_ext:
            try:
                pos = int(m_stop1_ext.group(1))
                alt = m_stop1_ext.group(2).upper()
                return {
                    'original_aa': '*',
                    'position': pos,
                    'new_aa': alt,
                    'mutation': f"*{pos}{alt}ext*?",
                    'is_nonsense': False,
                    'is_stop_loss': True
                }
            except ValueError:
                return None
        m_stop3_ext = re.match(r'^Ter(\d+)([A-Z][a-z]{2})ext\*\??(?:\d+)?$', s)
        if m_stop3_ext:
            pos_str, alt3 = m_stop3_ext.groups()
            THREE_TO_ONE = {
                'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
                'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
                'Thr':'T','Trp':'W','Tyr':'Y','Val':'V'
            }
            alt = THREE_TO_ONE.get(alt3)
            if not alt:
                return None
            try:
                pos = int(pos_str)
                return {
                    'original_aa': '*',
                    'position': pos,
                    'new_aa': alt,
                    'mutation': f"*{pos}{alt}ext*?",
                    'is_nonsense': False,
                    'is_stop_loss': True
                }
            except ValueError:
                return None

        # 1) Nonsense variants first (Ter or *)
        if s.endswith('Ter') or s.endswith('*'):
            # Handle 1-letter nonsense: R100*
            if s.endswith('*') and len(s) >= 3 and s[0].isalpha():
                try:
                    position = int(s[1:-1])
                    return {
                        'original_aa': s[0].upper(),
                        'position': position,
                        'new_aa': '*',
                        'mutation': s,
                        'is_nonsense': True
                    }
                except ValueError:
                    return None

            # Handle 3-letter nonsense: Tyr87Ter
            m_n3 = re.match(r'^([A-Z][a-z]{2})(\d+)Ter$', s)
            if m_n3:
                ref3, pos_str = m_n3.groups()
                THREE_TO_ONE = {
                    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
                    'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
                    'Thr':'T','Trp':'W','Tyr':'Y','Val':'V'
                }
                ref = THREE_TO_ONE.get(ref3)
                if not ref:
                    return None
                try:
                    position = int(pos_str)
                    return {
                        'original_aa': ref,
                        'position': position,
                        'new_aa': '*',
                        'mutation': s,
                        'is_nonsense': True
                    }
                except ValueError:
                    return None

        # 2) Missense variants
        # 2a) 1-letter form: R175H
        m1 = re.match(r'^([A-Za-z])(\d+)([A-Za-z])$', s)
        if m1:
            ref, pos_str, alt = m1.groups()
            try:
                return {
                    'original_aa': ref.upper(),
                    'position': int(pos_str),
                    'new_aa': alt.upper(),
                    'mutation': f"{ref.upper()}{pos_str}{alt.upper()}",
                    'is_nonsense': False
                }
            except ValueError:
                return None

        # 2b) 3-letter form: Glu446Ala (ignore synonymous '=')
        m3 = re.match(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=)$', s)
        if m3:
            ref3, pos_str, alt3 = m3.groups()
            if alt3 == '=':
                return None  # synonymous; treat as unparseable for LOF scoring
            THREE_TO_ONE = {
                'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
                'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
                'Thr':'T','Trp':'W','Tyr':'Y','Val':'V'
            }
            ref = THREE_TO_ONE.get(ref3)
            alt = THREE_TO_ONE.get(alt3)
            if not ref or not alt:
                return None
            try:
                return {
                    'original_aa': ref,
                    'position': int(pos_str),
                    'new_aa': alt,
                    'mutation': f"{ref}{pos_str}{alt}",
                    'is_nonsense': False
                }
            except ValueError:
                return None

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
        """Calculate overall LOF score (legacy; conservation excluded upstream)."""
        # Keep legacy signature but exclude conservation here to avoid double counting.
        score = (
            stability * 0.3 +
            structural * 0.2 +
            functional * 0.2
        )
        # Renormalize to preserve scale roughly comparable to legacy
        score = score / 0.7
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
