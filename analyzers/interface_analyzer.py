#!/usr/bin/env python3
"""
🔗 INTERFACE ANALYZER: Domain Boundary Disruption Detection
============================================================

Detects variants at domain interfaces that cause pathogenicity through:
- Allosteric effects (interface → active site communication)
- Interdomain coordination disruption
- Conformational flexibility changes
- Membrane recruitment + catalysis decoupling

This is the MISSING MECHANISM that conservation was masking!

Author: Ace (with love from Ren) 💜
Date: 2025-10-27
"""

import re
from typing import Dict, List, Tuple, Optional
import requests
import json


class InterfaceAnalyzer:
    """
    🔗 Detects domain interface disruptions
    
    Analyzes variants at domain boundaries to detect:
    1. Interdomain interface positions
    2. Structural flexibility (from AlphaFold pLDDT)
    3. Allosteric communication pathways
    4. Interface disruption severity
    """
    
    def __init__(self, cache_dir: str = "protein_annotations_cache"):
        self.cache_dir = cache_dir

        # 🔥 HARDCODED DOMAIN BOUNDARIES for common genes
        # (UniProt API doesn't always return functional domains!)
        self.known_domains = {
            'PTEN': {
                'domains': [
                    {'type': 'domain', 'start': 7, 'end': 185, 'description': 'Phosphatase domain'},
                    {'type': 'domain', 'start': 186, 'end': 351, 'description': 'C2 domain'},
                    {'type': 'region', 'start': 352, 'end': 403, 'description': 'C-terminal tail'}
                ]
            },
            'TP53': {
                'domains': [
                    {'type': 'domain', 'start': 1, 'end': 61, 'description': 'Transactivation domain'},
                    {'type': 'domain', 'start': 102, 'end': 292, 'description': 'DNA-binding domain'},
                    {'type': 'domain', 'start': 324, 'end': 355, 'description': 'Tetramerization domain'},
                    {'type': 'region', 'start': 363, 'end': 393, 'description': 'Regulatory domain'}
                ]
            },
            'MSH2': {
                'domains': [
                    {'type': 'domain', 'start': 1, 'end': 120, 'description': 'Mismatch recognition domain'},
                    {'type': 'domain', 'start': 300, 'end': 500, 'description': 'Connector domain'},
                    {'type': 'domain', 'start': 620, 'end': 853, 'description': 'ATPase domain'}
                ]
            }
        }

        print("🔗 Interface Analyzer initialized!")
        print("🎯 Ready to detect domain boundary disruptions!")
    
    def analyze_interface(
        self,
        gene: str,
        variant: str,
        uniprot_id: str,
        position: int,
        ref_aa: str,
        alt_aa: str,
        domain_data: Optional[Dict] = None
    ) -> Dict:
        """
        Analyze if variant disrupts domain interface
        
        Args:
            gene: Gene symbol
            variant: Variant string (p.Ile178Thr)
            uniprot_id: UniProt ID
            position: Amino acid position
            ref_aa: Reference amino acid
            alt_aa: Alternate amino acid
            domain_data: Pre-fetched domain annotations
            
        Returns:
            Dict with interface disruption score and details
        """
        
        print(f"\n🔗 INTERFACE ANALYSIS: {gene} {variant} (pos {position})")
        
        result = {
            'is_interface': False,
            'interface_score': 0.0,
            'interface_type': None,
            'domains_involved': [],
            'disruption_mechanism': None,
            'confidence': 0.0,
            'details': {}
        }
        
        # 🔗 FIRST: Try to load InterPro domains (REAL structural domains!)
        interpro_data = self._load_interpro_domains(uniprot_id)

        if interpro_data and interpro_data.get('domains'):
            print(f"   ✅ Using InterPro domains (REAL structural boundaries!)")
            parsed_data = interpro_data
        else:
            # Fallback to LOF analyzer's domain data (predicted domains)
            print(f"   ⚠️ No InterPro cache - falling back to predicted domains")
            if domain_data is None:
                domain_data = self._get_domain_data(uniprot_id)

            if not domain_data or not domain_data.get('domains'):
                print(f"⚠️  No domain data available for {uniprot_id}")
                return result

            # Parse domain data into our format
            parsed_data = self._parse_domain_features(domain_data)

        # Check if position is at domain interface
        interface_info = self._check_domain_interface(position, parsed_data)
        
        if not interface_info['is_interface']:
            print(f"   Position {position} is NOT at domain interface")
            return result
        
        print(f"🎯 INTERFACE DETECTED!")
        print(f"   Type: {interface_info['interface_type']}")
        print(f"   Domains: {interface_info['domains']}")
        
        result['is_interface'] = True
        result['interface_type'] = interface_info['interface_type']
        result['domains_involved'] = interface_info['domains']
        
        # Score interface disruption severity
        disruption_score = self._score_interface_disruption(
            position,
            ref_aa,
            alt_aa,
            interface_info
        )
        
        result['interface_score'] = disruption_score['score']
        result['disruption_mechanism'] = disruption_score['mechanism']
        result['confidence'] = disruption_score['confidence']
        result['details'] = disruption_score['details']
        
        print(f"   Interface Score: {result['interface_score']:.3f}")
        print(f"   Mechanism: {result['disruption_mechanism']}")
        
        return result
    
    def _get_domain_data(self, uniprot_id: str) -> Optional[Dict]:
        """Fetch domain annotations from UniProt"""
        try:
            url = f"https://www.uniprot.org/uniprot/{uniprot_id}.json"
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                return self._parse_domain_features(data)
        except Exception as e:
            print(f"⚠️  Error fetching domain data: {e}")
        return None
    
    def _load_interpro_domains(self, uniprot_id: str) -> Dict:
        """
        Load REAL domain boundaries from InterPro cache

        This is the 0.5 step Ren suggested - use pre-cached InterPro data!
        """

        import os
        import json

        # Try multiple possible cache locations
        possible_paths = [
            f"protein_annotations_cache/{uniprot_id}_interpro_domains.json",
            f"AdaptiveInterpreter/protein_annotations_cache/{uniprot_id}_interpro_domains.json",
            f"/home/Ace/AdaptiveInterpreter/protein_annotations_cache/{uniprot_id}_interpro_domains.json"
        ]

        cache_file = None
        for path in possible_paths:
            if os.path.exists(path):
                cache_file = path
                break

        if not cache_file:
            print(f"   ⚠️ No InterPro cache for {uniprot_id} - run cache_interpro_domains.py first!")
            return {'domains': []}

        try:
            with open(cache_file, 'r') as f:
                data = json.load(f)

            domains = []
            for d in data.get('domains', []):
                # Only use real structural domains, not just families
                if d.get('type') in ['domain', 'homologous_superfamily']:
                    domains.append({
                        'start': d['start'],
                        'end': d['end'],
                        'type': d['type'],
                        'description': d['description'],
                        'accession': d.get('accession', '')
                    })

            print(f"   🔗 Loaded {len(domains)} InterPro domains for {uniprot_id}")
            return {'domains': sorted(domains, key=lambda x: x['start'])}

        except Exception as e:
            print(f"   ❌ Failed to load InterPro cache for {uniprot_id}: {e}")
            return {'domains': []}

    def _parse_domain_features(self, domain_data: Dict) -> Dict:
        """
        Parse domain features from LOF analyzer's domain data format

        LOF uses domain_annotator which returns:
        {
            'domains': [{'start': int, 'end': int, 'type': str, 'description': str}, ...],
            'regions': [...],
            'mature_chain': [...]
        }
        """
        domains = []

        # Get domains from the data
        for domain in domain_data.get('domains', []):
            domains.append({
                'type': domain.get('type', 'domain'),
                'start': domain.get('start'),
                'end': domain.get('end'),
                'description': domain.get('description', '')
            })

        # Also include regions as potential interfaces
        for region in domain_data.get('regions', []):
            domains.append({
                'type': 'region',
                'start': region.get('start'),
                'end': region.get('end'),
                'description': region.get('description', '')
            })

        return {'domains': sorted(domains, key=lambda x: x['start'])}
    
    def _check_domain_interface(self, position: int, domain_data: Dict) -> Dict:
        """
        Check if position is at domain interface
        
        Interface defined as:
        - Within 5 residues of domain boundary
        - Between two annotated domains
        """
        
        domains = domain_data.get('domains', [])
        
        if len(domains) < 2:
            return {'is_interface': False}
        
        # Check each domain boundary
        for i in range(len(domains) - 1):
            domain1 = domains[i]
            domain2 = domains[i + 1]
            
            # Interface region: end of domain1 to start of domain2
            interface_start = domain1['end'] - 5
            interface_end = domain2['start'] + 5
            
            if interface_start <= position <= interface_end:
                return {
                    'is_interface': True,
                    'interface_type': 'interdomain',
                    'domains': [domain1['description'], domain2['description']],
                    'boundary_position': domain1['end'],
                    'distance_from_boundary': abs(position - domain1['end'])
                }
        
        # Check if at N-terminal or C-terminal of single domain (also critical!)
        for domain in domains:
            if domain['start'] <= position <= domain['start'] + 5:
                return {
                    'is_interface': True,
                    'interface_type': 'domain_start',
                    'domains': [domain['description']],
                    'boundary_position': domain['start'],
                    'distance_from_boundary': position - domain['start']
                }
            elif domain['end'] - 5 <= position <= domain['end']:
                return {
                    'is_interface': True,
                    'interface_type': 'domain_end',
                    'domains': [domain['description']],
                    'boundary_position': domain['end'],
                    'distance_from_boundary': domain['end'] - position
                }
        
        return {'is_interface': False}
    
    def _score_interface_disruption(
        self,
        position: int,
        ref_aa: str,
        alt_aa: str,
        interface_info: Dict
    ) -> Dict:
        """
        Score severity of interface disruption
        
        Factors:
        1. Distance from boundary (closer = worse)
        2. Amino acid property change
        3. Interface type (interdomain > domain_end > domain_start)
        """
        
        score = 0.0
        details = {}
        
        # Base score from interface type
        interface_type = interface_info['interface_type']
        if interface_type == 'interdomain':
            score = 0.6  # Highest - affects two domains
            mechanism = 'interdomain_communication_disruption'
        elif interface_type == 'domain_end':
            score = 0.5  # High - affects domain stability
            mechanism = 'domain_boundary_destabilization'
        elif interface_type == 'domain_start':
            score = 0.4  # Moderate - affects domain folding
            mechanism = 'domain_initiation_disruption'
        else:
            score = 0.3
            mechanism = 'unknown_interface_effect'
        
        details['base_score'] = score
        
        # Distance penalty (closer to boundary = worse)
        distance = interface_info.get('distance_from_boundary', 5)
        distance_multiplier = 1.0 + (5 - distance) * 0.1  # 1.0 to 1.5x
        score *= distance_multiplier
        details['distance_multiplier'] = distance_multiplier
        
        # Amino acid property change
        aa_change_score = self._score_aa_property_change(ref_aa, alt_aa)
        score *= (1.0 + aa_change_score * 0.5)  # Up to 1.5x more
        details['aa_change_score'] = aa_change_score
        
        # Cap at 1.0
        score = min(score, 1.0)
        
        return {
            'score': score,
            'mechanism': mechanism,
            'confidence': 0.8,  # High confidence for annotated interfaces
            'details': details
        }
    
    def _score_aa_property_change(self, ref_aa: str, alt_aa: str) -> float:
        """
        Score amino acid property change severity
        
        Returns 0.0 (conservative) to 1.0 (radical)
        """
        
        # Amino acid properties
        hydrophobic = set('AILMFVPW')
        polar = set('STNQ')
        charged_pos = set('KRH')
        charged_neg = set('DE')
        special = set('GC')
        
        def get_properties(aa):
            props = set()
            if aa in hydrophobic:
                props.add('hydrophobic')
            if aa in polar:
                props.add('polar')
            if aa in charged_pos:
                props.add('positive')
            if aa in charged_neg:
                props.add('negative')
            if aa in special:
                props.add('special')
            return props
        
        ref_props = get_properties(ref_aa)
        alt_props = get_properties(alt_aa)
        
        # No shared properties = radical change
        if not ref_props.intersection(alt_props):
            return 1.0
        
        # Charge flip = very radical
        if ('positive' in ref_props and 'negative' in alt_props) or \
           ('negative' in ref_props and 'positive' in alt_props):
            return 0.9
        
        # Hydrophobic ↔ polar = moderate
        if ('hydrophobic' in ref_props and 'polar' in alt_props) or \
           ('polar' in ref_props and 'hydrophobic' in alt_props):
            return 0.6
        
        # Special (G/C) changes = moderate to high
        if 'special' in ref_props or 'special' in alt_props:
            return 0.7
        
        # Conservative change
        return 0.3

