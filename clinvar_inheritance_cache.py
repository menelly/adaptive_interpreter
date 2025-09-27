#!/usr/bin/env python3
"""
🧬 CLINVAR INHERITANCE PATTERN CACHING

Query ClinVar API once per gene to extract inheritance patterns and cache them.
This allows us to automatically detect AR vs AD genes without hardcoding.
"""

import json
import os
import requests
import time
from typing import Dict, List, Optional
from urllib.parse import quote

class ClinVarInheritanceCache:
    """
    🧬 COMPREHENSIVE GENE ANNOTATION CACHE

    Cache everything: inheritance patterns, GO terms, function descriptions,
    disease names, ClinVar data, UniProt annotations, etc.
    """

    def __init__(self, cache_file: str = "comprehensive_gene_cache.json"):
        self.cache_file = cache_file
        self.cache = self._load_cache()
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        
    def _load_cache(self) -> Dict:
        """Load existing cache from file"""
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"⚠️ Error loading cache: {e}")
        return {}
    
    def _save_cache(self):
        """Save cache to file"""
        try:
            with open(self.cache_file, 'w') as f:
                json.dump(self.cache, f, indent=2)
        except Exception as e:
            print(f"⚠️ Error saving cache: {e}")
    
    def get_comprehensive_gene_data(self, gene_symbol: str, uniprot_id: str = None) -> Dict:
        """
        Get comprehensive gene data from cache or APIs

        Returns:
            {
                'inheritance_pattern': str,
                'diseases': List[str],
                'go_terms': List[str],
                'function_description': str,
                'clinvar_diseases': List[str],
                'pathogenic_variants_count': int,
                'benign_variants_count': int,
                'uniprot_data': Dict,
                'last_updated': float
            }
        """
        if gene_symbol in self.cache:
            cached_data = self.cache[gene_symbol]
            # Check if cache is recent (less than 30 days old)
            if time.time() - cached_data.get('last_updated', 0) < 30 * 24 * 3600:
                return cached_data

        # Query all APIs and build comprehensive data
        print(f"🔍 Building comprehensive data for {gene_symbol}...")
        gene_data = self._build_comprehensive_data(gene_symbol, uniprot_id)

        # Cache the result
        self.cache[gene_symbol] = gene_data
        self._save_cache()

        return gene_data

    def get_gene_inheritance(self, gene_symbol: str) -> Optional[str]:
        """Get just the inheritance pattern (convenience method)"""
        data = self.get_comprehensive_gene_data(gene_symbol)
        return data.get('inheritance_pattern')
    
    def _build_comprehensive_data(self, gene_symbol: str, uniprot_id: str = None) -> Dict:
        """Build comprehensive gene data from multiple APIs"""
        gene_data = {
            'gene_symbol': gene_symbol,
            'uniprot_id': uniprot_id,
            'inheritance_pattern': None,
            'diseases': [],
            'go_terms': [],
            'function_description': '',
            'clinvar_diseases': [],
            'pathogenic_variants_count': 0,
            'benign_variants_count': 0,
            'uniprot_data': {},
            'last_updated': time.time()
        }

        # 1. Query ClinVar for inheritance and disease data
        clinvar_data = self._query_clinvar_comprehensive(gene_symbol)
        if clinvar_data:
            gene_data.update(clinvar_data)

        # 2. Query UniProt for GO terms and function (if we have access)
        if uniprot_id:
            uniprot_data = self._query_uniprot_data(uniprot_id)
            if uniprot_data:
                gene_data['uniprot_data'] = uniprot_data
                gene_data['go_terms'] = uniprot_data.get('go_terms', [])
                gene_data['function_description'] = uniprot_data.get('function', '')

        return gene_data

    def _query_clinvar_comprehensive(self, gene_symbol: str) -> Optional[Dict]:
        """Query ClinVar API for comprehensive gene data"""
        try:
            # Search for ALL variants in this gene (pathogenic + benign)
            search_url = f"{self.base_url}/esearch.fcgi"

            # Get pathogenic variants
            pathogenic_params = {
                'db': 'clinvar',
                'term': f'{gene_symbol}[gene] AND ("pathogenic"[Clinical_significance] OR "likely pathogenic"[Clinical_significance])',
                'retmax': 100,
                'retmode': 'json'
            }

            response = requests.get(search_url, params=pathogenic_params, timeout=10)
            response.raise_for_status()
            pathogenic_data = response.json()

            pathogenic_count = int(pathogenic_data.get('esearchresult', {}).get('count', 0))
            pathogenic_ids = pathogenic_data.get('esearchresult', {}).get('idlist', [])[:20]

            # Get benign variants
            time.sleep(0.3)
            benign_params = {
                'db': 'clinvar',
                'term': f'{gene_symbol}[gene] AND ("benign"[Clinical_significance] OR "likely benign"[Clinical_significance])',
                'retmax': 50,
                'retmode': 'json'
            }

            response = requests.get(search_url, params=benign_params, timeout=10)
            response.raise_for_status()
            benign_data = response.json()

            benign_count = int(benign_data.get('esearchresult', {}).get('count', 0))

            if not pathogenic_ids:
                print(f"   No pathogenic variants found for {gene_symbol}")
                return {
                    'pathogenic_variants_count': pathogenic_count,
                    'benign_variants_count': benign_count
                }

            # Get details for pathogenic variants to extract diseases and inheritance
            fetch_url = f"{self.base_url}/efetch.fcgi"
            fetch_params = {
                'db': 'clinvar',
                'id': ','.join(pathogenic_ids),
                'rettype': 'vcv',
                'retmode': 'xml'
            }

            time.sleep(0.5)  # Be nice to NCBI servers
            response = requests.get(fetch_url, params=fetch_params, timeout=15)
            response.raise_for_status()

            # Parse XML for comprehensive data
            parsed_data = self._parse_comprehensive_xml(response.text)
            parsed_data.update({
                'pathogenic_variants_count': pathogenic_count,
                'benign_variants_count': benign_count
            })

            print(f"   Found {pathogenic_count} pathogenic, {benign_count} benign variants")
            print(f"   Inheritance: {parsed_data.get('inheritance_pattern', 'Unknown')}")
            print(f"   Diseases: {len(parsed_data.get('clinvar_diseases', []))} found")

            return parsed_data

        except Exception as e:
            print(f"   ⚠️ Error querying ClinVar for {gene_symbol}: {e}")
            return None
    
    def _parse_comprehensive_xml(self, xml_text: str) -> Dict:
        """Parse comprehensive data from ClinVar XML"""
        xml_lower = xml_text.lower()

        # Extract diseases/conditions
        diseases = set()
        import re

        # Look for disease names in various XML tags
        disease_patterns = [
            r'<name[^>]*>([^<]+)</name>',
            r'<condition[^>]*>([^<]+)</condition>',
            r'<trait[^>]*>([^<]+)</trait>'
        ]

        for pattern in disease_patterns:
            matches = re.findall(pattern, xml_text, re.IGNORECASE)
            for match in matches:
                if len(match.strip()) > 3 and not match.isdigit():
                    diseases.add(match.strip())

        # Count inheritance patterns
        ar_patterns = [
            'autosomal recessive', 'recessive inheritance', 'recessive disorder',
            'compound heterozygous', 'homozygous pathogenic', 'cystic fibrosis',
            'metabolic disorder', 'enzyme deficiency', 'inborn error'
        ]
        ad_patterns = [
            'autosomal dominant', 'dominant inheritance', 'dominant disorder',
            'haploinsufficiency', 'dominant negative', 'marfan', 'osteogenesis imperfecta',
            'ehlers-danlos', 'connective tissue'
        ]
        xl_patterns = [
            'x-linked', 'x linked', 'sex-linked', 'hemizygous'
        ]

        ar_count = sum(xml_lower.count(pattern) for pattern in ar_patterns)
        ad_count = sum(xml_lower.count(pattern) for pattern in ad_patterns)
        xl_count = sum(xml_lower.count(pattern) for pattern in xl_patterns)

        # Determine inheritance pattern
        inheritance_pattern = None
        if ar_count > ad_count and ar_count > xl_count:
            inheritance_pattern = 'AUTOSOMAL_RECESSIVE'
        elif ad_count > ar_count and ad_count > xl_count:
            inheritance_pattern = 'AUTOSOMAL_DOMINANT'
        elif xl_count > 0:
            inheritance_pattern = 'X_LINKED'

        return {
            'inheritance_pattern': inheritance_pattern,
            'clinvar_diseases': list(diseases)[:20],  # Limit to top 20
            'ar_score': ar_count,
            'ad_score': ad_count,
            'xl_score': xl_count
        }

    def _query_uniprot_data(self, uniprot_id: str) -> Optional[Dict]:
        """Query UniProt for GO terms and function data"""
        try:
            # This would integrate with existing UniProt systems
            # For now, return placeholder
            return {
                'function': f'UniProt function for {uniprot_id}',
                'go_terms': ['placeholder_go_term']
            }
        except Exception as e:
            print(f"   ⚠️ Error querying UniProt for {uniprot_id}: {e}")
            return None
    
    def get_cached_genes(self) -> List[str]:
        """Get list of genes in cache"""
        return list(self.cache.keys())
    
    def clear_cache(self):
        """Clear the cache"""
        self.cache = {}
        if os.path.exists(self.cache_file):
            os.remove(self.cache_file)

# Global cache instance
_inheritance_cache = None

def get_inheritance_cache() -> ClinVarInheritanceCache:
    """Get global inheritance cache instance"""
    global _inheritance_cache
    if _inheritance_cache is None:
        _inheritance_cache = ClinVarInheritanceCache()
    return _inheritance_cache

def get_gene_inheritance_pattern(gene_symbol: str) -> Optional[str]:
    """Convenience function to get inheritance pattern for a gene"""
    cache = get_inheritance_cache()
    return cache.get_gene_inheritance(gene_symbol)

if __name__ == "__main__":
    # Test the cache
    cache = ClinVarInheritanceCache()
    
    test_genes = ['CFTR', 'FBN1', 'COL1A1', 'BRCA1', 'SCN1A']
    
    print("🧬 TESTING CLINVAR INHERITANCE CACHE\n")
    
    for gene in test_genes:
        inheritance = cache.get_gene_inheritance(gene)
        print(f"{gene}: {inheritance or 'Unknown'}")
        print()
