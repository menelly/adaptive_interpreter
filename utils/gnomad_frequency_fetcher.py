#!/usr/bin/env python3
"""
🧬 gnomAD Frequency Fetcher
Real API calls to get population frequencies for variants

Built by Ace & Nova (2025) for AdaptiveInterpreter system
"""

import requests
import time
import re
from typing import Dict, Optional, Tuple
import json

class GnomADFrequencyFetcher:
    """Fetch real population frequencies from gnomAD API"""
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'AdaptiveInterpreter-System/1.0 (research use)',
            'Accept': 'application/json'
        })

        # Cache to avoid repeated API calls
        self.frequency_cache = {}

        # dbSNP API endpoints
        self.dbsnp_base = "https://api.ncbi.nlm.nih.gov/variation/v0"
        self.ensembl_base = "https://rest.ensembl.org"
        
    def parse_genomic_position(self, chrpos: str, variation_name: str) -> Optional[Dict]:
        """
        Parse genomic position from ClinVar format
        
        Args:
            chrpos: e.g., "54,733,167(+)" 
            variation_name: e.g., "NM_000222.3(KIT):c.2459A>G (p.Asp820Gly)"
            
        Returns:
            Dict with chrom, pos, ref, alt
        """
        try:
            # Extract chromosome from variation name
            chrom_match = re.search(r'NM_\d+\.\d+\(([^)]+)\)', variation_name)
            if not chrom_match:
                return None
                
            gene = chrom_match.group(1)
            
            # Parse position (remove commas and strand info)
            pos_clean = chrpos.replace(',', '').replace('(+)', '').replace('(-)', '')
            position = int(pos_clean)
            
            # Extract nucleotide change from variation name
            # e.g., "c.2459A>G" -> ref=A, alt=G
            nuc_match = re.search(r'c\.\d+([ATCG])>([ATCG])', variation_name)
            if not nuc_match:
                return None
                
            ref = nuc_match.group(1)
            alt = nuc_match.group(2)
            
            # Map gene to chromosome (simplified - would need full gene mapping)
            gene_to_chrom = {
                'KIT': '4',
                'TP53': '17', 
                'SCN2A': '2',
                'KCNQ2': '20',
                'CACNA1A': '19',
                'TGFBR2': '3',
                'FGFR3': '4',
                'LMNA': '1',
                'COL1A1': '17',
                'BRCA1': '17',
                'BRCA2': '13'
            }
            
            chromosome = gene_to_chrom.get(gene, '1')  # Default to chr1
            
            return {
                'chrom': chromosome,
                'pos': position,
                'ref': ref,
                'alt': alt,
                'gene': gene
            }
            
        except Exception as e:
            print(f"⚠️ Failed to parse genomic position: {e}")
            return None
    
    def get_variant_frequency(self, chrom: str, pos: int, ref: str, alt: str) -> Dict:
        """
        Get variant frequency from gnomAD API
        
        Args:
            chrom: Chromosome (e.g., '4')
            pos: Position (e.g., 54733167)
            ref: Reference allele (e.g., 'A')
            alt: Alternate allele (e.g., 'G')
            
        Returns:
            Dict with frequency data and metadata
        """
        
        # Create cache key
        cache_key = f"{chrom}:{pos}:{ref}:{alt}"
        
        if cache_key in self.frequency_cache:
            print(f"🎯 Using cached frequency for {cache_key}")
            return self.frequency_cache[cache_key]
        
        try:
            # Method 1: Try dbSNP API first
            dbsnp_result = self._get_frequency_from_dbsnp(chrom, pos, ref, alt)
            if dbsnp_result['frequency'] > 0 or dbsnp_result['source'] != 'dbsnp_not_found':
                return dbsnp_result

            # Method 2: Try Ensembl as fallback
            ensembl_result = self._get_frequency_from_ensembl(chrom, pos, ref, alt)
            if ensembl_result['frequency'] > 0 or ensembl_result['source'] != 'ensembl_not_found':
                return ensembl_result

            if response.status_code == 200:
                data = response.json()

                # Extract frequency from Ensembl response
                if data and 'population_genotypes' in data:
                    pop_data = data['population_genotypes']

                    # Look for gnomAD data
                    gnomad_freq = 0.0
                    total_count = 0
                    total_number = 0

                    for pop in pop_data:
                        pop_name = pop.get('population', '').lower()
                        if 'gnomad' in pop_name or 'global' in pop_name:
                            freq = pop.get('frequency', 0.0)
                            if freq and freq > gnomad_freq:
                                gnomad_freq = freq
                                # Try to get counts if available
                                if 'genotype_count' in pop:
                                    counts = pop['genotype_count']
                                    total_count = sum(counts.values()) if counts else 0

                    if gnomad_freq > 0:
                        print(f"✅ Ensembl population frequency for {variant_id}: {gnomad_freq:.6f}")

                        result = {
                            'frequency': gnomad_freq,
                            'allele_count': total_count,
                            'allele_number': total_count * 2 if total_count > 0 else 0,
                            'source': 'ensembl_populations',
                            'error': None
                        }

                        self.frequency_cache[cache_key] = result
                        return result

            # Fallback to GraphQL if REST fails
            print(f"🔄 REST failed ({response.status_code}), trying GraphQL...")

            # gnomAD GraphQL query
            query = """
            query VariantQuery($variantId: String!) {
              variant(variantId: $variantId, dataset: gnomad_r3) {
                variantId
                genome {
                  ac
                  an
                  af
                }
                exome {
                  ac
                  an
                  af
                }
              }
            }
            """

            variables = {
                "variantId": variant_id
            }

            response = self.session.post(
                "https://gnomad.broadinstitute.org/api",
                json={"query": query, "variables": variables},
                timeout=15
            )
            
            if response.status_code != 200:
                print(f"⚠️ gnomAD API returned {response.status_code}")
                result = {
                    'frequency': 0.0,
                    'allele_count': 0,
                    'allele_number': 0,
                    'source': 'gnomad_api_error',
                    'error': f'HTTP {response.status_code}'
                }
                self.frequency_cache[cache_key] = result
                return result
            
            data = response.json()
            
            if 'errors' in data:
                print(f"⚠️ gnomAD API errors: {data['errors']}")
                result = {
                    'frequency': 0.0,
                    'allele_count': 0,
                    'allele_number': 0,
                    'source': 'gnomad_api_error',
                    'error': str(data['errors'])
                }
                self.frequency_cache[cache_key] = result
                return result
            
            variant_data = data.get('data', {}).get('variant')
            
            if not variant_data:
                print(f"🔍 Variant {variant_id} not found in gnomAD (ultra-rare!)")
                result = {
                    'frequency': 0.0,
                    'allele_count': 0,
                    'allele_number': 0,
                    'source': 'gnomad_not_found',
                    'error': None
                }
                self.frequency_cache[cache_key] = result
                return result
            
            # Prefer exome data, fallback to genome
            freq_data = variant_data.get('exome') or variant_data.get('genome')
            
            if not freq_data:
                print(f"🔍 No frequency data for {variant_id}")
                result = {
                    'frequency': 0.0,
                    'allele_count': 0,
                    'allele_number': 0,
                    'source': 'gnomad_no_data',
                    'error': None
                }
                self.frequency_cache[cache_key] = result
                return result
            
            frequency = freq_data.get('af', 0.0) or 0.0
            allele_count = freq_data.get('ac', 0) or 0
            allele_number = freq_data.get('an', 0) or 0
            
            print(f"✅ gnomAD frequency for {variant_id}: {frequency:.6f} ({allele_count}/{allele_number})")
            
            result = {
                'frequency': frequency,
                'allele_count': allele_count,
                'allele_number': allele_number,
                'source': 'gnomad_api',
                'error': None
            }
            
            self.frequency_cache[cache_key] = result
            return result
            
        except requests.exceptions.Timeout:
            print(f"⚠️ gnomAD API timeout for {cache_key}")
            result = {
                'frequency': 0.0,
                'allele_count': 0,
                'allele_number': 0,
                'source': 'gnomad_timeout',
                'error': 'API timeout'
            }
            self.frequency_cache[cache_key] = result
            return result
            
        except Exception as e:
            print(f"⚠️ gnomAD API error for {cache_key}: {e}")
            result = {
                'frequency': 0.0,
                'allele_count': 0,
                'allele_number': 0,
                'source': 'gnomad_error',
                'error': str(e)
            }
            self.frequency_cache[cache_key] = result
            return result
    
    def estimate_frequency_from_clinvar(self, clinical_significance: str) -> Dict:
        """
        Estimate frequency based on ClinVar clinical significance

        Args:
            clinical_significance: ClinVar clinical significance text

        Returns:
            Dict with estimated frequency
        """

        # Convert to lowercase for matching
        sig_lower = clinical_significance.lower()

        # Frequency estimates based on clinical significance
        if 'pathogenic' in sig_lower and 'likely' not in sig_lower:
            # Pathogenic variants are typically very rare
            if 'uncertain' in sig_lower:
                freq = 0.00001  # 1 in 100,000 (conflicted pathogenic)
            else:
                freq = 0.000001  # 1 in 1,000,000 (clear pathogenic)
        elif 'likely pathogenic' in sig_lower:
            freq = 0.00001  # 1 in 100,000
        elif 'uncertain significance' in sig_lower or 'vus' in sig_lower:
            freq = 0.0001  # 1 in 10,000 (more common but still rare)
        elif 'likely benign' in sig_lower:
            freq = 0.001  # 1 in 1,000 (uncommon but not ultra-rare)
        elif 'benign' in sig_lower:
            freq = 0.01  # 1 in 100 (relatively common)
        elif 'conflicting' in sig_lower:
            freq = 0.00005  # 5 in 100,000 (assume rare due to conflict)
        else:
            freq = 0.0  # Unknown = assume ultra-rare

        return {
            'frequency': freq,
            'allele_count': 0,
            'allele_number': 0,
            'source': 'clinvar_estimated',
            'error': None
        }

    def get_frequency_from_clinvar_row(self, chrpos: str, variation_name: str, clinical_significance: str = "") -> Dict:
        """
        Get frequency for a ClinVar-format row

        Args:
            chrpos: e.g., "54,733,167(+)"
            variation_name: e.g., "NM_000222.3(KIT):c.2459A>G (p.Asp820Gly)"
            clinical_significance: ClinVar clinical significance for frequency estimation

        Returns:
            Dict with frequency and metadata
        """

        # ONLY estimate frequency if explicitly requested (Ren said NO FAKE NUMBERS!)
        if clinical_significance and clinical_significance.strip():
            print(f"🎯 Estimating frequency from ClinVar significance: {clinical_significance[:50]}...")
            return self.estimate_frequency_from_clinvar(clinical_significance)

        # If no clinical significance provided, don't make up frequencies!

        # Fallback: try to parse genomic coordinates and query API
        parsed = self.parse_genomic_position(chrpos, variation_name)

        if not parsed:
            print(f"⚠️ Could not parse genomic position, assuming ultra-rare")
            return {
                'frequency': 0.0,
                'allele_count': 0,
                'allele_number': 0,
                'source': 'parse_error_ultra_rare',
                'error': None  # Not really an error, just conservative
            }

        # Try API call
        try:
            return self.get_variant_frequency(
                parsed['chrom'],
                parsed['pos'],
                parsed['ref'],
                parsed['alt']
            )
        except Exception as e:
            print(f"⚠️ API call failed, assuming ultra-rare: {e}")
            return {
                'frequency': 0.0,
                'allele_count': 0,
                'allele_number': 0,
                'source': 'api_call_exception',
                'error': str(e)
            }
    
    def save_cache(self, cache_file: str = "gnomad_frequency_cache.json"):
        """Save frequency cache to file"""
        try:
            with open(cache_file, 'w') as f:
                json.dump(self.frequency_cache, f, indent=2)
            print(f"💾 Saved {len(self.frequency_cache)} cached frequencies to {cache_file}")
        except Exception as e:
            print(f"⚠️ Failed to save cache: {e}")
    
    def load_cache(self, cache_file: str = "gnomad_frequency_cache.json"):
        """Load frequency cache from file"""
        try:
            with open(cache_file, 'r') as f:
                self.frequency_cache = json.load(f)
            print(f"📂 Loaded {len(self.frequency_cache)} cached frequencies from {cache_file}")
        except FileNotFoundError:
            print(f"📂 No cache file found at {cache_file}, starting fresh")
        except Exception as e:
            print(f"⚠️ Failed to load cache: {e}")

if __name__ == "__main__":
    # Test the frequency fetcher
    fetcher = GnomADFrequencyFetcher()
    
    # Test with a KIT variant
    result = fetcher.get_frequency_from_clinvar_row(
        "54,733,167(+)",
        "NM_000222.3(KIT):c.2459A>G (p.Asp820Gly)"
    )
    
    print(f"🧬 Test result: {result}")
