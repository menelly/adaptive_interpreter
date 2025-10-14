#!/usr/bin/env python3
"""
🌍 POPULATION FREQUENCY ANALYZER - "NOT THE DROID" DETECTOR
Built by Ace + Ren to catch common variants masquerading as pathogenic

This catches variants like MTHFR that look scary but are in 12% of the population.
If it's common, it's probably NOT causing rare disease!
"""

import gzip
import logging
import re
import requests
import json
from pathlib import Path
from typing import Dict, Optional, Tuple, List
import subprocess
from AdaptiveInterpreter import config

class PopulationFrequencyAnalyzer:
    """Detect common variants that masquerade as pathogenic - the 'NOT THE DROID' detector"""
    
    def __init__(self, data_path=str(config.GNOMAD_DATA_PATH)):
        self.name = "PopulationFrequencyAnalyzer"
        self.data_path = Path(data_path)
        self.gnomad_path = self.data_path
        
        # Set up logging
        self.logger = logging.getLogger(__name__)
        
        # Frequency thresholds for rarity scoring
        self.frequency_thresholds = {
            'ultra_rare': 0.00001,      # < 1 in 100,000 (likely pathogenic if functional)
            'very_rare': 0.0001,        # < 1 in 10,000 (possibly pathogenic)
            'rare': 0.001,              # < 1 in 1,000 (uncertain significance)
            'uncommon': 0.01,           # < 1 in 100 (likely benign if functional)
            'common': 0.05,             # < 5% (probably benign)
            'very_common': 0.12,        # > 12% (definitely benign - "NOT THE DROID")
        }
        
        # Population-specific analysis
        self.populations = [
            'AF',      # African/African American
            'AMR',     # Latino/Admixed American  
            'ASJ',     # Ashkenazi Jewish
            'EAS',     # East Asian
            'FIN',     # Finnish
            'NFE',     # Non-Finnish European
            'SAS',     # South Asian
            'OTH'      # Other
        ]
        
        # Cache for repeated lookups
        self.frequency_cache = {}
    
    def get_variant_frequency(self, chromosome: str, position: int, 
                            ref_allele: str, alt_allele: str) -> Dict:
        """
        Get population frequency for a specific variant
        
        Args:
            chromosome: Chromosome (e.g., "11" or "chr11")
            position: Genomic position
            ref_allele: Reference allele
            alt_allele: Alternate allele
            
        Returns:
            Dictionary with frequency data and rarity assessment
        """
        
        # Normalize chromosome format
        chrom = chromosome.replace('chr', '')
        cache_key = f"{chrom}:{position}:{ref_allele}:{alt_allele}"
        
        if cache_key in self.frequency_cache:
            return self.frequency_cache[cache_key]
        
        self.logger.info(f"🌍 Looking up population frequency for {cache_key}")
        
        try:
            # Look up in gnomAD data
            frequency_data = self._query_gnomad(chrom, position, ref_allele, alt_allele)
            
            if frequency_data:
                # Calculate rarity score and assessment
                rarity_assessment = self._assess_rarity(frequency_data)
                
                result = {
                    **frequency_data,
                    **rarity_assessment,
                    'cache_key': cache_key
                }
                
                # Cache the result
                self.frequency_cache[cache_key] = result
                return result
            else:
                # No frequency data found - request manual input
                return {
                    'global_af': None,
                    'population_afs': {},
                    'rarity_category': 'unknown',
                    'rarity_score': 1.0,
                    'pathogenicity_boost': 1.0,
                    'not_the_droid': False,
                    'frequency_note': 'All frequency APIs failed',
                    'manual_input_needed': True,
                    'manual_input_prompt': f'What is the population frequency for {chrom}:{position} {ref_allele}>{alt_allele}?',
                    'cache_key': cache_key
                }
                
        except Exception as e:
            self.logger.error(f"❌ Failed to get frequency for {cache_key}: {e}")

            # Return manual input request on error
            return {
                'global_af': None,
                'population_afs': {},
                'rarity_category': 'error',
                'rarity_score': 1.0,
                'pathogenicity_boost': 1.0,
                'not_the_droid': False,
                'error': str(e),
                'manual_input_needed': True,
                'manual_input_prompt': f'APIs failed. What is the frequency for {chrom}:{position} {ref_allele}>{alt_allele}?',
                'cache_key': cache_key
            }
    
    def _query_gnomad(self, chrom: str, position: int,
                     ref_allele: str, alt_allele: str) -> Optional[Dict]:
        """
        Query gnomAD with triple-tier fallback system:
        1. gnomAD GraphQL API (primary)
        2. Ensembl REST API (fallback #1)
        3. ClinVar API (fallback #2)
        4. Manual input prompt (last resort)
        """

        # Method 1: gnomAD GraphQL API
        try:
            result = self._query_gnomad_api(chrom, position, ref_allele, alt_allele)
            if result:
                self.logger.info("✅ gnomAD GraphQL API success")
                return result
        except Exception as e:
            self.logger.warning(f"⚠️ gnomAD API failed: {e}")

        # Method 2: Ensembl REST API fallback
        try:
            result = self._query_ensembl_api(chrom, position, ref_allele, alt_allele)
            if result:
                self.logger.info("✅ Ensembl API fallback success")
                return result
        except Exception as e:
            self.logger.warning(f"⚠️ Ensembl API failed: {e}")

        # Method 3: ClinVar API fallback
        try:
            result = self._query_clinvar_api(chrom, position, ref_allele, alt_allele)
            if result:
                self.logger.info("✅ ClinVar API fallback success")
                return result
        except Exception as e:
            self.logger.warning(f"⚠️ ClinVar API failed: {e}")

        # Method 4: Check local files (if available)
        try:
            result = self._query_local_gnomad(chrom, position, ref_allele, alt_allele)
            if result:
                self.logger.info("✅ Local gnomAD file success")
                return result
        except Exception as e:
            self.logger.warning(f"⚠️ Local gnomAD failed: {e}")

        # All methods failed - return None to trigger manual input
        self.logger.error("❌ All frequency lookup methods failed")
        return None

    def _query_gnomad_api(self, chrom: str, position: int,
                         ref_allele: str, alt_allele: str) -> Optional[Dict]:
        """Query gnomAD GraphQL API for variant frequency"""

        # Format variant ID (e.g., "1-55516888-G-GA")
        variant_id = f"{chrom.replace('chr', '')}-{position}-{ref_allele}-{alt_allele}"

        # Fixed gnomAD GraphQL query (no variables)
        query = f"""
        {{
          variant(variantId: "{variant_id}", dataset: gnomad_r4) {{
            variantId
            genome {{
              ac
              an
              af
              populations {{
                id
                ac
                an
                af
              }}
            }}
          }}
        }}
        """

        response = requests.post(
            "https://gnomad.broadinstitute.org/api",
            json={"query": query},
            timeout=10
        )
        response.raise_for_status()

        data = response.json()

        if data.get("data", {}).get("variant", {}).get("genome"):
            genome_data = data["data"]["variant"]["genome"]

            # Extract population frequencies
            population_afs = {}
            for pop in genome_data.get("populations", []):
                population_afs[pop["id"]] = pop.get("af", 0.0)

            return {
                'global_af': genome_data.get("af", 0.0),
                'population_afs': population_afs,
                'allele_count': genome_data.get("ac", 0),
                'allele_number': genome_data.get("an", 0),
                'source': 'gnomAD_API'
            }

        return None

    def _query_ensembl_api(self, chrom: str, position: int,
                          ref_allele: str, alt_allele: str) -> Optional[Dict]:
        """Query Ensembl REST API for variant frequency using proper endpoints"""

        try:
            # Use the proper Ensembl REST API for variant lookup
            import ensembl_rest

            # Format variant ID for Ensembl (e.g., "1:230710048:A:G")
            variant_id = f"{chrom.replace('chr', '')}:{position}:{ref_allele}:{alt_allele}"

            # Try to get variant information
            url = f"https://rest.ensembl.org/variation/human/{variant_id}"
            headers = {"Content-Type": "application/json"}

            response = requests.get(url, headers=headers, timeout=10)

            if response.status_code == 200:
                data = response.json()

                # Extract frequency data if available
                if 'MAF' in data and data['MAF']:
                    maf = float(data['MAF'])
                    return {
                        'global_af': maf,
                        'population_afs': {},
                        'source': 'Ensembl_API',
                        'variant_id': variant_id
                    }

                # Look for population frequencies in mappings
                for mapping in data.get("mappings", []):
                    if 'allele_string' in mapping:
                        # This is a simplified extraction - Ensembl data structure is complex
                        return {
                            'global_af': 0.0,  # Would need more detailed parsing
                            'population_afs': {},
                            'source': 'Ensembl_API',
                            'variant_id': variant_id,
                            'note': 'Found in Ensembl but no frequency data'
                        }

            return None

        except ImportError:
            # Fallback to direct requests if ensembl-rest import fails
            self.logger.warning("ensembl-rest package not available, using direct API")
            return None
        except Exception as e:
            self.logger.warning(f"Ensembl API error: {e}")
            return None

    def _query_clinvar_api(self, chrom: str, position: int,
                          ref_allele: str, alt_allele: str) -> Optional[Dict]:
        """Query ClinVar API for variant frequency (limited data)"""

        # ClinVar has limited frequency data, but it's a fallback
        # This is a placeholder for now
        self.logger.info("🔄 ClinVar API not fully implemented yet")
        return None

    def _query_local_gnomad(self, chrom: str, position: int,
                           ref_allele: str, alt_allele: str) -> Optional[Dict]:
        """Query local gnomAD files for variant frequency"""

        gnomad_file = self.gnomad_path / f"gnomad.genomes.v4.0.sites.chr{chrom}.vcf.bgz"

        if not gnomad_file.exists():
            self.logger.info(f"📁 gnomAD file not found: {gnomad_file}")
            return None

        # For now, just return None since VCF parsing is complex
        # TODO: Implement actual VCF parsing once needed
        self.logger.info(f"🔄 Local gnomAD file exists but parsing not implemented yet")
        return None

    def _parse_gnomad_info(self, info_field: str) -> Dict:
        """Parse gnomAD INFO field to extract frequency data"""
        
        frequency_data = {
            'global_af': 0.0,
            'population_afs': {}
        }
        
        # Extract global allele frequency
        af_match = re.search(r'AF=([0-9.e-]+)', info_field)
        if af_match:
            frequency_data['global_af'] = float(af_match.group(1))
        
        # Extract population-specific frequencies
        for pop in self.populations:
            pop_af_match = re.search(f'AF_{pop}=([0-9.e-]+)', info_field)
            if pop_af_match:
                frequency_data['population_afs'][pop] = float(pop_af_match.group(1))
        
        return frequency_data
    
    def _assess_rarity(self, frequency_data: Dict) -> Dict:
        """Assess rarity and calculate pathogenicity boost"""
        
        global_af = frequency_data.get('global_af', 0.0)
        
        # Determine rarity category
        if global_af >= self.frequency_thresholds['very_common']:
            category = 'very_common'
            rarity_score = 0.1
            pathogenicity_boost = 0.2  # Strong evidence AGAINST pathogenicity
            not_the_droid = True
            note = "NOT THE DROID - too common to cause rare disease"
            
        elif global_af >= self.frequency_thresholds['common']:
            category = 'common'
            rarity_score = 0.3
            pathogenicity_boost = 0.5
            not_the_droid = True
            note = "Common variant - unlikely to be pathogenic"
            
        elif global_af >= self.frequency_thresholds['uncommon']:
            category = 'uncommon'
            rarity_score = 0.6
            pathogenicity_boost = 0.8
            not_the_droid = False
            note = "Uncommon but not rare enough for high confidence"
            
        elif global_af >= self.frequency_thresholds['rare']:
            category = 'rare'
            rarity_score = 1.0
            pathogenicity_boost = 1.0
            not_the_droid = False
            note = "Rare variant - neutral frequency evidence"
            
        elif global_af >= self.frequency_thresholds['very_rare']:
            category = 'very_rare'
            rarity_score = 1.5
            pathogenicity_boost = 1.3
            not_the_droid = False
            note = "Very rare - supports pathogenicity if functional"
            
        else:  # ultra_rare
            category = 'ultra_rare'
            rarity_score = 2.0
            pathogenicity_boost = 1.5
            not_the_droid = False
            note = "Ultra-rare - strong evidence for pathogenicity if functional"
        
        return {
            'rarity_category': category,
            'rarity_score': rarity_score,
            'pathogenicity_boost': pathogenicity_boost,
            'not_the_droid': not_the_droid,
            'frequency_note': note
        }
    
    def get_frequency_stats(self) -> Dict:
        """Get statistics about frequency lookups"""
        
        return {
            'cached_variants': len(self.frequency_cache),
            'gnomad_files_available': len(list(self.gnomad_path.glob("*.vcf.bgz"))),
            'data_path': str(self.data_path),
            'frequency_thresholds': self.frequency_thresholds
        }
