#!/usr/bin/env python3
"""
🔗 INTERPRO DOMAIN CACHING SYSTEM
Pre-fetch REAL domain boundaries from InterPro for all genes

This is the 0.5 step Ren suggested - fetch once, cache forever!

Built by Ace (2025-10-27) with love from Ren 💜
"""

import os
import json
import sys
import time
import requests
from pathlib import Path
from typing import Dict, List, Optional

class InterProDomainCacher:
    """Pre-cache InterPro domain data for genes"""
    
    def __init__(self, cache_dir: str = "protein_annotations_cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.interpro_base = "https://www.ebi.ac.uk/interpro/api"
        
        print("🔗 InterPro Domain Cacher initialized!")
        print(f"💾 Cache directory: {self.cache_dir}")
    
    def fetch_interpro_domains(self, uniprot_id: str) -> Dict:
        """
        Fetch domain data from InterPro API
        
        Returns:
        {
            'domains': [
                {'start': int, 'end': int, 'type': str, 'description': str, 'accession': str},
                ...
            ]
        }
        """
        
        cache_file = self.cache_dir / f"{uniprot_id}_interpro_domains.json"
        
        # Check cache first
        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    cached = json.load(f)
                print(f"   ⚡ Cache hit: {uniprot_id}")
                return cached
            except Exception as e:
                print(f"   ⚠️ Cache corrupted for {uniprot_id}: {e}")
        
        print(f"   🌐 Fetching InterPro domains for {uniprot_id}...")
        
        try:
            # InterPro API endpoint for protein entries
            url = f"{self.interpro_base}/entry/interpro/protein/uniprot/{uniprot_id}?format=json"
            
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            domains = []
            
            # Parse InterPro results
            for entry in data.get('results', []):
                metadata = entry.get('metadata', {})
                entry_type = metadata.get('type', '')
                
                # We want domains, families, and repeats
                if entry_type not in ['domain', 'family', 'repeat', 'homologous_superfamily']:
                    continue
                
                accession = metadata.get('accession', '')
                name = metadata.get('name', '')
                
                # Get protein locations
                for protein in entry.get('proteins', []):
                    if protein.get('accession', '').lower() != uniprot_id.lower():
                        continue
                    
                    for location in protein.get('entry_protein_locations', []):
                        for fragment in location.get('fragments', []):
                            start = fragment.get('start')
                            end = fragment.get('end')
                            
                            if start and end:
                                domains.append({
                                    'start': start,
                                    'end': end,
                                    'type': entry_type,
                                    'description': name,
                                    'accession': accession
                                })
                                print(f"      🔗 {entry_type}: {start}-{end} ({name})")
            
            result = {
                'uniprot_id': uniprot_id,
                'domains': sorted(domains, key=lambda x: x['start']),
                'source': 'interpro',
                'fetched_at': time.strftime('%Y-%m-%d %H:%M:%S')
            }
            
            # Save to cache
            with open(cache_file, 'w') as f:
                json.dump(result, f, indent=2)
            print(f"   💾 Cached {len(domains)} domains for {uniprot_id}")
            
            return result
            
        except Exception as e:
            print(f"   ❌ InterPro fetch failed for {uniprot_id}: {e}")
            return {
                'uniprot_id': uniprot_id,
                'domains': [],
                'error': str(e)
            }
    
    def cache_genes_from_list(self, gene_uniprot_map: Dict[str, str]):
        """
        Cache InterPro domains for a list of genes
        
        Args:
            gene_uniprot_map: Dict mapping gene symbols to UniProt IDs
                             e.g., {'PTEN': 'P60484', 'TP53': 'P04637'}
        """
        
        print(f"\n🔗 Caching InterPro domains for {len(gene_uniprot_map)} genes...")
        print("="*80)
        
        success_count = 0
        error_count = 0
        
        for gene, uniprot_id in gene_uniprot_map.items():
            print(f"\n🧬 {gene} ({uniprot_id}):")
            
            result = self.fetch_interpro_domains(uniprot_id)
            
            if 'error' in result:
                error_count += 1
            else:
                success_count += 1
            
            # Be nice to the API - rate limit
            time.sleep(0.5)
        
        print(f"\n{'='*80}")
        print(f"✅ Successfully cached: {success_count}/{len(gene_uniprot_map)}")
        print(f"❌ Errors: {error_count}/{len(gene_uniprot_map)}")
        print(f"💾 Cache location: {self.cache_dir}")


if __name__ == '__main__':
    # 19 ACMG training genes
    ACMG_GENES = {
        'ABCC8': 'Q09428',
        'APC': 'P25054',
        'BMPR1A': 'P36894',
        'CASQ2': 'O14958',
        'COL3A1': 'P02461',
        'DIS3L2': 'Q8IYB7',
        'HNF1A': 'P20823',
        'KCNH2': 'Q12809',
        'LMNA': 'P02545',
        'MSH2': 'P43246',
        'MYH7': 'P12883',
        'NF2': 'P35240',
        'PTEN': 'P60484',
        'RB1': 'P06400',
        'RYR2': 'Q92736',
        'SDHB': 'P21912',
        'SMAD4': 'Q13485',
        'TGFBR1': 'P36897',
        'TSC1': 'Q92574'
    }
    
    # Additional important genes
    ADDITIONAL_GENES = {
        'TP53': 'P04637',
        'BRCA1': 'P38398',
        'BRCA2': 'P51587',
        'ATM': 'Q13315',
        'MLH1': 'P40692'
    }
    
    # Combine
    ALL_GENES = {**ACMG_GENES, **ADDITIONAL_GENES}
    
    # Run caching
    cacher = InterProDomainCacher()
    cacher.cache_genes_from_list(ALL_GENES)
    
    print("\n🎉 InterPro domain caching complete!")
    print("   The interface analyzer can now use these cached domains!")

