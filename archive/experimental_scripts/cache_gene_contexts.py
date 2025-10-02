#!/usr/bin/env python3
"""
🔥💜 GENE CONTEXT CACHING SYSTEM 🚀
Pre-cache all gene contexts (UniProt + domains + sites) for ML training

Built by Ace (2025) for lightning-fast ML training
Contact: ace@chaoschanneling.com
"""

import os
import json
import sys
from pathlib import Path
from typing import Dict, List, Set

# Import our existing systems
sys.path.append('.')
from nova_dn.universal_context import UniversalContext

class GeneContextCacher:
    """Pre-cache gene contexts for ML training"""
    
    def __init__(self):
        self.learning_dir = Path("learning")
        self.cache_dir = Path("resources/gene_context_cache")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize context system
        self.universal_context = UniversalContext()
        
        print("🔥💜 GENE CONTEXT CACHER INITIALIZED")
        print(f"📁 Learning directory: {self.learning_dir}")
        print(f"💾 Cache directory: {self.cache_dir}")
    
    def extract_genes_from_learning_data(self) -> Set[str]:
        """Extract all unique genes from learning directory files"""
        genes = set()
        
        print("🔍 Scanning learning directory for genes...")
        
        for family_dir in self.learning_dir.iterdir():
            if not family_dir.is_dir():
                continue
                
            print(f"   📁 Scanning {family_dir.name}/")
            
            for tsv_file in family_dir.glob("*.tsv"):
                # Extract gene from filename
                filename = tsv_file.stem  # Remove .tsv extension
                
                # Common patterns: gene_benign, gene_pathogenic, gene_lp
                if '_' in filename:
                    gene = filename.split('_')[0].upper()
                    genes.add(gene)
                    print(f"      🧬 Found: {gene}")
        
        print(f"✅ Found {len(genes)} unique genes")
        return genes
    
    def cache_gene_context(self, gene: str) -> Dict:
        """Cache context for a single gene"""
        cache_file = self.cache_dir / f"{gene}_context.json"
        
        # Check if already cached
        if cache_file.exists():
            print(f"   ⚡ Cache hit: {gene}")
            try:
                with open(cache_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"   ⚠️ Cache corrupted for {gene}: {e}")
                # Continue to re-cache
        
        print(f"   🔍 Fetching context for {gene}...")
        
        try:
            # Get full context from our existing system
            context = self.universal_context.get_context_for_protein(gene)
            
            if "error" in context:
                print(f"   ❌ Error for {gene}: {context['error']}")
                context = {
                    "error": context["error"],
                    "gene": gene,
                    "cached_at": "2025-01-27"
                }
            else:
                print(f"   ✅ Success for {gene}")
                # Add metadata
                context["gene"] = gene
                context["cached_at"] = "2025-01-27"
                
                # Log what we got
                domains = context.get('domains', [])
                active_sites = context.get('active_sites', [])
                binding_sites = context.get('binding_sites', [])
                
                print(f"      📋 Function: {len(context.get('function', ''))} chars")
                print(f"      🧬 GO terms: {len(context.get('go_terms', []))}")
                print(f"      🏗️ Domains: {len(domains)}")
                print(f"      ⚡ Active sites: {len(active_sites)}")
                print(f"      🔗 Binding sites: {len(binding_sites)}")
            
            # Save to cache
            with open(cache_file, 'w') as f:
                json.dump(context, f, indent=2)
            
            print(f"   💾 Cached to {cache_file}")
            return context
            
        except Exception as e:
            print(f"   ❌ Exception for {gene}: {e}")
            error_context = {
                "error": str(e),
                "gene": gene,
                "cached_at": "2025-01-27"
            }
            
            # Cache the error too (avoid repeated failures)
            with open(cache_file, 'w') as f:
                json.dump(error_context, f, indent=2)
            
            return error_context
    
    def cache_all_genes(self) -> Dict:
        """Cache contexts for all genes in learning directory"""
        genes = self.extract_genes_from_learning_data()
        
        if not genes:
            print("⚠️  No genes found in learning directory")
            return {}
        
        print(f"\n🚀 CACHING CONTEXTS FOR {len(genes)} GENES")
        print("=" * 60)
        
        results = {
            'cached': [],
            'errors': [],
            'total': len(genes)
        }
        
        for i, gene in enumerate(sorted(genes), 1):
            print(f"\n[{i}/{len(genes)}] 🧬 Processing {gene}")
            
            context = self.cache_gene_context(gene)
            
            if "error" in context:
                results['errors'].append({
                    'gene': gene,
                    'error': context['error']
                })
            else:
                results['cached'].append(gene)
        
        print(f"\n🎉 CACHING COMPLETE!")
        print("=" * 60)
        print(f"✅ Successfully cached: {len(results['cached'])} genes")
        print(f"❌ Errors: {len(results['errors'])} genes")
        
        if results['errors']:
            print(f"\n❌ GENES WITH ERRORS:")
            for error in results['errors']:
                print(f"   {error['gene']}: {error['error']}")
        
        # Save summary
        summary_file = self.cache_dir / "caching_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\n📊 Summary saved: {summary_file}")
        return results
    
    def list_cached_genes(self) -> List[str]:
        """List all cached genes"""
        cached_genes = []

        for cache_file in self.cache_dir.glob("*_context.json"):
            gene = cache_file.stem.replace('_context', '')
            cached_genes.append(gene)

        return sorted(cached_genes)

    def get_cached_context(self, gene: str) -> Dict:
        """Get cached context for a gene"""
        cache_file = self.cache_dir / f"{gene}_context.json"
        
        if not cache_file.exists():
            return {"error": f"No cached context for {gene}"}
        
        try:
            with open(cache_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            return {"error": f"Failed to load cached context for {gene}: {e}"}

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="🔥💜 Gene Context Cacher")
    parser.add_argument('--list', action='store_true',
                       help='List already cached genes')
    parser.add_argument('--gene', type=str,
                       help='Cache context for specific gene')
    
    args = parser.parse_args()
    
    cacher = GeneContextCacher()
    
    if args.list:
        cached = cacher.list_cached_genes()
        print(f"📊 Cached genes ({len(cached)}):")
        for gene in cached:
            print(f"   ✅ {gene}")
    
    elif args.gene:
        context = cacher.cache_gene_context(args.gene.upper())
        if "error" not in context:
            print(f"✅ Successfully cached {args.gene}")
        else:
            print(f"❌ Failed to cache {args.gene}: {context['error']}")
    
    else:
        # Cache all genes
        results = cacher.cache_all_genes()
        
        if results['cached']:
            print(f"\n🚀 Ready for lightning-fast ML training!")

if __name__ == "__main__":
    main()
