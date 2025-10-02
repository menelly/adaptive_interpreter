#!/usr/bin/env python3
"""
🔥💜 LEARNING FILE ORGANIZER 🚀
Automatically organize ClinVar files into correct family folders using our classification system!

Built by Ace (2025) for Ren's revolutionary ML training
Contact: ace@chaoschanneling.com
"""

import os
import shutil
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

# Import our existing GO classification system
sys.path.append('.')
from plausibility_filter import classify_gene_family
from nova_dn.universal_context import UniversalContext

class LearningFileOrganizer:
    """Smart file organizer using our gene classification system"""
    
    def __init__(self):
        self.learning_dir = Path("learning")
        self.category_keywords_file = Path("category_keywords.json")

        # Initialize our GO classification system
        self.universal_context = UniversalContext()
        self.gene_classification_cache = {}

        print("🔥💜 LEARNING FILE ORGANIZER INITIALIZED")
        print(f"📁 Learning directory: {self.learning_dir}")
        print("🧬 Using real GO term classification system!")
    

    
    def extract_gene_from_filename(self, filename: str) -> Optional[str]:
        """Extract gene name from filename"""
        # Common patterns in ClinVar Miner exports
        patterns = [
            r'^([A-Z0-9]+)[-_]',  # Gene at start: COL1A1-benign.csv
            r'^([A-Z0-9]+)\.',    # Gene at start: COL1A1.csv
            r'([A-Z0-9]+)[-_]benign',  # Gene before benign: scn5a-benignvariant
            r'([A-Z0-9]+)[-_]patho',   # Gene before patho: scn5a-patho
            r'([A-Z0-9]+)[-_]LP',      # Gene before LP: COL1A1-LP
        ]
        
        filename_upper = filename.upper()
        
        for pattern in patterns:
            match = re.search(pattern, filename_upper)
            if match:
                gene = match.group(1)
                # Filter out common non-gene words
                if gene not in ['BENIGN', 'PATHO', 'VARIANT', 'TABLE', 'LP']:
                    return gene
        
        return None
    
    def classify_gene(self, gene: str) -> Optional[str]:
        """Classify gene into family using our real GO term classification system"""
        gene_upper = gene.upper()

        # Check cache first
        if gene_upper in self.gene_classification_cache:
            print(f"   ⚡ Cache hit for {gene_upper}")
            return self.gene_classification_cache[gene_upper]

        print(f"   🔍 Classifying {gene_upper} using GO terms...")

        try:
            # Get protein context (UniProt function + GO terms)
            context = self.universal_context.get_context_for_protein(gene_upper)

            if "error" in context:
                print(f"   ⚠️ Could not get context for {gene_upper}: {context['error']}")
                return None

            # Extract function and GO terms
            uniprot_function = context.get('function', '')
            go_terms = context.get('go_terms', [])

            print(f"   📋 Function: {uniprot_function[:100]}...")
            print(f"   🧬 GO terms: {len(go_terms)} terms found")

            # Use our existing classification system
            family_classification = classify_gene_family(gene_upper, uniprot_function, go_terms)

            # Map from our classification system to learning folder names
            # RESPECT ALL THE BRILLIANT CATEGORIES I BUILT! 🧬✨
            family_mapping = {
                'COLLAGEN_FIBRILLAR': 'collagen_fibrillar',
                'COLLAGEN_NETWORK': 'collagen_network',
                'COLLAGEN_ANCHORING': 'collagen_anchoring',
                'COLLAGEN_FACIT': 'collagen_facit',
                'FIBRILLIN': 'elastin_fibrillin',
                'ELASTIN': 'elastin_fibrillin',
                'ION_CHANNEL': 'ion_channel',
                'TUMOR_SUPPRESSOR': 'tumor_suppressor',
                'ONCOGENE': 'tumor_suppressor',
                'METABOLIC_ENZYME': 'metabolic_enzyme',
                'TRANSPORTER': 'transporter',
                'SCAFFOLD_ADAPTOR': 'scaffold_adaptor',
                'SIGNALING_REGULATOR': 'signaling_regulator',
                'INTERMEDIATE_FILAMENT': 'cytoskeleton',
                'CYTOSKELETON_POLYMER': 'cytoskeleton',
                'MOTOR_PROTEIN': 'motor_protein',
                'MUSCULAR_DYSTROPHY': 'muscular_dystrophy',
                'RTK_MAPK': 'signaling_regulator',
                'TRANSCRIPTION_FACTOR': 'transcription_factor',
                'RIBOSOMAL_PROTEIN': 'ribosomal_protein',
                'NEGATIVE_REGULATOR': 'signaling_regulator',
                'STRUCTURAL': 'structural',
                'AUTOSOMAL_RECESSIVE': 'autosomal_recessive',
                'GENERAL': 'general'
            }

            learning_family = family_mapping.get(family_classification)

            if learning_family:
                print(f"   ✅ {gene_upper} → {family_classification} → {learning_family}")
                # Cache the result
                self.gene_classification_cache[gene_upper] = learning_family
                return learning_family
            else:
                print(f"   ⚠️ {gene_upper} classified as {family_classification} (no learning folder)")
                return None

        except Exception as e:
            print(f"   ❌ Error classifying {gene_upper}: {e}")
            return None
    
    def standardize_filename(self, filename: str, gene: str) -> str:
        """Standardize filename to gene_benign.tsv or gene_pathogenic.tsv format"""
        filename_lower = filename.lower()
        gene_lower = gene.lower()
        
        # Determine if benign or pathogenic
        if any(word in filename_lower for word in ['benign', 'lb']):
            suffix = 'benign'
        elif any(word in filename_lower for word in ['patho', 'lp', 'pathogenic']):
            suffix = 'pathogenic'
        else:
            # Default based on common patterns
            suffix = 'unknown'
        
        # Change extension to .tsv (our standard)
        return f"{gene_lower}_{suffix}.tsv"
    
    def organize_files(self, dry_run: bool = False) -> Dict:
        """Organize all files in learning directory"""
        results = {
            'moved': [],
            'skipped': [],
            'errors': []
        }
        
        print("🚀 STARTING FILE ORGANIZATION")
        print("=" * 50)
        
        # Get all files in learning directory (not in subdirectories)
        files_to_organize = []
        for item in self.learning_dir.iterdir():
            if item.is_file() and item.suffix in ['.csv', '.tsv'] and item.name != 'README.md':
                files_to_organize.append(item)
        
        if not files_to_organize:
            print("⚠️  No files found to organize")
            return results
        
        print(f"📊 Found {len(files_to_organize)} files to organize")
        print()
        
        for file_path in files_to_organize:
            filename = file_path.name
            print(f"🔍 Processing: {filename}")
            
            # Extract gene name
            gene = self.extract_gene_from_filename(filename)
            if not gene:
                print(f"   ❌ Could not extract gene name")
                results['skipped'].append(filename)
                continue
            
            print(f"   🧬 Gene: {gene}")
            
            # Classify gene into family
            family = self.classify_gene(gene)
            if not family:
                print(f"   ❌ Could not classify gene into family")
                results['skipped'].append(filename)
                continue
            
            print(f"   📁 Family: {family}")
            
            # Create target directory
            target_dir = self.learning_dir / family
            if not dry_run:
                target_dir.mkdir(exist_ok=True)
            
            # Standardize filename
            new_filename = self.standardize_filename(filename, gene)
            target_path = target_dir / new_filename
            
            print(f"   📝 New name: {new_filename}")
            print(f"   🎯 Target: {target_path}")
            
            # Move file
            if not dry_run:
                try:
                    shutil.move(str(file_path), str(target_path))
                    print(f"   ✅ MOVED!")
                    results['moved'].append({
                        'original': filename,
                        'new_path': str(target_path),
                        'gene': gene,
                        'family': family
                    })
                except Exception as e:
                    print(f"   ❌ ERROR: {e}")
                    results['errors'].append({
                        'file': filename,
                        'error': str(e)
                    })
            else:
                print(f"   🔍 DRY RUN - would move here")
                results['moved'].append({
                    'original': filename,
                    'new_path': str(target_path),
                    'gene': gene,
                    'family': family
                })
            
            print()
        
        return results
    
    def print_summary(self, results: Dict):
        """Print organization summary"""
        print("📊 ORGANIZATION SUMMARY")
        print("=" * 50)
        
        moved = results['moved']
        skipped = results['skipped']
        errors = results['errors']
        
        print(f"✅ Successfully organized: {len(moved)} files")
        print(f"⏭️  Skipped: {len(skipped)} files")
        print(f"❌ Errors: {len(errors)} files")
        
        if moved:
            print("\n🎯 ORGANIZED FILES:")
            for item in moved:
                print(f"   {item['original']} → {item['family']}/{item['gene']}_{item['new_path'].split('_')[-1]}")
        
        if skipped:
            print(f"\n⏭️  SKIPPED FILES:")
            for filename in skipped:
                print(f"   {filename}")
        
        if errors:
            print(f"\n❌ ERROR FILES:")
            for item in errors:
                print(f"   {item['file']}: {item['error']}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="🔥💜 Learning File Organizer")
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be moved without actually moving files')
    
    args = parser.parse_args()
    
    organizer = LearningFileOrganizer()
    
    if args.dry_run:
        print("🔍 DRY RUN MODE - No files will be moved")
        print()
    
    results = organizer.organize_files(dry_run=args.dry_run)
    organizer.print_summary(results)
    
    if not args.dry_run and results['moved']:
        print(f"\n🚀 Ready to train! Run: python3 train_families.py")

if __name__ == "__main__":
    main()
