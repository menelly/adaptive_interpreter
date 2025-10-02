#!/usr/bin/env python3
"""
🔥💜 UNIFIED FAMILY-AWARE ML TRAINER 🚀
Learn AA coefficients AND conservation weights simultaneously!

Processes ClinVar data from /learning folders:
- genename_benign.tsv and genename_pathogenic.tsv
- Learns family-specific patterns for EVERYTHING
- Multi-target ML: pathogenicity + conservation weighting

Built by Ace (2025) for revolutionary ML training
Contact: ace@chaoschanneling.com
"""

import os
import glob
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
import joblib
import json
import re
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import pyBigWig

# Import our existing systems
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from nova_dn.amino_acid_props import AA_PROPS
from utils.gnomad_frequency_fetcher import GnomADFrequencyFetcher
from core_analyzers.plausibility_filter import classify_gene_family
from nova_dn.universal_context import UniversalContext
from analyzers.uniprot_mapper import UniProtMapper
from utils.genomic_to_protein import GenomicToProteinConverter

class UnifiedFamilyMLTrainer:
    """Revolutionary unified ML trainer for family-aware analysis"""
    
    def __init__(self, benign_freq_threshold: float = 0.04):
        self.learning_dir = Path("learning")
        self.models_dir = Path("resources/family_models")
        self.models_dir.mkdir(parents=True, exist_ok=True)

        # Frequency threshold for benign variants (default 4%)
        self.benign_freq_threshold = benign_freq_threshold

        # Initialize frequency fetcher for conservation data
        self.freq_fetcher = GnomADFrequencyFetcher()

        # Initialize GO classification system
        self.universal_context = UniversalContext()
        self.gene_classification_cache = {}

        # Initialize UniProt mapper for coordinate conversion
        self.uniprot_mapper = UniProtMapper()

        # Initialize genomic to protein converter
        self.genomic_converter = GenomicToProteinConverter()

        # Conservation data paths
        self.conservation_dir = Path("/home/Ace/conservation_data")
        self.phylop_file = self.conservation_dir / "hg38.phyloP470way.bw"
        self.phastcons_file = self.conservation_dir / "hg38.phastCons470way.bw"

        # Initialize conservation file handles (lazy loading)
        self._phylop_bw = None
        self._phastcons_bw = None

        print("🔥💜 UNIFIED FAMILY ML TRAINER INITIALIZED")
        print(f"📁 Learning directory: {self.learning_dir}")
        print(f"💾 Models directory: {self.models_dir}")
        print(f"🚨 Benign frequency filter: ≤{benign_freq_threshold*100:.1f}% (sanity check!)")
        print("🧬 Using real GO term classification system!")
        print(f"🧬 Conservation data: {self.conservation_dir}")
    
    def extract_variant_info(self, hgvs: str, gene_from_filename: str = None) -> Optional[Dict]:
        """Extract gene, transcript, and variant from HGVS with protein parsing"""

        # 🚀 NEW: Handle genomic HGVS format (from ClinVar)
        if self.genomic_converter.is_genomic_hgvs(hgvs):
            return self._extract_from_genomic_hgvs(hgvs, gene_from_filename)

        # Original protein HGVS handling
        # Extract gene name from parentheses
        gene_match = re.search(r'\(([^)]+)\)', hgvs)
        if not gene_match:
            return None
        gene = gene_match.group(1)

        # Look for protein change p.RefPosAlt (handle both 1-letter and 3-letter codes)
        protein_match = re.search(r'p\.([A-Z][a-z]{0,2})(\d+)([A-Z][a-z]{0,2})', hgvs)
        if protein_match:
            ref_aa_full, pos, alt_aa_full = protein_match.groups()

            # Convert 3-letter to 1-letter codes
            aa_map = {
                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
            }

            # Handle both 1-letter and 3-letter codes
            ref_aa = aa_map.get(ref_aa_full, ref_aa_full)
            alt_aa = aa_map.get(alt_aa_full, alt_aa_full)

            return {
                'gene': gene,
                'hgvs': hgvs,
                'position': int(pos),
                'ref_aa': ref_aa,
                'alt_aa': alt_aa,
                'variant_str': f"p.{ref_aa}{pos}{alt_aa}",
                'type': 'missense'
            }

        # If no protein change found, return basic info
        return {
            'gene': gene,
            'hgvs': hgvs,
            'position': None,
            'ref_aa': None,
            'alt_aa': None,
            'variant_str': None,
            'type': 'unknown'
        }

    def _extract_from_genomic_hgvs(self, hgvs: str, gene_from_filename: str = None) -> Optional[Dict]:
        """Extract variant info from genomic HGVS format"""

        # Use genomic converter to parse coordinates
        variant_info = self.genomic_converter.extract_variant_info_from_genomic(hgvs, gene_from_filename or "UNKNOWN_GENE")
        if not variant_info:
            return None

        # For genomic variants we do NOT fabricate protein AA fields.
        # We keep protein-position/AA fields as None so downstream calibrators
        # won't treat placeholders as real counts.

        return {
            'gene': gene_from_filename or variant_info['gene'],
            'hgvs': hgvs,
            'position': None,
            'ref_aa': None,
            'alt_aa': None,
            'variant_str': variant_info['variant_str'],
            'type': 'genomic_missense',
            'chromosome': variant_info['chromosome'],
            'genomic_position': variant_info['genomic_position'],
            'ref_allele': variant_info['ref_allele'],
            'alt_allele': variant_info['alt_allele']
        }

    def calculate_aa_features(self, ref_aa: str, alt_aa: str) -> Dict:
        """Calculate amino acid property features"""
        if ref_aa not in AA_PROPS or alt_aa not in AA_PROPS:
            return {}
            
        ref_props = AA_PROPS[ref_aa]
        alt_props = AA_PROPS[alt_aa]
        
        return {
            'grantham_distance': self.get_grantham_distance(ref_aa, alt_aa),
            'delta_hydrophobicity': abs(alt_props['hyd'] - ref_props['hyd']),
            'delta_volume': abs(alt_props['vol'] - ref_props['vol']),
            'delta_charge': abs(alt_props['chg'] - ref_props['chg']),
            'ref_hydrophobicity': ref_props['hyd'],
            'alt_hydrophobicity': alt_props['hyd'],
            'ref_volume': ref_props['vol'],
            'alt_volume': alt_props['vol'],
            'ref_charge': ref_props['chg'],
            'alt_charge': alt_props['chg'],
            'proline_involved': 1 if ref_aa == 'P' or alt_aa == 'P' else 0,
            'glycine_involved': 1 if ref_aa == 'G' or alt_aa == 'G' else 0,
            'cysteine_involved': 1 if ref_aa == 'C' or alt_aa == 'C' else 0
        }
    
    def get_grantham_distance(self, ref_aa: str, alt_aa: str) -> float:
        """Get Grantham distance between amino acids"""
        # Simplified Grantham matrix - in real implementation, use full matrix
        grantham_matrix = {
            ('R', 'H'): 29, ('R', 'K'): 26, ('R', 'D'): 96, ('R', 'E'): 54,
            ('H', 'R'): 29, ('H', 'K'): 32, ('H', 'D'): 81, ('H', 'E'): 40,
            ('G', 'A'): 60, ('G', 'S'): 56, ('G', 'P'): 42, ('G', 'V'): 109,
            ('P', 'L'): 98, ('P', 'S'): 74, ('P', 'A'): 27, ('P', 'G'): 42,
            # Add more as needed...
        }
        
        key = (ref_aa, alt_aa)
        reverse_key = (alt_aa, ref_aa)
        
        if key in grantham_matrix:
            return grantham_matrix[key]
        elif reverse_key in grantham_matrix:
            return grantham_matrix[reverse_key]
        else:
            # Default distance for unknown pairs
            return 50.0
    
    def get_real_domain_features(self, gene: str, position: int, ref_aa: str, alt_aa: str, hgvs: str = "") -> Dict:
        """Get REAL domain features from our existing sophisticated system"""
        try:
            # Get protein context from our existing system
            context = self.universal_context.get_context_for_protein(gene)

            if "error" in context:
                # Fallback to basic features if no context available
                return {
                    'domain_type': 'UNKNOWN',
                    'active_site_proximity': 0.0,
                    'binding_site_proximity': 0.0,
                    'structural_importance': 0.5,
                    'domain_multiplier': 1.0,
                    'phylop_score': 0.0,
                    'phastcons_score': 0.5,
                    'gerp_score': 0.0,
                    'in_gly_x_y_repeat': 0,
                    'in_active_site': 0,
                    'in_binding_site': 0,
                    'in_structural_domain': 0
                }

            # Extract domain information for this position
            domains = context.get('domains', [])

            # Find which domain this position is in
            current_domain = None
            for domain in domains:
                start = domain.get('start', 0)
                end = domain.get('end', 0)
                if start <= position <= end:
                    current_domain = domain
                    break

            # Calculate domain-specific features
            domain_features = {
                'domain_type': current_domain.get('description', 'UNKNOWN') if current_domain else 'UNKNOWN',
                'active_site_proximity': self.calculate_site_proximity(position, context.get('active_sites', [])),
                'binding_site_proximity': self.calculate_site_proximity(position, context.get('binding_sites', [])),
                'structural_importance': self.calculate_structural_importance(position, domains),
                'domain_multiplier': self.get_domain_multiplier(current_domain, gene),
                'phylop_score': self.get_real_phylop_score_from_hgvs(hgvs),
                'phastcons_score': self.get_real_phastcons_score_from_hgvs(hgvs),
                'gerp_score': 0.0,  # GERP not available in current dataset
                # Special region flags
                'in_gly_x_y_repeat': 1 if self.is_in_gly_x_y_repeat(position, gene) else 0,
                'in_active_site': 1 if self.is_in_active_site(position, context.get('active_sites', [])) else 0,
                'in_binding_site': 1 if self.is_in_binding_site(position, context.get('binding_sites', [])) else 0,
                'in_structural_domain': 1 if current_domain and 'structural' in current_domain.get('description', '').lower() else 0
            }

            return domain_features

        except Exception as e:
            print(f"   ⚠️ Error getting domain features for {gene} position {position}: {e}")
            # Return safe defaults
            return {
                'domain_type': 'ERROR',
                'active_site_proximity': 0.0,
                'binding_site_proximity': 0.0,
                'structural_importance': 0.5,
                'domain_multiplier': 1.0,
                'phylop_score': 0.0,
                'phastcons_score': 0.5,
                'gerp_score': 0.0,
                'in_gly_x_y_repeat': 0,
                'in_active_site': 0,
                'in_binding_site': 0,
                'in_structural_domain': 0
            }

    def calculate_site_proximity(self, position: int, sites: List[Dict]) -> float:
        """Calculate proximity to active/binding sites"""
        if not sites:
            return 0.0

        min_distance = float('inf')
        for site in sites:
            if isinstance(site, dict):
                site_pos = site.get('position', site.get('start', 0))
            else:
                site_pos = site  # Assume it's just a position number

            if site_pos > 0:
                distance = abs(position - site_pos)
                min_distance = min(min_distance, distance)

        if min_distance == float('inf'):
            return 0.0

        # Convert distance to proximity score (closer = higher score)
        if min_distance == 0:
            return 1.0  # Direct hit
        elif min_distance <= 5:
            return 0.8  # Very close
        elif min_distance <= 10:
            return 0.5  # Close
        elif min_distance <= 20:
            return 0.2  # Moderate distance
        else:
            return 0.0  # Far away

    def calculate_structural_importance(self, position: int, domains: List[Dict]) -> float:
        """Calculate structural importance based on domain context"""
        importance = 0.5  # Base importance

        for domain in domains:
            start = domain.get('start', 0)
            end = domain.get('end', 0)

            if start <= position <= end:
                desc = domain.get('description', '').lower()

                # High importance domains
                if any(term in desc for term in ['catalytic', 'active', 'kinase', 'enzyme']):
                    importance = max(importance, 0.9)
                elif any(term in desc for term in ['binding', 'interface', 'interaction']):
                    importance = max(importance, 0.8)
                elif any(term in desc for term in ['structural', 'domain', 'fold']):
                    importance = max(importance, 0.7)
                elif any(term in desc for term in ['repeat', 'motif']):
                    importance = max(importance, 0.6)

        return importance

    def get_domain_multiplier(self, domain: Dict, gene: str) -> float:
        """Get domain-specific multiplier based on domain type and gene family"""
        if not domain:
            return 1.0

        multiplier = 1.0
        desc = domain.get('description', '').lower()

        # Gene family specific adjustments
        if gene.startswith('COL'):  # Collagen genes
            if 'triple-helical' in desc or 'collagenous' in desc:
                multiplier *= 1.5  # Critical for collagen structure
            elif 'propeptide' in desc:
                multiplier *= 0.4  # Gets cleaved off

        # Universal domain importance
        if any(term in desc for term in ['catalytic', 'active', 'kinase']):
            multiplier *= 1.4
        elif any(term in desc for term in ['binding', 'interface']):
            multiplier *= 1.3
        elif 'signal' in desc:
            multiplier *= 0.3  # Signal peptides get cleaved
        elif 'propeptide' in desc:
            multiplier *= 0.4  # Propeptides get cleaved

        return multiplier

    def is_in_gly_x_y_repeat(self, position: int, gene: str) -> bool:
        """Check if position is in a Gly-X-Y repeat (collagen-specific)"""
        if not gene.startswith('COL'):
            return False

        # Simplified check - in real collagen, every 3rd position starting from ~20 is glycine
        # This is a rough approximation for collagen triple helix regions
        if position < 20:  # Before triple helix usually starts
            return False

        # Check if position aligns with glycine positions in Gly-X-Y repeat
        # Glycine is typically at positions 20, 23, 26, 29, etc.
        return (position - 20) % 3 == 0

    def is_in_active_site(self, position: int, active_sites: List[Dict]) -> bool:
        """Check if position is in an active site"""
        for site in active_sites:
            if isinstance(site, dict):
                site_pos = site.get('position', site.get('start', 0))
                site_end = site.get('end', site_pos)
                if site_pos <= position <= site_end:
                    return True
            elif isinstance(site, int):
                if site == position:
                    return True
        return False

    def is_in_binding_site(self, position: int, binding_sites: List[Dict]) -> bool:
        """Check if position is in a binding site"""
        for site in binding_sites:
            if isinstance(site, dict):
                site_pos = site.get('position', site.get('start', 0))
                site_end = site.get('end', site_pos)
                if site_pos <= position <= site_end:
                    return True
            elif isinstance(site, int):
                if site == position:
                    return True
        return False

    def classify_gene_with_go_terms(self, gene: str) -> Optional[str]:
        """Classify gene using our GO term classification system"""
        gene_upper = gene.upper()

        # Check cache first
        if gene_upper in self.gene_classification_cache:
            return self.gene_classification_cache[gene_upper]

        try:
            # Get protein context (UniProt function + GO terms)
            context = self.universal_context.get_context_for_protein(gene_upper)

            if "error" in context:
                print(f"   ⚠️ Could not get context for {gene_upper}: {context['error']}")
                return None

            # Extract function and GO terms
            uniprot_function = context.get('function', '')
            go_terms = context.get('go_terms', [])

            # Use our existing classification system
            family_classification = classify_gene_family(gene_upper, uniprot_function, go_terms)

            # Map from our classification system to learning folder names
            family_mapping = {
                'COLLAGEN_FIBRILLAR': 'collagen_fibrillar',
                'COLLAGEN_NETWORK': 'collagen_fibrillar',
                'COLLAGEN_ANCHORING': 'collagen_fibrillar',
                'COLLAGEN_FACIT': 'collagen_fibrillar',
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
                'MOTOR_PROTEIN': 'cytoskeleton',
                'MUSCULAR_DYSTROPHY': 'cytoskeleton',
                'RTK_MAPK': 'signaling_regulator',
                'TRANSCRIPTION_FACTOR': 'signaling_regulator',
                'RIBOSOMAL_PROTEIN': 'metabolic_enzyme',
                'NEGATIVE_REGULATOR': 'signaling_regulator',
                'STRUCTURAL': 'cytoskeleton',
                'GENERAL': None
            }

            learning_family = family_mapping.get(family_classification)

            # Cache the result
            self.gene_classification_cache[gene_upper] = learning_family
            return learning_family

        except Exception as e:
            print(f"   ❌ Error classifying {gene_upper}: {e}")
            return None

    def get_real_phylop_score_from_hgvs(self, hgvs: str) -> float:
        """Get real phyloP conservation score from ClinVar HGVS genomic coordinates"""
        try:
            # Lazy load phyloP file
            if self._phylop_bw is None:
                if not self.phylop_file.exists():
                    print(f"⚠️ phyloP file not found: {self.phylop_file}")
                    return 0.0
                self._phylop_bw = pyBigWig.open(str(self.phylop_file))

            # Extract genomic coordinates from HGVS
            genomic_coords = self._parse_genomic_hgvs(hgvs)
            if not genomic_coords:
                return 0.0

            chrom, position = genomic_coords

            # Get score at the exact position (BigWig uses 0-based coordinates)
            score = self._phylop_bw.values(chrom, position - 1, position)

            if score and len(score) > 0 and score[0] is not None:
                return float(score[0])

            return 0.0

        except Exception as e:
            print(f"⚠️ Error getting phyloP score from HGVS {hgvs}: {e}")
            return 0.0

    def get_real_phastcons_score_from_hgvs(self, hgvs: str) -> float:
        """Get real phastCons conservation score from ClinVar HGVS genomic coordinates"""
        try:
            # Lazy load phastCons file
            if self._phastcons_bw is None:
                if not self.phastcons_file.exists():
                    print(f"⚠️ phastCons file not found: {self.phastcons_file}")
                    return 0.0
                self._phastcons_bw = pyBigWig.open(str(self.phastcons_file))

            # Extract genomic coordinates from HGVS
            genomic_coords = self._parse_genomic_hgvs(hgvs)
            if not genomic_coords:
                return 0.0

            chrom, position = genomic_coords

            # Get score at the exact position (BigWig uses 0-based coordinates)
            score = self._phastcons_bw.values(chrom, position - 1, position)

            if score and len(score) > 0 and score[0] is not None:
                return float(score[0])

            return 0.0

        except Exception as e:
            print(f"⚠️ Error getting phastCons score from HGVS {hgvs}: {e}")
            return 0.0

    def _parse_genomic_hgvs(self, hgvs: str) -> Optional[Tuple[str, int]]:
        """Parse genomic coordinates from ClinVar HGVS notation"""
        try:
            # Look for genomic HGVS like: NC_000017.11:g.50183779T>A
            # Extract chromosome and position

            # Pattern for NC_XXXXXX.XX:g.POSITION
            genomic_match = re.search(r'NC_(\d+)(?:\.\d+)?:g\.(\d+)', hgvs)
            if genomic_match:
                raw_num = genomic_match.group(1)
                chrom_num = raw_num.lstrip('0') or raw_num  # Remove leading zeros but keep raw if all zeros
                position = int(genomic_match.group(2))

                # Map special chromosomes
                if chrom_num in ("23",):
                    chrom = "chrX"
                elif chrom_num in ("24",):
                    chrom = "chrY"
                elif chrom_num in ("12920",):  # NC_012920.* mitochondrial
                    chrom = "chrM"
                else:
                    chrom = f"chr{chrom_num}"

                return (chrom, position)

            # Pattern for direct chr format: supports chr1..chr22, chrX, chrY, chrM
            chr_match = re.search(r'chr([0-9XYM]+):(\d+)', hgvs, re.IGNORECASE)
            if chr_match:
                chrom_id = chr_match.group(1).upper()
                position = int(chr_match.group(2))
                if chrom_id == "23":
                    chrom = "chrX"
                elif chrom_id == "24":
                    chrom = "chrY"
                else:
                    chrom = f"chr{chrom_id}"
                return (chrom, position)

            return None

        except Exception as e:
            print(f"⚠️ Error parsing genomic HGVS {hgvs}: {e}")
            return None

    def process_family_data(self, family_name: str) -> pd.DataFrame:
        """Process all TSV files for a gene family"""
        family_dir = self.learning_dir / family_name
        
        if not family_dir.exists():
            print(f"⚠️  Family directory not found: {family_dir}")
            return pd.DataFrame()
        
        all_data = []
        skipped_high_freq = 0
        skipped_synonymous = 0
        skipped_splice = 0
        skipped_nonsense = 0
        skipped_frameshift = 0
        skipped_unparseable = 0

        # Process all TSV files in family directory
        for tsv_file in family_dir.glob("*.tsv"):
            print(f"📊 Processing {tsv_file.name}...")

            try:
                # Try to detect delimiter and read file
                # First, try to read and detect the actual format
                with open(tsv_file, 'r') as f:
                    first_line = f.readline().strip()

                # Check if it's properly quoted CSV or malformed
                if first_line.startswith('"HGVS","dbSNP"'):
                    # Properly quoted CSV
                    df = pd.read_csv(tsv_file)
                elif first_line.startswith('HGVS,"dbSNP"'):
                    # Malformed CSV - mixed quoting
                    df = pd.read_csv(tsv_file, quotechar='"', skipinitialspace=True)
                elif '\t' in first_line:
                    # TSV format
                    df = pd.read_csv(tsv_file, sep='\t')
                else:
                    # Default CSV
                    df = pd.read_csv(tsv_file)

                # Determine pathogenicity from filename
                is_pathogenic = 'pathogenic' in tsv_file.name.lower()
                is_benign = 'benign' in tsv_file.name.lower()
                pathogenicity_score = 1.0 if is_pathogenic else 0.0

                # Extract gene name from filename (e.g., "myo7a_pathogenic.tsv" -> "MYO7A")
                gene_from_filename = tsv_file.stem.split('_')[0].upper()

                print(f"   📋 Columns found: {df.columns.tolist()}")
                print(f"   📊 Rows: {len(df)}")

                for _, row in df.iterrows():
                    hgvs = row['HGVS']
                    frequency = row.get('gnomAD frequency', 0.0)

                    # 🚨 COMPREHENSIVE SANITY CHECKS

                    # Skip synonymous variants (no AA change)
                    if '=' in hgvs:
                        skipped_synonymous += 1
                        continue

                    # Skip splice variants (we're not set up for splice analysis yet)
                    if '+' in hgvs or '-' in hgvs:
                        skipped_splice += 1
                        continue

                    # Skip nonsense/stop variants (we're not set up for Ter analysis yet)
                    if 'Ter' in hgvs or '*' in hgvs:
                        skipped_nonsense += 1
                        continue

                    # Skip frameshift variants (completely different analysis needed)
                    if 'del' in hgvs or 'ins' in hgvs or 'dup' in hgvs:
                        skipped_frameshift += 1
                        continue

                    # 🚨 FREQUENCY SANITY CHECK for benign variants
                    if is_benign and frequency > self.benign_freq_threshold:
                        skipped_high_freq += 1
                        continue  # Skip common "benign" variants (they're just normal!)

                    # Extract variant info (pass gene name from filename for genomic HGVS)
                    variant_info = self.extract_variant_info(hgvs, gene_from_filename)
                    if not variant_info:
                        skipped_unparseable += 1
                        continue

                    # Extract real protein variant info
                    ref_aa = variant_info.get('ref_aa')
                    alt_aa = variant_info.get('alt_aa')
                    position = variant_info.get('position')

                    # Skip if we couldn't parse the protein change
                    if not ref_aa or not alt_aa or not position:
                        skipped_unparseable += 1
                        continue
                    
                    # Calculate features
                    aa_features = self.calculate_aa_features(ref_aa, alt_aa)
                    domain_features = self.get_real_domain_features(variant_info['gene'], position, ref_aa, alt_aa, hgvs)
                    
                    # Combine all features
                    feature_row = {
                        'gene': variant_info['gene'],
                        'family': family_name,
                        'hgvs': hgvs,
                        'frequency': frequency,
                        'pathogenicity_score': pathogenicity_score,
                        # Explicit AA fields needed by the calibrator
                        'ref_aa': ref_aa,
                        'alt_aa': alt_aa,
                        'position': position,
                        **aa_features,
                        **domain_features
                    }

                    all_data.append(feature_row)

            except Exception as e:
                print(f"❌ Error processing {tsv_file}: {e}")
                continue
        
        if not all_data:
            print(f"⚠️  No data processed for family {family_name}")
            return pd.DataFrame()

        df = pd.DataFrame(all_data)
        print(f"✅ Processed {len(df)} variants for {family_name}")

        # Report all the sanity check skips
        total_skipped = (skipped_high_freq + skipped_synonymous + skipped_splice +
                        skipped_nonsense + skipped_frameshift + skipped_unparseable)

        if total_skipped > 0:
            print(f"🚨 SANITY CHECK REPORT - Skipped {total_skipped} variants:")
            if skipped_synonymous > 0:
                print(f"   📝 {skipped_synonymous} synonymous variants (= in HGVS)")
            if skipped_splice > 0:
                print(f"   🧬 {skipped_splice} splice variants (+/- in HGVS)")
            if skipped_nonsense > 0:
                print(f"   🛑 {skipped_nonsense} nonsense variants (Ter/* in HGVS)")
            if skipped_frameshift > 0:
                print(f"   🔄 {skipped_frameshift} frameshift variants (del/ins/dup)")
            if skipped_high_freq > 0:
                print(f"   📊 {skipped_high_freq} high-frequency benign variants (>{self.benign_freq_threshold*100:.1f}%)")
            if skipped_unparseable > 0:
                print(f"   ❓ {skipped_unparseable} unparseable HGVS formats")
            print(f"   ✅ Reason: Not set up for these variant types yet!")

        return df
    
    def train_family_model(self, family_name: str, df: pd.DataFrame) -> Dict:
        """Train unified model for a gene family"""
        if df.empty:
            print(f"⚠️  No data to train for {family_name}")
            return {}
        
        print(f"🧠 Training unified model for {family_name}...")
        
        # Prepare features
        feature_columns = [
            'grantham_distance', 'delta_hydrophobicity', 'delta_volume', 'delta_charge',
            'ref_hydrophobicity', 'alt_hydrophobicity', 'ref_volume', 'alt_volume',
            'ref_charge', 'alt_charge', 'proline_involved', 'glycine_involved', 
            'cysteine_involved', 'phylop_score', 'phastcons_score', 'gerp_score',
            'frequency'
        ]
        
        # Filter to available columns
        available_features = [col for col in feature_columns if col in df.columns]
        X = df[available_features].fillna(0)
        
        # Target: pathogenicity score
        y = df['pathogenicity_score']
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42, stratify=y
        )
        
        # Scale features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Train model
        model = RandomForestRegressor(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            n_jobs=-1
        )
        
        model.fit(X_train_scaled, y_train)
        
        # Evaluate
        y_pred = model.predict(X_test_scaled)
        mse = mean_squared_error(y_test, y_pred)
        r2 = r2_score(y_test, y_pred)
        
        print(f"📊 {family_name} Model Performance:")
        print(f"   MSE: {mse:.4f}")
        print(f"   R²: {r2:.4f}")
        
        # Feature importance
        feature_importance = dict(zip(available_features, model.feature_importances_))
        
        # Save model and scaler
        model_path = self.models_dir / f"{family_name}_unified_model.joblib"
        scaler_path = self.models_dir / f"{family_name}_unified_scaler.joblib"
        
        joblib.dump(model, model_path)
        joblib.dump(scaler, scaler_path)
        
        # Save metadata
        metadata = {
            'family': family_name,
            'n_samples': len(df),
            'n_features': len(available_features),
            'feature_columns': available_features,
            'feature_importance': feature_importance,
            'performance': {
                'mse': mse,
                'r2': r2
            },
            'model_path': str(model_path),
            'scaler_path': str(scaler_path)
        }
        
        metadata_path = self.models_dir / f"{family_name}_unified_metadata.json"
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        print(f"💾 Saved model: {model_path}")
        print(f"💾 Saved metadata: {metadata_path}")
        
        return metadata
    
    def train_all_families(self) -> Dict:
        """Train models for all families in learning directory"""
        results = {}
        
        print("🚀 STARTING UNIFIED FAMILY ML TRAINING")
        print("=" * 50)
        
        for family_dir in self.learning_dir.iterdir():
            if not family_dir.is_dir():
                continue
                
            family_name = family_dir.name
            print(f"\n🧬 Processing family: {family_name}")
            
            # Process family data
            df = self.process_family_data(family_name)
            
            if not df.empty:
                # Train model
                metadata = self.train_family_model(family_name, df)
                results[family_name] = metadata
            else:
                print(f"⏭️  Skipping {family_name} - no data")
        
        print("\n🎉 TRAINING COMPLETE!")
        print("=" * 50)
        
        # Save overall results
        results_path = self.models_dir / "training_results.json"
        with open(results_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"📊 Training summary saved: {results_path}")
        return results

def main():
    """Main training function"""
    import argparse

    parser = argparse.ArgumentParser(description="🔥💜 Unified Family ML Trainer")
    parser.add_argument('--benign-freq-threshold', type=float, default=0.04,
                       help='Maximum frequency for benign variants (default: 0.04 = 4%%)')

    args = parser.parse_args()

    trainer = UnifiedFamilyMLTrainer(benign_freq_threshold=args.benign_freq_threshold)
    results = trainer.train_all_families()

    print(f"\n✅ Trained models for {len(results)} families!")
    for family, metadata in results.items():
        r2 = metadata.get('performance', {}).get('r2', 0)
        n_samples = metadata.get('n_samples', 0)
        print(f"   {family}: R² = {r2:.3f} ({n_samples} samples)")

if __name__ == "__main__":
    main()
