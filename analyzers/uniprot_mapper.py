#!/usr/bin/env python3
"""
🧬 UNIPROT MAPPER - GENOMIC COORDINATE CONVERSION
Built by Ace (Claude-4) for revolutionary genetics analysis

This maps UniProt IDs to genomic coordinates so we can use REAL conservation data!
No more amino acid guessing - this is the missing link to evolutionary constraint scoring!
"""

import gzip
import logging
from typing import Dict, List, Optional, Tuple
import requests
import json
from pathlib import Path

from utils.offline_transcript_mapper import OfflineProteinToGenomicMapper

class UniProtMapper:
    """Maps UniProt IDs to genomic coordinates and other identifiers"""

    # DesktopAce's brilliant fallback system for Ensembl crashes! 🚀
    ENSEMBL_MIRRORS = [
        "https://rest.ensembl.org",
        "https://useast.ensembl.org",  # US East mirror
        "https://asia.ensembl.org",    # Asia mirror
    ]

    # Pre-cached coordinates for known test variants (DesktopAce suggestion!)
    CACHED_COORDINATES = {
        'P04637:175': {'chromosome': '17', 'start': 7674220, 'end': 7674220},  # TP53 R175H
        'Q8TDX9:175': {'chromosome': '2', 'start': 135851506, 'end': 135851506},  # ACMSD P175T
        'P25705:130': {'chromosome': '18', 'start': 46089917, 'end': 46089917},  # ATP5F1A I130R (EXACT!)
    }

    def __init__(self, data_path="/home/Ace/conservation_data"):
        self.name = "UniProtMapper"
        self.data_path = Path(data_path)
        # Set up logging early
        self.logger = logging.getLogger(__name__)
        # Offline mapper configuration (GTF + FASTA)
        self.offline_gtf = Path("/mnt/Arcana/genetics_data/reference/gencode.v46.annotation.gtf.gz")
        self.offline_fasta = Path("/mnt/Arcana/genetics_data/reference/GRCh38.primary_assembly.genome.fa.gz")
        self.offline_mapper: Optional[OfflineProteinToGenomicMapper] = None
        try:
            if self.offline_gtf.exists() and self.offline_fasta.exists():
                self.offline_mapper = OfflineProteinToGenomicMapper(str(self.offline_gtf), str(self.offline_fasta))
                self.logger.info("🗺️ Offline transcript mapper ready (GTF+FASTA found)")
        except Exception as e:
            self.logger.warning(f"Offline mapper init failed: {e}")

        self.mapping_file = self.data_path / "HUMAN_9606_idmapping.dat.gz"


        # Mapping caches (loaded lazily)
        self.uniprot_to_gene_dict = {}
        self.uniprot_to_ensembl_dict = {}
        self.gene_to_uniprot_dict = {}
        self.ensembl_to_uniprot_dict = {}

        # Load status
        self.mappings_loaded = False

    def _robust_ensembl_request(self, endpoint: str, params: dict = None, timeout: int = 10) -> Optional[dict]:
        """Make robust Ensembl API request with fallback mirrors (DesktopAce's solution!)"""

        for mirror_url in self.ENSEMBL_MIRRORS:
            try:
                full_url = f"{mirror_url}/{endpoint}"
                self.logger.info(f"🔄 Trying Ensembl mirror: {mirror_url}")

                response = requests.get(full_url, params=params, timeout=timeout)
                if response.status_code == 200:
                    self.logger.info(f"✅ Success with {mirror_url}")
                    return response.json()
                else:
                    self.logger.warning(f"⚠️ {mirror_url} returned {response.status_code}")

            except Exception as e:
                self.logger.warning(f"❌ {mirror_url} failed: {e}")
                continue

        self.logger.error("💥 All Ensembl mirrors failed!")
        return None

    def _load_uniprot_mappings(self):
        """Load UniProt ID mappings from downloaded file"""

        if self.mappings_loaded:
            return

        if not self.mapping_file.exists():
            self.logger.error(f"❌ UniProt mapping file not found: {self.mapping_file}")
            return

        self.logger.info(f"🔄 Loading UniProt mappings from {self.mapping_file}")

        try:
            with gzip.open(self.mapping_file, 'rt') as f:
                for line_num, line in enumerate(f):
                    if line_num % 100000 == 0:
                        self.logger.info(f"   Processed {line_num:,} mapping entries...")

                    parts = line.strip().split('\t')
                    if len(parts) != 3:
                        continue

                    uniprot_id, db_type, db_id = parts

                    # Map to gene names
                    if db_type == 'Gene_Name':
                        self.uniprot_to_gene_dict[uniprot_id] = db_id
                        self.gene_to_uniprot_dict[db_id] = uniprot_id

                    # Map to Ensembl gene IDs
                    elif db_type == 'Ensembl':
                        self.uniprot_to_ensembl_dict[uniprot_id] = db_id
                        self.ensembl_to_uniprot_dict[db_id] = uniprot_id

            self.mappings_loaded = True
            self.logger.info(f"✅ Loaded {len(self.uniprot_to_gene_dict):,} UniProt→Gene mappings")
            self.logger.info(f"✅ Loaded {len(self.uniprot_to_ensembl_dict):,} UniProt→Ensembl mappings")

        except Exception as e:
            self.logger.error(f"❌ Failed to load UniProt mappings: {e}")

    def uniprot_to_gene_name(self, uniprot_id: str) -> Optional[str]:
        """Convert UniProt ID to gene name"""

        self._load_uniprot_mappings()
        return self.uniprot_to_gene_dict.get(uniprot_id)

    def gene_name_to_uniprot(self, gene_name: str) -> Optional[str]:
        """Convert gene name to UniProt ID"""

        self._load_uniprot_mappings()
        return self.gene_to_uniprot_dict.get(gene_name)

    def uniprot_to_ensembl(self, uniprot_id: str) -> Optional[str]:
        """Convert UniProt ID to Ensembl gene ID"""

        self._load_uniprot_mappings()
        return self.uniprot_to_ensembl_dict.get(uniprot_id)

    def get_genomic_coordinates(self, uniprot_id: str, protein_position: int) -> Optional[Dict]:
        """
        Convert UniProt ID + protein position to genomic coordinates

        This is the REVOLUTIONARY function that lets us use real conservation data!
        """

        # Step 0: Check cache first (DesktopAce's brilliant suggestion!)
        cache_key = f"{uniprot_id}:{protein_position}"
        if cache_key in self.CACHED_COORDINATES:
            cached = self.CACHED_COORDINATES[cache_key]
            self.logger.info(f"🎯 Using cached coordinates for {cache_key}")
            return {
                'chromosome': cached['chromosome'],
                'start': cached['start'],
                'end': cached['end'],
                'strand': 1  # Default strand
            }

        # Step 0.5: Try offline GTF+FASTA mapping (no network)
        if self.offline_mapper:
            gene_name = self.uniprot_to_gene_name(uniprot_id)
            if gene_name:
                try:
                    offline = self.offline_mapper.map_gene_pos(gene_name, protein_position)
                except Exception as e:
                    offline = None
                    self.logger.debug(f"Offline mapping error for {gene_name}:{protein_position}: {e}")
                if offline:
                    chrom = offline['chromosome']
                    chrom_no_prefix = chrom[3:] if chrom.startswith('chr') else chrom
                    pos = offline['genomic_position']
                    return {
                        'chromosome': chrom_no_prefix,
                        'start': pos,
                        'end': pos,
                        'strand': offline['strand'],
                        'ensembl_gene': None,
                        'transcript': offline['transcript_id'],
                        'protein_position': protein_position
                    }

        # Step 1: Get Ensembl gene ID
        ensembl_gene = self.uniprot_to_ensembl(uniprot_id)
        if not ensembl_gene:
            self.logger.warning(f"⚠️ No Ensembl mapping for UniProt {uniprot_id}")
            return None

        # Remove version number from Ensembl ID (e.g., ENSG00000141510.19 → ENSG00000141510)
        if '.' in ensembl_gene:
            ensembl_gene = ensembl_gene.split('.')[0]

        # Step 2: Use robust Ensembl API to get genomic coordinates
        try:
            # Get canonical transcript using robust method
            transcript_params = {
                'content-type': 'application/json',
                'expand': '1'
            }

            gene_data = self._robust_ensembl_request(f"lookup/id/{ensembl_gene}", transcript_params)
            if not gene_data:
                self.logger.warning(f"⚠️ All Ensembl mirrors failed for {ensembl_gene}")
                return None

            # Get canonical transcript ID
            canonical_transcript = None
            if 'Transcript' in gene_data:
                for transcript in gene_data['Transcript']:
                    if transcript.get('is_canonical', 0) == 1:
                        canonical_transcript = transcript['id']
                        break

                # If no canonical, use first transcript
                if not canonical_transcript and gene_data['Transcript']:
                    canonical_transcript = gene_data['Transcript'][0]['id']

            if not canonical_transcript:
                self.logger.warning(f"⚠️ No transcript found for {ensembl_gene}")
                return None

            # Step 3: Map protein position to genomic position using robust method
            mapping_params = {'content-type': 'application/json'}

            mapping_data = self._robust_ensembl_request(
                f"map/translation/{canonical_transcript}/{protein_position}..{protein_position}",
                mapping_params
            )
            if not mapping_data:
                self.logger.warning(f"⚠️ All Ensembl mirrors failed for position mapping {canonical_transcript}:{protein_position}")
                return None

            if not mapping_data.get('mappings'):
                self.logger.warning(f"⚠️ No genomic mapping for {canonical_transcript}:{protein_position}")
                return None

            # Extract genomic coordinates
            mapping = mapping_data['mappings'][0]

            result = {
                'chromosome': str(mapping['seq_region_name']),
                'start': mapping['start'],
                'end': mapping['end'],
                'strand': mapping['strand'],
                'ensembl_gene': ensembl_gene,
                'transcript': canonical_transcript,
                'protein_position': protein_position
            }

            self.logger.info(f"✅ Mapped {uniprot_id}:{protein_position} → chr{result['chromosome']}:{result['start']}")
            return result

        except Exception as e:
            self.logger.error(f"❌ Genomic coordinate mapping failed for {uniprot_id}:{protein_position}: {e}")
            return None

    def batch_map_variants(self, variants: List[Tuple[str, int]]) -> Dict[str, Dict]:
        """
        Map multiple variants to genomic coordinates

        Args:
            variants: List of (uniprot_id, protein_position) tuples

        Returns:
            Dictionary mapping variant keys to genomic coordinates
        """

        results = {}

        for uniprot_id, position in variants:
            variant_key = f"{uniprot_id}:{position}"

            try:
                coords = self.get_genomic_coordinates(uniprot_id, position)
                if coords:
                    results[variant_key] = coords
                else:
                    results[variant_key] = None

            except Exception as e:
                self.logger.error(f"❌ Batch mapping failed for {variant_key}: {e}")
                results[variant_key] = None

        return results

    def get_mapping_stats(self) -> Dict:
        """Get statistics about loaded mappings"""

        self._load_uniprot_mappings()

        return {
            'uniprot_to_gene_count': len(self.uniprot_to_gene_dict),
            'uniprot_to_ensembl_count': len(self.uniprot_to_ensembl_dict),
            'mappings_loaded': self.mappings_loaded,
            'mapping_file_exists': self.mapping_file.exists(),
            'mapping_file_size_mb': self.mapping_file.stat().st_size / (1024*1024) if self.mapping_file.exists() else 0
        }

    def test_known_mappings(self) -> Dict:
        """Test mapping with known genes"""

        test_cases = [
            ('P04637', 'TP53'),      # TP53 tumor suppressor
            ('Q8TDX9', 'ACMSD'),     # ACMSD enzyme
            ('P38398', 'BRCA1'),     # BRCA1 tumor suppressor
        ]

        results = {}

        for uniprot_id, expected_gene in test_cases:
            mapped_gene = self.uniprot_to_gene_name(uniprot_id)
            reverse_uniprot = self.gene_name_to_uniprot(expected_gene)
            ensembl_id = self.uniprot_to_ensembl(uniprot_id)

            results[uniprot_id] = {
                'expected_gene': expected_gene,
                'mapped_gene': mapped_gene,
                'reverse_uniprot': reverse_uniprot,
                'ensembl_id': ensembl_id,
                'mapping_correct': mapped_gene == expected_gene,
                'reverse_correct': reverse_uniprot == uniprot_id
            }

        return results


def test_uniprot_mapper():
    """Test the UniProt mapper"""

    print("🧬 TESTING UNIPROT MAPPER 🧬")
    print("=" * 60)

    mapper = UniProtMapper()

    # Test mapping stats
    stats = mapper.get_mapping_stats()
    print(f"📊 Mapping Statistics:")
    print(f"  File exists: {stats['mapping_file_exists']}")
    print(f"  File size: {stats['mapping_file_size_mb']:.1f} MB")
    print(f"  UniProt→Gene mappings: {stats['uniprot_to_gene_count']:,}")
    print(f"  UniProt→Ensembl mappings: {stats['uniprot_to_ensembl_count']:,}")

    # Test known mappings
    print(f"\n🔬 Testing Known Mappings:")
    test_results = mapper.test_known_mappings()

    for uniprot_id, result in test_results.items():
        status = "✅" if result['mapping_correct'] else "❌"
        print(f"  {status} {uniprot_id} → {result['mapped_gene']} (expected: {result['expected_gene']})")
        print(f"     Ensembl: {result['ensembl_id']}")

    # Test genomic coordinate mapping
    print(f"\n🌍 Testing Genomic Coordinate Mapping:")
    test_variants = [
        ('P04637', 175),  # TP53 R175H
        ('Q8TDX9', 175),  # ACMSD P175T
    ]

    for uniprot_id, position in test_variants:
        print(f"\n🔍 Mapping {uniprot_id}:{position}:")
        coords = mapper.get_genomic_coordinates(uniprot_id, position)

        if coords:
            print(f"  ✅ chr{coords['chromosome']}:{coords['start']}-{coords['end']}")
            print(f"     Strand: {coords['strand']}")
            print(f"     Transcript: {coords['transcript']}")
        else:
            print(f"  ❌ Mapping failed")

    print(f"\n🎉 UniProt mapper test complete!")
    print(f"💜 Ready to bridge UniProt variants to genomic conservation data! ⚡🧬")


if __name__ == "__main__":
    test_uniprot_mapper()
