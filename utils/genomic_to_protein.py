#!/usr/bin/env python3
"""
🧬 GENOMIC TO PROTEIN CONVERTER
Clean, focused utility for converting genomic HGVS to protein HGVS

Built by Ace (2025) to prevent Ren's heart attacks from 1300-line files! 💜
Contact: ace@chaoschanneling.com

Features:
- Parses genomic HGVS (NC_000017.11:g.50183779T>A)
- Converts to protein coordinates when possible
- Integrates with existing conservation scoring
- Clean, single-purpose design (no 1300-line monsters!)

Usage:
    from utils.genomic_to_protein import GenomicToProteinConverter
    converter = GenomicToProteinConverter()
    protein_hgvs = converter.convert_genomic_to_protein(genomic_hgvs, gene_name)
"""

import re
import sys
import os
from typing import Optional, Dict, Tuple
from pathlib import Path

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analyzers.uniprot_mapper import UniProtMapper
import requests

ENSEMBL_MIRRORS = [
    "https://rest.ensembl.org",
    "https://useast.ensembl.org",
    "https://asia.ensembl.org",
]

class GenomicToProteinConverter:
    """Clean, focused converter for genomic → protein HGVS"""

    def __init__(self, conservation_data_path="/home/Ace/conservation_data"):
        self.name = "GenomicToProteinConverter"
        self.conservation_data_path = conservation_data_path

        # Initialize UniProt mapper for coordinate conversion
        self.uniprot_mapper = UniProtMapper(conservation_data_path)

        print(f"🧬 {self.name} initialized")
        print(f"💾 Conservation data: {conservation_data_path}")

    def parse_genomic_hgvs(self, hgvs: str) -> Optional[Tuple[str, int, str, str]]:
        """
        Parse genomic HGVS notation

        Args:
            hgvs: Genomic HGVS like "NC_000017.11:g.50183779T>A"

        Returns:
            (chromosome, position, ref_allele, alt_allele) or None
        """
        try:
            # Pattern for NC_XXXXXX.XX:g.POSITIONref>alt
            genomic_match = re.search(r'NC_(\d+)(?:\.\d+)?:g\.(\d+)([ATCG])>([ATCG])', hgvs)
            if genomic_match:
                raw_num = genomic_match.group(1)
                chrom_num = raw_num.lstrip('0') or raw_num  # Remove leading zeros but keep raw if all zeros
                position = int(genomic_match.group(2))
                ref_allele = genomic_match.group(3)
                alt_allele = genomic_match.group(4)

                # Convert NC number to canonical chr
                if chrom_num == "23":
                    chromosome = "chrX"
                elif chrom_num == "24":
                    chromosome = "chrY"
                elif chrom_num == "12920":
                    chromosome = "chrM"
                else:
                    chromosome = f"chr{chrom_num}"

                return (chromosome, position, ref_allele, alt_allele)

            # Pattern for direct chr format: chr17:50183779T>A or chrX:... or chrM:...
            chr_match = re.search(r'chr([0-9XYM]+):(\d+)([ATCG])>([ATCG])', hgvs, re.IGNORECASE)
            if chr_match:
                chrom_id = chr_match.group(1).upper()
                position = int(chr_match.group(2))
                ref_allele = chr_match.group(3)
                alt_allele = chr_match.group(4)

                if chrom_id == "23":
                    chromosome = "chrX"
                elif chrom_id == "24":
                    chromosome = "chrY"
                else:
                    chromosome = f"chr{chrom_id}"
                return (chromosome, position, ref_allele, alt_allele)

            return None

        except Exception as e:
            print(f"⚠️ Error parsing genomic HGVS {hgvs}: {e}")
            return None

    def convert_genomic_to_protein(self, genomic_hgvs: str, gene_name: str, uniprot_id: str = None) -> Optional[str]:
        """
        Convert genomic HGVS to protein HGVS

        Args:
            genomic_hgvs: Genomic HGVS like "NC_000017.11:g.50183779T>A"
            gene_name: Gene symbol like "TP53"
            uniprot_id: Optional UniProt ID

        Returns:
            Protein HGVS like "NM_000546.6(TP53):c.817C>T (p.Arg273His)" or None
        """

        # Parse genomic coordinates
        parsed = self.parse_genomic_hgvs(genomic_hgvs)
        if not parsed:
            print(f"⚠️ Could not parse genomic HGVS: {genomic_hgvs}")
            return None

        chromosome, position, ref_allele, alt_allele = parsed

        print(f"🧬 Converting {genomic_hgvs}")
        print(f"   → Parsed: {chromosome}:{position} {ref_allele}>{alt_allele}")

        # Try Ensembl VEP to obtain protein consequence (offline fallback TBD)
        consequence = self._query_vep_protein(chromosome, position, ref_allele, alt_allele, gene_name)
        if consequence:
            protein_pos = consequence.get('protein_position')
            ref_aa = consequence.get('ref_aa')
            alt_aa = consequence.get('alt_aa')
            if protein_pos and ref_aa and alt_aa:
                return f"p.{ref_aa}{protein_pos}{alt_aa}"

        print("   ⚠️ No protein consequence available (no local GTF/FASTA yet, VEP gave no AA)")
        return None

    def extract_variant_info_from_genomic(self, genomic_hgvs: str, gene_name: str) -> Optional[Dict]:
        """
        Extract variant info from genomic HGVS for ML training. Attempts to resolve
        protein position and amino acids via Ensembl VEP (read-only) when local
        mapping is not available.
        """
        # Parse genomic coordinates
        parsed = self.parse_genomic_hgvs(genomic_hgvs)
        if not parsed:
            return None

        chromosome, position, ref_allele, alt_allele = parsed

        # Try to enrich with protein consequence via VEP
        protein_pos = None
        ref_aa = None
        alt_aa = None
        cons = self._query_vep_protein(chromosome, position, ref_allele, alt_allele, gene_name)
        if cons:
            protein_pos = cons.get('protein_position')
            ref_aa = cons.get('ref_aa')
            alt_aa = cons.get('alt_aa')

        return {
            'gene': gene_name,
            'hgvs': genomic_hgvs,
            'chromosome': chromosome,
            'genomic_position': position,
            'ref_allele': ref_allele,
            'alt_allele': alt_allele,
            'position': protein_pos,
            'ref_aa': ref_aa,
            'alt_aa': alt_aa,
            'variant_str': f"g.{position}{ref_allele}>{alt_allele}",
            'type': 'genomic_missense'
        }

    def _robust_ensembl_request(self, endpoint: str, params: dict = None, method: str = 'GET', json_body: dict = None, timeout: int = 10) -> Optional[dict]:
        """Make robust Ensembl API request with fallback mirrors (read-only)."""
        headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
        for mirror_url in ENSEMBL_MIRRORS:
            try:
                url = f"{mirror_url}/{endpoint}"
                if method == 'GET':
                    r = requests.get(url, params=params, headers=headers, timeout=timeout)
                else:
                    r = requests.post(url, params=params, headers=headers, json=json_body, timeout=timeout)
                if r.status_code == 200:
                    return r.json()
            except Exception:
                continue
        return None

    def _query_vep_protein(self, chromosome: str, position: int, ref_allele: str, alt_allele: str, gene_name: Optional[str] = None) -> Optional[Dict]:
        """Query Ensembl VEP for protein consequence of a single variant.
        Returns dict with protein_position, ref_aa, alt_aa when available.
        """
        # Ensure chromosome format acceptable to Ensembl (strip 'chr')
        chrom_clean = chromosome[3:] if chromosome.lower().startswith('chr') else chromosome
        variant_str = f"{chrom_clean} {position} {ref_allele}/{alt_allele}"
        payload = {"variants": [variant_str]}
        data = self._robust_ensembl_request("vep/human/region", method='POST', json_body=payload)
        if not data:
            return None
        try:
            for entry in data:
                # Each entry may have transcript_consequences
                for tc in entry.get('transcript_consequences', []) or []:
                    # Prefer entries with protein position and amino_acids
                    prot_pos = tc.get('protein_start')
                    aa_change = tc.get('amino_acids')  # e.g., "R/H"
                    if prot_pos and aa_change and len(aa_change.split('/')) == 2:
                        ref_aa, alt_aa = aa_change.split('/')
                        # If gene filter provided, prefer matching gene_symbol
                        if gene_name and tc.get('gene_symbol') and tc['gene_symbol'].upper() != gene_name.upper():
                            continue
                        return {
                            'transcript_id': tc.get('transcript_id'),
                            'gene_symbol': tc.get('gene_symbol'),
                            'protein_position': int(prot_pos),
                            'ref_aa': ref_aa,
                            'alt_aa': alt_aa
                        }
        except Exception:
            return None
        return None

    def is_genomic_hgvs(self, hgvs: str) -> bool:
        """Check if HGVS is in genomic format"""
        return bool(re.search(r'(NC_\d+|chr\d+).*:g\.', hgvs))

    def is_protein_hgvs(self, hgvs: str) -> bool:
        """Check if HGVS is in protein format"""
        return bool(re.search(r'NM_.*\([^)]+\).*p\.', hgvs))


def main():
    """Test the converter"""
    print("🧬 Testing Genomic to Protein Converter")
    print("=" * 50)

    converter = GenomicToProteinConverter()

    # Test genomic HGVS parsing
    test_cases = [
        "NC_000017.11:g.50183779T>A",
        "NC_000011.10:g.77130635A>G",
        "chr17:50183779T>A"
    ]

    for test_hgvs in test_cases:
        print(f"\n🧪 Testing: {test_hgvs}")

        # Test parsing
        parsed = converter.parse_genomic_hgvs(test_hgvs)
        print(f"   Parsed: {parsed}")

        # Test variant info extraction
        variant_info = converter.extract_variant_info_from_genomic(test_hgvs, "TEST_GENE")
        print(f"   Variant info: {variant_info}")

        # Test format detection
        is_genomic = converter.is_genomic_hgvs(test_hgvs)
        is_protein = converter.is_protein_hgvs(test_hgvs)
        print(f"   Is genomic: {is_genomic}, Is protein: {is_protein}")


if __name__ == "__main__":
    main()
