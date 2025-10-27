#!/usr/bin/env python3
"""
🔥💜 CLINVAR BULK EXTRACTOR - NO MORE DIALUP SOUNDS! 🚀

Extract variants for specific genes from ClinVar bulk downloads.
Nova's brilliant solution to avoid ClinVar Miner's 20-second page loads!

Usage:
    python clinvar_bulk_extractor.py --genes COL1A1,COL3A1,SCN5A --output extracted_variants/
"""

import argparse
import gzip
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set
import pandas as pd
from collections import defaultdict

class ClinVarBulkExtractor:
    """🧬 Extract gene-specific variants from ClinVar bulk data"""

    def __init__(self, clinvar_dir: str = "/mnt/Arcana/clinvar"):
        self.clinvar_dir = Path(clinvar_dir)
        self.vcf_file = self.clinvar_dir / "clinvar.vcf.gz"

        # Significance mappings
        self.significance_map = {
            'Pathogenic': 'pathogenic',
            'Likely_pathogenic': 'pathogenic',
            'Benign': 'benign',
            'Likely_benign': 'benign',
            'Uncertain_significance': 'vus',
            'Conflicting_interpretations_of_pathogenicity': 'conflicting',
            'drug_response': 'drug_response',
            'association': 'association',
            'risk_factor': 'risk_factor',
            'protective': 'protective',
            'other': 'other'
        }

        print(f"🔍 ClinVar directory: {self.clinvar_dir}")
        print(f"🧬 VCF file: {self.vcf_file}")

    def extract_gene_variants(self, target_genes: Set[str]) -> Dict[str, List[Dict]]:
        """Extract variants for target genes from ClinVar VCF"""
        print(f"🎯 Extracting variants for {len(target_genes)} genes: {', '.join(sorted(target_genes))}")

        if not self.vcf_file.exists():
            raise FileNotFoundError(f"ClinVar VCF not found: {self.vcf_file}")

        gene_variants = defaultdict(list)
        total_variants = 0
        extracted_variants = 0

        print("📊 Processing ClinVar VCF...")

        with gzip.open(self.vcf_file, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                if line_num % 100000 == 0:
                    print(f"   📋 Processed {line_num:,} lines, extracted {extracted_variants:,} variants")

                # Skip header lines
                if line.startswith('#'):
                    continue

                total_variants += 1

                # Parse VCF line
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue

                chrom, pos, variant_id, ref, alt, qual, filter_status, info = fields[:8]

                # Extract gene symbols from INFO field
                gene_symbols = self._extract_gene_symbols(info)

                # Check if any target genes are present
                matching_genes = gene_symbols.intersection(target_genes)
                if not matching_genes:
                    continue

                # Extract variant information
                variant_info = self._parse_variant_info(chrom, pos, variant_id, ref, alt, info)
                if not variant_info:
                    continue

                # Add to results for each matching gene
                for gene in matching_genes:
                    variant_info['gene'] = gene
                    gene_variants[gene].append(variant_info.copy())
                    extracted_variants += 1

        print(f"✅ Extraction complete!")
        print(f"   📊 Total variants processed: {total_variants:,}")
        print(f"   🎯 Variants extracted: {extracted_variants:,}")
        print(f"   🧬 Genes with variants: {len(gene_variants)}")

        return dict(gene_variants)

    def _extract_gene_symbols(self, info_field: str) -> Set[str]:
        """Extract gene symbols from VCF INFO field"""
        gene_symbols = set()

        # Look for GENEINFO field: GENEINFO=GENE1:123|GENE2:456
        geneinfo_match = re.search(r'GENEINFO=([^;]+)', info_field)
        if geneinfo_match:
            geneinfo = geneinfo_match.group(1)
            # Split by | and extract gene names before :
            for gene_entry in geneinfo.split('|'):
                if ':' in gene_entry:
                    gene_name = gene_entry.split(':')[0]
                    if gene_name and gene_name != '.':
                        gene_symbols.add(gene_name)

        return gene_symbols

    def _parse_variant_info(self, chrom: str, pos: str, variant_id: str,
                          ref: str, alt: str, info: str) -> Optional[Dict]:
        """Parse variant information from VCF fields"""
        try:
            # Extract clinical significance
            clnsig_match = re.search(r'CLNSIG=([^;]+)', info)
            if not clnsig_match:
                return None

            clnsig_raw = clnsig_match.group(1)
            significance = self._map_significance(clnsig_raw)

            # Extract HGVS notation
            hgvs_list = self._extract_hgvs(info)

            # Extract allele frequency if available
            frequency = self._extract_frequency(info)

            # Extract review status
            review_status = self._extract_review_status(info)

            # Extract molecular consequence if present
            molecular_consequence = self._extract_molecular_consequence(info, hgvs_list)

            return {
                'chromosome': chrom,
                'position': int(pos),
                'variant_id': variant_id,
                'ref_allele': ref,
                'alt_allele': alt,
                'significance': significance,
                'significance_raw': clnsig_raw,
                'hgvs_list': hgvs_list,
                'frequency': frequency,
                'review_status': review_status,
                'molecular_consequence': molecular_consequence,
            }

        except Exception as e:
            print(f"⚠️ Error parsing variant {variant_id}: {e}")
            return None

    def _map_significance(self, clnsig_raw: str) -> str:
        """Map ClinVar significance to our categories"""
        # Handle multiple significances separated by |
        significances = clnsig_raw.split('|')

        # Priority order: pathogenic > benign > vus > others
        for priority_sig in ['Pathogenic', 'Likely_pathogenic', 'Benign', 'Likely_benign']:
            if any(priority_sig in sig for sig in significances):
                return self.significance_map.get(priority_sig, 'other')

        # Default to first significance if no priority match
        first_sig = significances[0] if significances else 'other'
        return self.significance_map.get(first_sig, 'other')

    def _extract_hgvs(self, info: str) -> List[str]:
        """Extract HGVS notations from INFO field"""
        hgvs_list = []

        # Look for CLNHGVS field
        hgvs_match = re.search(r'CLNHGVS=([^;]+)', info)
        if hgvs_match:
            hgvs_raw = hgvs_match.group(1)
            # Split by | and clean up
            for hgvs in hgvs_raw.split('|'):
                if hgvs and hgvs != '.':
                    hgvs_list.append(hgvs)

        return hgvs_list

    def _extract_frequency(self, info: str) -> Optional[float]:
        """Extract allele frequency if available"""
        # Look for AF field (allele frequency)
        af_match = re.search(r'AF=([^;]+)', info)
        if af_match:
            try:
                return float(af_match.group(1))
            except ValueError:
                pass

        return None

    def _extract_review_status(self, info: str) -> str:
        """Extract review status"""
        # Look for CLNREVSTAT field
        revstat_match = re.search(r'CLNREVSTAT=([^;]+)', info)
        if revstat_match:
            return revstat_match.group(1)
        return 'unknown'

    def _extract_molecular_consequence(self, info: str, hgvs_list: List[str]) -> str:
        """Extract and normalize molecular consequence.
        Prefer MC= from INFO (Sequence Ontology), fallback to heuristics on HGVS.
        Returns a simple tag like: missense, synonymous, stop_gained, frameshift, splice_acceptor,
        splice_donor, intronic, utr, other, unknown.
        """
        try:
            # Try MC field: e.g., MC=SO:0001583|missense_variant,SO:0001587|stop_gained
            m = re.search(r'MC=([^;]+)', info)
            if m:
                raw = m.group(1)
                # MC can be comma-separated entries
                terms = []
                for tok in raw.split(','):
                    parts = tok.split('|')
                    if len(parts) == 2:
                        terms.append(parts[1].strip().lower())
                    else:
                        terms.append(parts[0].strip().lower())
                # Preference order
                pref = [
                    'stop_gained', 'stop_lost', 'frameshift_variant',
                    'splice_acceptor_variant', 'splice_donor_variant', 'missense_variant',
                    'synonymous_variant', 'intron_variant', '5_prime_utr_variant', '3_prime_utr_variant'
                ]
                for p in pref:
                    if any(p in t for t in terms):
                        return self._normalize_mc_label(p)
                # Fallback to first
                if terms:
                    return self._normalize_mc_label(terms[0])
            # Heuristic fallbacks using HGVS list
            p_hgvs = next((h for h in hgvs_list if 'p.' in h), None)
            if p_hgvs:
                # Synonymous: p.Pro123Pro or p.Ala359=
                if re.search(r"p\.([A-Z][a-z]{2}|[A-Z])\d+(=|\1)", p_hgvs):
                    return 'synonymous'
                # Stop gained: p.Trp26* or p.W26*
                if re.search(r"p\.[A-Z][a-z]{2}\d+\*", p_hgvs) or re.search(r"p\.[A-Z]\d+\*", p_hgvs):
                    return 'stop_gained'
                # Frameshift: contains fs
                if 'fs' in p_hgvs:
                    return 'frameshift'
                # Default for protein-level without stop/fs markers
                return 'missense'
            # No protein HGVS; if only genomic present, likely intronic/UTR/unknown
            # Try CLNVC (variant class) for hints
            vcm = re.search(r'CLNVC=([^;]+)', info)
            if vcm:
                vclass = vcm.group(1).lower()
                if 'single_nucleotide' in vclass and hgvs_list and any('c.' in h for h in hgvs_list):
                    # coding SNV without protein annotation — unknown coding
                    return 'other'
            return 'unknown'
        except Exception:
            return 'unknown'

    def _normalize_mc_label(self, term: str) -> str:
        t = term.lower()
        if 'missense' in t: return 'missense'
        if 'synonymous' in t: return 'synonymous'
        if 'stop_gained' in t: return 'stop_gained'
        if 'stop_lost' in t: return 'stop_lost'
        if 'frameshift' in t: return 'frameshift'
        if 'splice_acceptor' in t: return 'splice_acceptor'
        if 'splice_donor' in t: return 'splice_donor'
        if 'intron' in t: return 'intronic'
        if 'utr' in t: return 'utr'
        return 'other'

    def save_gene_variants(self, gene_variants: Dict[str, List[Dict]],
                          output_dir: str, format: str = 'tsv'):
        """Save extracted variants to files"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print(f"💾 Saving variants to {output_path}")

        for gene, variants in gene_variants.items():
            if not variants:
                continue

            print(f"   📝 {gene}: {len(variants)} variants")

            # Group by significance
            by_significance = defaultdict(list)
            for variant in variants:
                by_significance[variant['significance']].append(variant)

            # Save each significance category
            for significance, sig_variants in by_significance.items():
                if significance in ['pathogenic', 'benign']:  # Only save P/B for ML training
                    filename = f"{gene.lower()}_{significance}.{format}"
                    filepath = output_path / filename

                    if format == 'tsv':
                        self._save_tsv(sig_variants, filepath)
                    elif format == 'json':
                        self._save_json(sig_variants, filepath)

                    print(f"      ✅ {significance}: {len(sig_variants)} variants → {filename}")

        print("🎉 All variants saved!")

    def _save_tsv(self, variants: List[Dict], filepath: Path):
        """Save variants as TSV in our standard format"""
        rows = []

        for variant in variants:
            # Find the best HGVS notation (prefer protein changes)
            best_hgvs = self._select_best_hgvs(variant['hgvs_list'])

            if best_hgvs:
                rows.append({
                    'HGVS': best_hgvs,
                    'dbSNP': variant['variant_id'] if variant['variant_id'].startswith('rs') else '',
                    'gnomAD frequency': variant['frequency'] or 0.0,
                    'ClinVar_significance': variant['significance_raw'],
                    'Review_status': variant['review_status']
                })

        if rows:
            df = pd.DataFrame(rows)
            df.to_csv(filepath, sep='\t', index=False)

    def _save_json(self, variants: List[Dict], filepath: Path):
        """Save variants as JSON"""
        with open(filepath, 'w') as f:
            json.dump(variants, f, indent=2)

    def _select_best_hgvs(self, hgvs_list: List[str]) -> Optional[str]:
        """Select the best HGVS notation for ML training"""
        if not hgvs_list:
            return None

        # Prefer protein changes (p.) over cDNA (c.) over genomic (g.)
        protein_hgvs = [h for h in hgvs_list if 'p.' in h]
        if protein_hgvs:
            return protein_hgvs[0]

        cdna_hgvs = [h for h in hgvs_list if 'c.' in h]
        if cdna_hgvs:
            return cdna_hgvs[0]

        # Return first available
        return hgvs_list[0]


def main():
    parser = argparse.ArgumentParser(description="Extract gene variants from ClinVar bulk data")
    parser.add_argument('--genes', required=True,
                       help='Comma-separated list of gene symbols (e.g., COL1A1,COL3A1,SCN5A)')
    parser.add_argument('--output', required=True,
                       help='Output directory for extracted variants')
    parser.add_argument('--format', choices=['tsv', 'json'], default='tsv',
                       help='Output format (default: tsv)')
    parser.add_argument('--clinvar-dir', default='/mnt/Arcana/clinvar',
                       help='Directory containing ClinVar files')

    args = parser.parse_args()

    # Parse target genes
    target_genes = set(gene.strip().upper() for gene in args.genes.split(','))

    print("🔥💜 CLINVAR BULK EXTRACTOR - NOVA'S BRILLIANT SOLUTION! 🚀")
    print("=" * 60)

    try:
        # Initialize extractor
        extractor = ClinVarBulkExtractor(args.clinvar_dir)

        # Extract variants
        gene_variants = extractor.extract_gene_variants(target_genes)

        # Save results
        extractor.save_gene_variants(gene_variants, args.output, args.format)

        print("\n🎉 SUCCESS! No more dialup sounds! 📞💀")

    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
