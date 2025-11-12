#!/usr/bin/env python3
"""
🧬 OFFLINE GENOMIC TO PROTEIN CONVERTER
Convert genomic coordinates to protein changes using local GTF and FASTA files

Built by Ace (2025) - No more flaky APIs! 💜
Uses local GENCODE GTF and reference genome for accurate, fast conversions

Features:
- Parses GENCODE GTF for gene/exon coordinates
- Uses reference genome FASTA for DNA sequences
- Converts genomic variants to protein changes
- Caches gene models for speed
- No API dependencies!

Usage:
    from utils.offline_genomic_to_protein import OfflineGenomicToProteinConverter
    converter = OfflineGenomicToProteinConverter()
    protein_change = converter.convert("chr19", 46755451, "A", "C", "FKRP")
"""

import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import pickle

# Genetic code for translation
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


class OfflineGenomicToProteinConverter:
    """Convert genomic coordinates to protein changes using local files"""
    
    def __init__(self, 
                 gtf_path: str = "/mnt/Arcana/genetics_data/reference/gencode.v46.annotation.gtf.gz",
                 fasta_path: str = "/mnt/Arcana/genetics_data/reference/GRCh38.primary_assembly.genome.fa.gz",
                 cache_dir: str = "/home/Ace/AdaptiveInterpreter/.cache"):
        
        self.gtf_path = Path(gtf_path)
        self.fasta_path = Path(fasta_path)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Gene models cache
        self.gene_models = {}
        self.gene_cache_file = self.cache_dir / "gene_models.pkl"
        
        # Chromosome sequences cache (loaded on demand)
        self.chr_sequences = {}
        
        print(f"🧬 Offline Genomic-to-Protein Converter initialized")
        print(f"   GTF: {self.gtf_path}")
        print(f"   FASTA: {self.fasta_path}")
        print(f"   Cache: {self.cache_dir}")
        
        # Load or build gene models
        self._load_gene_models()
    
    def _load_gene_models(self):
        """Load gene models from cache or build from GTF"""
        if self.gene_cache_file.exists():
            print("📦 Loading cached gene models...")
            with open(self.gene_cache_file, 'rb') as f:
                self.gene_models = pickle.load(f)
            print(f"   ✅ Loaded {len(self.gene_models)} genes from cache")
        else:
            print("🔨 Building gene models from GTF (this will take a few minutes)...")
            self._build_gene_models()
            print(f"   ✅ Built models for {len(self.gene_models)} genes")
            print("💾 Saving to cache...")
            with open(self.gene_cache_file, 'wb') as f:
                pickle.dump(self.gene_models, f)
            print("   ✅ Cache saved!")
    
    def _build_gene_models(self):
        """Parse GTF and build gene models"""
        gene_transcripts = defaultdict(list)
        
        print("   📖 Parsing GTF file...")
        line_count = 0
        
        with gzip.open(self.gtf_path, 'rt') as f:
            for line in f:
                line_count += 1
                if line_count % 100000 == 0:
                    print(f"      Processed {line_count:,} lines...")
                
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                feature_type = fields[2]
                if feature_type != 'CDS':  # Only care about coding sequences
                    continue
                
                # Parse attributes
                attrs = {}
                for attr in fields[8].split(';'):
                    attr = attr.strip()
                    if not attr:
                        continue
                    match = re.match(r'(\w+)\s+"([^"]+)"', attr)
                    if match:
                        attrs[match.group(1)] = match.group(2)
                
                gene_name = attrs.get('gene_name')
                transcript_id = attrs.get('transcript_id')
                
                if not gene_name or not transcript_id:
                    continue
                
                # Store CDS exon
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                
                gene_transcripts[gene_name].append({
                    'transcript_id': transcript_id,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand
                })
        
        print(f"   📊 Found {len(gene_transcripts)} genes with CDS annotations")
        
        # Build gene models (use longest transcript for each gene)
        for gene_name, cds_list in gene_transcripts.items():
            # Group by transcript
            transcripts = defaultdict(list)
            for cds in cds_list:
                transcripts[cds['transcript_id']].append(cds)
            
            # Find longest transcript
            longest_transcript = None
            longest_length = 0
            
            for transcript_id, exons in transcripts.items():
                total_length = sum(e['end'] - e['start'] + 1 for e in exons)
                if total_length > longest_length:
                    longest_length = total_length
                    longest_transcript = exons
            
            if longest_transcript:
                # Sort exons by position
                longest_transcript.sort(key=lambda x: x['start'])
                
                self.gene_models[gene_name] = {
                    'chrom': longest_transcript[0]['chrom'],
                    'strand': longest_transcript[0]['strand'],
                    'exons': [(e['start'], e['end']) for e in longest_transcript]
                }
    
    def _load_chromosome_sequence(self, chrom: str) -> Optional[str]:
        """Load chromosome sequence from FASTA (cached)"""
        if chrom in self.chr_sequences:
            return self.chr_sequences[chrom]
        
        print(f"   📖 Loading {chrom} sequence from FASTA...")
        
        # Normalize chromosome name
        chrom_normalized = chrom if chrom.startswith('chr') else f'chr{chrom}'
        
        sequence = []
        in_target_chr = False
        
        with gzip.open(self.fasta_path, 'rt') as f:
            for line in f:
                if line.startswith('>'):
                    # Check if this is our target chromosome
                    header = line[1:].split()[0]
                    if header == chrom_normalized or header == chrom:
                        in_target_chr = True
                        print(f"      Found {chrom_normalized}!")
                    else:
                        if in_target_chr:
                            # We've finished reading our chromosome
                            break
                        in_target_chr = False
                elif in_target_chr:
                    sequence.append(line.strip().upper())
        
        if sequence:
            full_sequence = ''.join(sequence)
            self.chr_sequences[chrom] = full_sequence
            print(f"      ✅ Loaded {len(full_sequence):,} bp")
            return full_sequence
        
        print(f"      ⚠️ Chromosome {chrom} not found in FASTA")
        return None
    
    def convert(self, chrom: str, position: int, ref_allele: str, alt_allele: str, gene_name: str) -> Optional[str]:
        """
        Convert genomic variant to protein change
        
        Args:
            chrom: Chromosome (e.g., "chr19" or "19")
            position: Genomic position (1-based)
            ref_allele: Reference allele (e.g., "A")
            alt_allele: Alternate allele (e.g., "C")
            gene_name: Gene symbol (e.g., "FKRP")
        
        Returns:
            Protein change (e.g., "p.L276I") or None if not in coding region
        """
        # Normalize chromosome name
        chrom = chrom if chrom.startswith('chr') else f'chr{chrom}'
        
        # Check if gene is in our models
        if gene_name not in self.gene_models:
            print(f"   ⚠️ Gene {gene_name} not found in gene models")
            return None
        
        gene_model = self.gene_models[gene_name]
        
        # Check if variant is in coding region
        cds_position = self._genomic_to_cds_position(position, gene_model)
        if cds_position is None:
            print(f"   ⚠️ Position {position} not in coding region of {gene_name}")
            return None
        
        # Load chromosome sequence
        chr_seq = self._load_chromosome_sequence(chrom)
        if not chr_seq:
            return None
        
        # Get reference codon and mutated codon
        aa_position = (cds_position - 1) // 3 + 1
        codon_start_cds = ((aa_position - 1) * 3) + 1
        
        # Get the codon sequence from genomic coordinates
        ref_codon, genomic_codon_positions = self._get_codon_sequence(codon_start_cds, gene_model, chr_seq)
        if not ref_codon:
            return None
        
        # Apply mutation
        alt_codon = self._apply_mutation(ref_codon, genomic_codon_positions, position, ref_allele, alt_allele)
        if not alt_codon:
            return None
        
        # Translate codons
        ref_aa = CODON_TABLE.get(ref_codon, '?')
        alt_aa = CODON_TABLE.get(alt_codon, '?')
        
        if ref_aa == '?' or alt_aa == '?':
            print(f"   ⚠️ Unknown codon: {ref_codon} or {alt_codon}")
            return None
        
        return f"p.{ref_aa}{aa_position}{alt_aa}"
    
    def _genomic_to_cds_position(self, genomic_pos: int, gene_model: Dict) -> Optional[int]:
        """Convert genomic position to CDS position"""
        exons = gene_model['exons']
        strand = gene_model['strand']
        
        # Check if position is in any exon
        cds_pos = 0
        for start, end in exons:
            if start <= genomic_pos <= end:
                # Position is in this exon
                if strand == '+':
                    cds_pos += (genomic_pos - start + 1)
                else:
                    cds_pos += (end - genomic_pos + 1)
                return cds_pos
            else:
                # Add this exon's length to CDS position
                if strand == '+':
                    if genomic_pos > end:
                        cds_pos += (end - start + 1)
                else:
                    if genomic_pos < start:
                        cds_pos += (end - start + 1)
        
        return None  # Not in coding region
    
    def _get_codon_sequence(self, cds_position: int, gene_model: Dict, chr_seq: str) -> Tuple[Optional[str], Optional[List[int]]]:
        """Get codon sequence and genomic positions for a CDS position

        Returns:
            (codon_sequence, [genomic_pos1, genomic_pos2, genomic_pos3])
        """
        exons = gene_model['exons']
        strand = gene_model['strand']

        # Get the 3 nucleotides for this codon
        codon_cds_positions = [cds_position, cds_position + 1, cds_position + 2]
        codon_genomic_positions = []
        codon_nucleotides = []

        for cds_pos in codon_cds_positions:
            genomic_pos = self._cds_to_genomic_position(cds_pos, gene_model)
            if genomic_pos is None:
                return None, None

            codon_genomic_positions.append(genomic_pos)

            # Get nucleotide from reference genome (1-based to 0-based)
            try:
                nucleotide = chr_seq[genomic_pos - 1]
                codon_nucleotides.append(nucleotide)
            except IndexError:
                print(f"   ⚠️ Position {genomic_pos} out of range for chromosome")
                return None, None

        codon_seq = ''.join(codon_nucleotides)

        # Reverse complement if on minus strand
        if strand == '-':
            codon_seq = self._reverse_complement(codon_seq)

        return codon_seq, codon_genomic_positions

    def _cds_to_genomic_position(self, cds_pos: int, gene_model: Dict) -> Optional[int]:
        """Convert CDS position to genomic position"""
        exons = gene_model['exons']
        strand = gene_model['strand']

        current_cds_pos = 0

        if strand == '+':
            # Forward strand: walk through exons left to right
            for start, end in exons:
                exon_length = end - start + 1
                if current_cds_pos + exon_length >= cds_pos:
                    # Position is in this exon
                    offset = cds_pos - current_cds_pos - 1
                    return start + offset
                current_cds_pos += exon_length
        else:
            # Reverse strand: walk through exons right to left
            for start, end in reversed(exons):
                exon_length = end - start + 1
                if current_cds_pos + exon_length >= cds_pos:
                    # Position is in this exon
                    offset = cds_pos - current_cds_pos - 1
                    return end - offset
                current_cds_pos += exon_length

        return None

    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    
    def _apply_mutation(self, ref_codon: str, genomic_positions: List[int], 
                       mut_pos: int, ref_allele: str, alt_allele: str) -> Optional[str]:
        """Apply mutation to codon"""
        # Find which position in the codon corresponds to the mutation
        try:
            codon_idx = genomic_positions.index(mut_pos)
            alt_codon = list(ref_codon)
            alt_codon[codon_idx] = alt_allele
            return ''.join(alt_codon)
        except (ValueError, IndexError):
            return None


def convert_tsv_file(input_file: str, output_file: str):
    """
    Convert a TSV file with genomic HGVS to protein HGVS using OFFLINE converter

    Args:
        input_file: Path to input TSV file
        output_file: Path to output TSV file
    """
    import pandas as pd
    import sys

    print(f"🔬 Converting {input_file}")
    print(f"📝 Output: {output_file}")

    # Read input file
    df = pd.read_csv(input_file, sep='\t')
    print(f"📊 Loaded {len(df)} variants")

    # Initialize converter
    converter = OfflineGenomicToProteinConverter()

    # Process each variant
    converted = 0
    skipped = 0
    failed = 0

    for idx, row in df.iterrows():
        if idx % 100 == 0:
            print(f"   Processing {idx}/{len(df)}...", end='\r')

        hgvs_p = row.get('hgvs_p', '')
        hgvs_g = row.get('hgvs_g', '')
        gene = row.get('gene', '')

        # Convert NaN to empty string
        if pd.isna(hgvs_p):
            hgvs_p = ''
        if pd.isna(hgvs_g):
            hgvs_g = ''

        # If already has protein HGVS, skip
        if hgvs_p and str(hgvs_p).strip():
            skipped += 1
            continue

        # Parse genomic HGVS (e.g., NC_000001.11:g.236686164C>T)
        if hgvs_g:
            # Try NC_ format first
            match = re.match(r'NC_(\d+)\.\d+:g\.(\d+)([ACGT])>([ACGT])', str(hgvs_g))
            if match:
                nc_num = int(match.group(1))
                pos = int(match.group(2))
                ref = match.group(3)
                alt = match.group(4)

                # Convert NC_ to chromosome (NC_000001 = chr1, etc.)
                # NC_000023 = chrX, NC_000024 = chrY
                if nc_num == 23:
                    chrom = "chrX"
                elif nc_num == 24:
                    chrom = "chrY"
                else:
                    chrom = f"chr{nc_num}"

                try:
                    result = converter.convert(chrom, pos, ref, alt, str(gene))
                    if result:
                        df.at[idx, 'hgvs_p'] = result
                        converted += 1
                    else:
                        failed += 1
                except Exception as e:
                    failed += 1
            else:
                failed += 1

    print(f"\n✅ Converted {converted}/{len(df)} variants")
    print(f"⏭️  Skipped {skipped} (already had protein HGVS)")
    print(f"❌ Failed {failed}")

    # Write output file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"💾 Saved to {output_file}")


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 3:
        # File conversion mode
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        convert_tsv_file(input_file, output_file)
    else:
        # Test mode
        print("🧬 Testing Offline Genomic-to-Protein Converter")
        print("=" * 60)

        converter = OfflineGenomicToProteinConverter()

        # Test with FKRP variant
        result = converter.convert("chr19", 46755451, "A", "C", "FKRP")
        print(f"\n✅ Result: {result}")

