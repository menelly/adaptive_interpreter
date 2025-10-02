#!/usr/bin/env python3
"""
Offline transcript/CDS mapper using local GTF + FASTA to map
(gene_name, protein_position) -> genomic position (codon center) for GRCh38.

- Parses GENCODE-style GTF (expects gene_name/transcript_id/exon_number attributes)
- Uses pyfaidx for random access to GRCh38 FASTA (supports .fa or .fa.gz)
- Strand-aware, CDS-ordered concatenation
- Returns a minimal dict suitable for conservation lookups

Notes:
- We choose the transcript with the longest CDS that still covers the requested
  amino acid position (>= 3*aa_pos bases). If multiple qualify, we prefer the
  one with most CDS exons.
- We return the middle base of the codon for conservation (more stable than
  first/third and symmetric across strand handling). If desired, callers can
  extend this to use all three bases.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

from pathlib import Path

# Lazy import so modules that don't need offline won't fail if pyfaidx is missing
try:
    from pyfaidx import Fasta
except Exception:  # pragma: no cover
    Fasta = None  # type: ignore


@dataclass
class CdsExon:
    chrom: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    strand: str  # '+' or '-'
    exon_number: Optional[int]


class OfflineProteinToGenomicMapper:
    def __init__(self, gtf_path: str, fasta_path: str):
        self.gtf_path = Path(gtf_path)
        self.fasta_path = Path(fasta_path)
        # pyfaidx is optional for coordinate mapping (needed only to access sequence)
        # If unavailable, we still perform coordinate mapping using GTF alone.
        if not self.gtf_path.exists():
            raise FileNotFoundError(f"GTF not found: {self.gtf_path}")
        if not self.fasta_path.exists():
            # Allow running without FASTA if only coordinates are needed
            self._fasta = None
        else:
            self._fasta = None  # opened lazily when/if sequence is requested
        # Caches
        self._gene_transcripts_cache: Dict[str, Dict[str, List[CdsExon]]] = {}

    def _open_fasta(self):
        if self._fasta is None:
            # as_raw=True avoids sequence validation overhead; we don't need it
            self._fasta = Fasta(str(self.fasta_path), as_raw=True)

    @staticmethod
    def _parse_attributes(attr_field: str) -> Dict[str, str]:
        out = {}
        # GENCODE style: key "value"; pairs separated by ;
        for part in attr_field.strip().split(';'):
            part = part.strip()
            if not part:
                continue
            if '"' in part:
                key, val = part.split(' ', 1)
                val = val.strip().strip('"')
                out[key] = val
            else:
                # key value (no quotes)
                toks = part.split()
                if len(toks) >= 2:
                    out[toks[0]] = toks[1]
        return out

    def _iter_gtf(self):
        opener = gzip.open if self.gtf_path.suffix.endswith('gz') else open
        with opener(self.gtf_path, 'rt', encoding='utf-8', newline='') as fh:
            for line in fh:
                if not line or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                chrom, source, feature, start, end, score, strand, frame, attrs = parts
                if feature != 'CDS':
                    continue
                attr = self._parse_attributes(attrs)
                gene_name = attr.get('gene_name') or attr.get('gene_id')
                transcript_id = attr.get('transcript_id')
                exon_number = None
                try:
                    exon_number = int(attr.get('exon_number')) if attr.get('exon_number') else None
                except Exception:
                    exon_number = None
                try:
                    yield gene_name, transcript_id, CdsExon(
                        chrom=chrom,
                        start=int(start),
                        end=int(end),
                        strand=strand,
                        exon_number=exon_number,
                    )
                except Exception:
                    # Skip malformed rows
                    continue

    @lru_cache(maxsize=1024)
    def _get_transcripts_for_gene(self, gene_name: str) -> Dict[str, List[CdsExon]]:
        tx_map: Dict[str, List[CdsExon]] = {}
        for gname, tx_id, cds in self._iter_gtf():
            if gname != gene_name:
                continue
            if not tx_id:
                continue
            tx_map.setdefault(tx_id, []).append(cds)
        # Normalize exon ordering within each transcript
        for tx_id, exons in tx_map.items():
            if not exons:
                continue
            strand = exons[0].strand
            if strand == '+':
                # Sort by genomic start; break ties by exon_number if present
                exons.sort(key=lambda e: (e.start, e.end, e.exon_number or 0))
            else:
                # Negative strand: CDS order goes from high to low genomic coordinates
                exons.sort(key=lambda e: (-(e.end), -(e.start), e.exon_number or 0))
        return tx_map

    def _build_cds_map(self, exons: List[CdsExon]) -> Tuple[List[Tuple[str, int]], str]:
        """
        Build a linear map of CDS genomic positions in coding order.
        Returns (positions, strand) where positions is a list of (chrom, pos) 1-based.
        """
        positions: List[Tuple[str, int]] = []
        if not exons:
            return positions, '+'
        strand = exons[0].strand
        for exon in exons:
            chrom = exon.chrom
            if strand == '+':
                for p in range(exon.start, exon.end + 1):
                    positions.append((chrom, p))
            else:
                # For negative strand, coding order is from exon.end down to exon.start
                for p in range(exon.end, exon.start - 1, -1):
                    positions.append((chrom, p))
        return positions, strand

    def map_gene_pos(self, gene_name: str, aa_pos: int) -> Optional[Dict]:
        """
        Map gene_name + amino acid position (1-based) to genomic position (codon center).
        Returns dict with chromosome (as in GTF, e.g., 'chr17'), genomic_position (1-based),
        strand ('+'/'-'), transcript_id, and codon_positions (list of 3 genomic positions).
        """
        if aa_pos <= 0:
            return None
        tx_map = self._get_transcripts_for_gene(gene_name)
        if not tx_map:
            return None
        # Choose transcript that covers the AA position; prefer longest CDS
        candidates = []
        for tx_id, exons in tx_map.items():
            positions, strand = self._build_cds_map(exons)
            if len(positions) >= aa_pos * 3:
                candidates.append((len(positions), len(exons), tx_id, positions, strand))
        if not candidates:
            return None
        candidates.sort(reverse=True)  # longest CDS first, then most exons
        _, _, tx_id, positions, strand = candidates[0]
        idx0 = (aa_pos - 1) * 3
        codon_pos = positions[idx0: idx0 + 3]
        if len(codon_pos) < 3:
            return None
        # Choose middle base of codon for conservation
        chrom_mid, genomic_mid = codon_pos[1]
        return {
            'chromosome': chrom_mid,
            'genomic_position': genomic_mid,
            'strand': 1 if strand == '+' else -1,
            'transcript_id': tx_id,
            'codon_positions': codon_pos,
        }

