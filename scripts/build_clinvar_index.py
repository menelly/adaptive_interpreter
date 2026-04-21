#!/usr/bin/env python3
"""Build a gene → ClinVar disease-submission index from a local ClinVar VCF.

Produces a compact JSON mapping each gene symbol to:
  - total submission count
  - P / LP / VUS / B submission counts
  - top 5 associated diseases by frequency

This index powers utils/clinical_enricher.py's gene-classification tier.

Usage:
    python3 scripts/build_clinvar_index.py \\
        --vcf /mnt/arcana/clinvar/clinvar.vcf.gz \\
        --output /home/Ace/clinvar_gene_index.json
"""
from __future__ import annotations
import argparse
import gzip
import json
import re
import time
from collections import defaultdict, Counter
from pathlib import Path

DEFAULT_VCF = '/mnt/arcana/clinvar/clinvar.vcf.gz'
DEFAULT_OUT = '/home/Ace/clinvar_gene_index.json'

RE_GENE = re.compile(r'GENEINFO=([^;]+)')
RE_SIG  = re.compile(r'CLNSIG=([^;]+)')
RE_DN   = re.compile(r'CLNDN=([^;]+)')


def build(vcf_path: Path, out_path: Path, verbose: bool = True) -> int:
    """Parse ClinVar VCF, aggregate per-gene, write JSON. Returns gene count."""
    t0 = time.time()
    gene_data: Dict[str, Dict] = defaultdict(lambda: {
        'total': 0, 'pathogenic': 0, 'likely_pathogenic': 0,
        'vus': 0, 'benign': 0, 'diseases': Counter(),
    })

    n = 0
    with gzip.open(vcf_path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#'):
                continue
            n += 1
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            info = parts[7]
            m_gene = RE_GENE.search(info)
            if not m_gene:
                continue
            m_sig = RE_SIG.search(info)
            m_dn  = RE_DN.search(info)
            sig = (m_sig.group(1) if m_sig else '').lower()
            dn_raw = m_dn.group(1) if m_dn else ''

            is_p  = 'pathogenic' in sig and 'likely' not in sig and 'conflicting' not in sig
            is_lp = 'likely_pathogenic' in sig
            is_vus = 'uncertain_significance' in sig
            is_b  = 'benign' in sig and 'likely' not in sig
            diseases = [
                d.replace('_', ' ').strip()
                for d in dn_raw.split('|')
                if d and d not in ('not_provided', 'not_specified')
            ]

            genes = [g.split(':')[0] for g in m_gene.group(1).split('|')]
            for g in genes:
                if not g:
                    continue
                rec = gene_data[g]
                rec['total'] += 1
                if is_p:     rec['pathogenic'] += 1
                elif is_lp:  rec['likely_pathogenic'] += 1
                elif is_vus: rec['vus'] += 1
                elif is_b:   rec['benign'] += 1
                for dis in diseases:
                    rec['diseases'][dis] += 1
            if verbose and n % 500_000 == 0:
                print(f'  {n:,} variants processed, {len(gene_data):,} genes so far')

    out: Dict[str, Dict] = {}
    for g, rec in gene_data.items():
        out[g] = {
            'total_submissions':  rec['total'],
            'P_submissions':       rec['pathogenic'],
            'LP_submissions':      rec['likely_pathogenic'],
            'VUS_submissions':     rec['vus'],
            'B_submissions':       rec['benign'],
            'top_diseases': [
                (dis, n_) for dis, n_ in rec['diseases'].most_common(5)
                if dis
            ],
        }

    with open(out_path, 'w') as f:
        json.dump(out, f, separators=(',', ':'))

    if verbose:
        print(f'✅ {len(out):,} genes indexed in {time.time()-t0:.1f}s → {out_path}')
        print(f'   File size: {out_path.stat().st_size / 1024 / 1024:.2f} MB')
    return len(out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--vcf', default=DEFAULT_VCF,
                    help=f'Path to ClinVar VCF (default {DEFAULT_VCF})')
    ap.add_argument('--output', default=DEFAULT_OUT,
                    help=f'Output JSON path (default {DEFAULT_OUT})')
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args()
    build(Path(args.vcf), Path(args.output), verbose=not args.quiet)
