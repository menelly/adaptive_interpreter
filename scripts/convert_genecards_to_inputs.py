#!/usr/bin/env python3
"""
Convert GeneCards TSV exports to AdaptiveInterpreter cascade input format.

Input  (GeneCards paste, may or may not have header):
  Accession | rsID | Clinical sig | Chrpos | Variation Name | Ref/Alt | AA Chg | Type

Output (cascade_batch_processor expects):
  gene | hgvs_p | hgvs_g | gnomad_freq | clinical_sig | molecular_consequence | review_status

cascade_batch_processor accepts hgvs_p directly (see cascade_batch_processor.py:496),
so we don't need to convert to genomic coordinates — the AA change is sufficient.
"""
import os
import re
import sys
import glob

GENECARDS_HEADERS = ['Accession', 'rsID', 'Clinical significance and condition',
                     'Chrpos', 'Variation Name', 'Ref/Alt', 'AA Chg', 'Type']

INPUT_GLOB = '/home/Ace/analysis/Untitled spreadsheet - *.tsv'
OUT_DIR = '/home/Ace/AdaptiveInterpreter/inputs_missense_only'


def normalize_clinical_sig(sig_text: str) -> str:
    """Reduce verbose ClinVar string to canonical bucket."""
    s = sig_text.lower()
    if 'pathogenic' in s and 'likely pathogenic' in s:
        return 'Pathogenic/Likely_pathogenic'
    if 'pathogenic' in s and 'likely' not in s.split('pathogenic')[0][-15:]:
        return 'Pathogenic'
    if 'likely pathogenic' in s:
        return 'Likely_pathogenic'
    if 'likely benign' in s:
        return 'Likely_benign'
    if 'benign' in s:
        return 'Benign'
    if 'uncertain' in s:
        return 'Uncertain_significance'
    if 'conflict' in s:
        return 'Conflicting_classifications_of_pathogenicity'
    return 'not_provided'


def aa_to_hgvs_p(aa_chg: str) -> str:
    """Convert 'A168D' to 'p.A168D'."""
    s = (aa_chg or '').strip()
    if not s:
        return ''
    if s.startswith('p.'):
        return s
    return f'p.{s}'


def convert_one(infile: str, gene: str, outfile: str) -> tuple:
    """Convert one GeneCards file. Returns (n_total, n_missense, n_unique_written)."""
    # Detect header
    with open(infile) as f:
        first = f.readline().strip()
    has_header = first.lower().startswith('accession\trsid')

    rows_out = []
    seen = set()
    n_total = 0
    n_missense = 0

    out_header = ['gene', 'hgvs_p', 'hgvs_g', 'gnomad_freq', 'clinical_sig',
                  'molecular_consequence', 'review_status']

    with open(infile) as f:
        if has_header:
            f.readline()  # skip header
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 8:
                continue
            n_total += 1
            vtype = parts[7].strip().lower()
            if 'missense' not in vtype:
                continue
            n_missense += 1
            aa_chg = parts[6].strip()
            hgvs_p = aa_to_hgvs_p(aa_chg)
            if not hgvs_p:
                continue
            if hgvs_p in seen:
                continue
            seen.add(hgvs_p)
            clinical_sig = normalize_clinical_sig(parts[2])
            row = [gene, hgvs_p, '', '0.0', clinical_sig, 'missense', 'criteria_provided,_multiple_submitters']
            rows_out.append('\t'.join(row))

    with open(outfile, 'w') as f:
        f.write('\t'.join(out_header) + '\n')
        f.write('\n'.join(rows_out) + '\n')

    return n_total, n_missense, len(seen)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    files = sorted(glob.glob(INPUT_GLOB))
    if not files:
        print(f"No files matching {INPUT_GLOB}")
        sys.exit(1)

    print(f"{'Gene':<10} {'Total':>7} {'Missense':>9} {'Unique':>7}  -> output")
    print('-' * 65)
    for fp in files:
        name = os.path.basename(fp)
        # Extract gene from "Untitled spreadsheet - GENE.tsv"
        m = re.search(r'-\s*([A-Z0-9]+)\.tsv$', name)
        if not m:
            print(f"  SKIP (can't parse gene): {name}")
            continue
        gene = m.group(1)
        outfile = os.path.join(OUT_DIR, f'{gene}_input.tsv')
        n_total, n_miss, n_uniq = convert_one(fp, gene, outfile)
        print(f"{gene:<10} {n_total:>7} {n_miss:>9} {n_uniq:>7}  -> {outfile}")


if __name__ == '__main__':
    main()
