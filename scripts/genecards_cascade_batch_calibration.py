#!/usr/bin/env python3
"""
Cascade analysis on multiple GeneCards TSV exports for CumBurSum calibration.
Loops through 4 AR genes (ATP7B, CFTR, DMD, POLG) and writes one cascade TSV per gene.

Uses CascadeAnalyzer class directly (the top-level analyze_variant_cascade function
was refactored away in earlier versions of the codebase).
"""
import pandas as pd
import sys
import os
import re
from pathlib import Path

sys.path.insert(0, '/home/Ace')
from AdaptiveInterpreter.analyzers.cascade_analyzer import CascadeAnalyzer

GENECARDS_HEADERS = ['Accession', 'rsID', 'Clinical significance and condition',
                     'Chrpos', 'Variation Name', 'Ref/Alt', 'AA Chg', 'Type']

INPUT_FILES = {
    'ATP7B': '/home/Ace/analysis/Untitled spreadsheet - ATP7B.tsv',
    'CFTR':  '/home/Ace/analysis/Untitled spreadsheet - CFTR.tsv',
    'DMD':   '/home/Ace/analysis/Untitled spreadsheet - DMD.tsv',
    'POLG':  '/home/Ace/analysis/Untitled spreadsheet - POLG.tsv',
}

OUTPUT_DIR = '/home/Ace/analysis/calibration_extension'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Single global analyzer instance — initialization is expensive
_ANALYZER = None
def get_analyzer():
    global _ANALYZER
    if _ANALYZER is None:
        _ANALYZER = CascadeAnalyzer()
    return _ANALYZER


def load_genecards_tsv(path: str) -> pd.DataFrame:
    """Load GeneCards TSV — handle files with or without header row."""
    with open(path) as f:
        first = f.readline().strip()
    has_header = first.lower().startswith('accession\trsid')

    if has_header:
        df = pd.read_csv(path, sep='\t')
    else:
        df = pd.read_csv(path, sep='\t', header=None, names=GENECARDS_HEADERS)
    return df


def aa_chg_to_hgvs_p(aa_chg: str) -> str:
    """Convert 'A168D' style to 'p.A168D' if not already prefixed."""
    s = str(aa_chg).strip()
    if not s:
        return ''
    if s.startswith('p.'):
        return s
    return f'p.{s}'


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


def run_gene(gene: str, input_path: str) -> tuple:
    """Run cascade on all missense variants in a GeneCards file."""
    print(f"\n{'='*70}")
    print(f"Processing {gene}")
    print(f"{'='*70}")

    df = load_genecards_tsv(input_path)
    print(f"Loaded {len(df)} rows from {os.path.basename(input_path)}")

    if 'Type' in df.columns:
        before = len(df)
        df = df[df['Type'].astype(str).str.contains('missense', case=False, na=False)]
        print(f"Filtered to {len(df)} missense (dropped {before - len(df)} non-missense)")

    analyzer = get_analyzer()
    results = []
    failed = 0
    seen_variants = set()  # dedupe — GeneCards often has duplicates

    for idx, row in df.iterrows():
        aa_chg = row.get('AA Chg', '')
        hgvs_p = aa_chg_to_hgvs_p(aa_chg)
        if not hgvs_p:
            continue
        if hgvs_p in seen_variants:
            continue
        seen_variants.add(hgvs_p)

        clin_text = str(row.get('Clinical significance and condition', ''))
        clinical_sig = normalize_clinical_sig(clin_text)

        try:
            result = analyzer.analyze_cascade(
                gene=gene,
                variant=hgvs_p,
                gnomad_freq=0.0,
                variant_type='missense',
                expected_clinvar=clinical_sig,
            )
            # Flatten any nested dicts to scalar columns
            flat = {}
            for k, v in result.items():
                if isinstance(v, (dict, list)):
                    flat[k] = str(v)[:500]
                else:
                    flat[k] = v
            flat['gene'] = gene
            flat['variant'] = hgvs_p
            flat['clinvar_sig_text'] = clin_text[:200]
            flat['clinical_sig_normalized'] = clinical_sig
            flat['rsid'] = str(row.get('rsID', ''))
            flat['variation_name'] = str(row.get('Variation Name', ''))
            results.append(flat)
        except Exception as e:
            failed += 1
            if failed <= 3:
                print(f"  Error on {hgvs_p}: {e}")

        if (len(results) + failed) % 200 == 0 and (len(results) + failed) > 0:
            print(f"  Progress: success={len(results)}, failed={failed}, dedup-skipped={len(df) - len(seen_variants)}")

    print(f"\n{gene} complete: {len(results)} success, {failed} failed, {len(seen_variants)} unique variants tried")

    if results:
        out_path = f"{OUTPUT_DIR}/{gene}_cascade.tsv"
        results_df = pd.DataFrame(results)
        results_df.to_csv(out_path, sep='\t', index=False)
        print(f"  Saved: {out_path}")

    return len(df), len(results), failed


def main():
    print(f"GeneCards cascade calibration extension")
    print(f"Output: {OUTPUT_DIR}")

    summary = []
    for gene, path in INPUT_FILES.items():
        if not os.path.exists(path):
            print(f"\n{gene}: SKIPPED (file not found)")
            continue
        n_total, n_success, n_failed = run_gene(gene, path)
        summary.append((gene, n_total, n_success, n_failed))

    print(f"\n{'='*70}")
    print(f"FINAL SUMMARY")
    print(f"{'='*70}")
    print(f"{'Gene':<8} {'Input':>7} {'Success':>9} {'Failed':>7}")
    for gene, total, success, failed in summary:
        print(f"{gene:<8} {total:>7} {success:>9} {failed:>7}")


if __name__ == '__main__':
    main()
