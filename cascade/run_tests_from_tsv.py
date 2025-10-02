#!/usr/bin/env python3
import csv
import json
import re
import sys
from pathlib import Path

from DNModeling.cascade.cascade_analyzer import CascadeAnalyzer

HEADER_HINTS = [
    'Clinical significance and condition', 'Chrpos', 'Variation Name', 'AA Chg'
]

# Extract gene symbol from string like: NM_000222.3(KIT):c.2459A>G (p.Asp820Gly)
GENE_RE = re.compile(r"\(([^)]+)\)")
PROT_RE = re.compile(r"p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}")
AA3_TO_1 = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Glu':'E','Gln':'Q','Gly':'G','His':'H',
    'Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'
}

AA_CHANGE_RE = re.compile(r"([A-Z])(\d+)([A-Z])$")
PROT_HGVS_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")


def parse_gene_and_protein(variation_name: str, aa_chg: str):
    gene = None
    prot = None

    # Gene
    m = GENE_RE.search(variation_name)
    if m:
        gene = m.group(1).strip()

    # Protein change
    # Prefer explicit p.Xxx###Yyy in variation_name
    m2 = PROT_HGVS_RE.search(variation_name)
    if m2:
        ref3, pos, alt3 = m2.groups()
        ref = AA3_TO_1.get(ref3, None)
        alt = AA3_TO_1.get(alt3, None)
        if ref and alt:
            prot = f"p.{ref}{pos}{alt}"
    # Else fall back to AA Chg like D820G
    if not prot and aa_chg:
        m3 = AA_CHANGE_RE.match(aa_chg.strip())
        if m3:
            ref, pos, alt = m3.groups()
            prot = f"p.{ref}{pos}{alt}"

    return gene, prot


def run_file(path: Path):
    ca = CascadeAnalyzer()
    results = []

    with path.open('r', encoding='utf-8') as f:
        # Sniff delimiter: TSV expected, but some files might be comma within quotes
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None)
        if header is None:
            print(f"Empty file: {path}")
            return []

        # Basic sanity: if the header looks comma-separated, re-read as CSV
        if len(header) == 1 and any(h in header[0] for h in HEADER_HINTS):
            f.seek(0)
            reader = csv.reader(f)
            header = next(reader)

        # Map columns
        cols = {name: idx for idx, name in enumerate(header)}
        var_idx = cols.get('Variation Name', None)
        aa_idx = cols.get('AA Chg', None)
        if var_idx is None:
            print(f"Cannot find 'Variation Name' column in {path}")
            return []

        for row in reader:
            if not row or len(row) <= var_idx:
                continue
            variation_name = row[var_idx]
            aa_chg = row[aa_idx] if aa_idx is not None and aa_idx < len(row) else ''

            gene, prot = parse_gene_and_protein(variation_name, aa_chg)
            if not gene or not prot:
                continue

            try:
                res = ca.analyze_cascade_biological(gene, prot, gnomad_freq=0.0, variant_type='missense', sequence=None)
                results.append({
                    'gene': gene,
                    'variant': prot,
                    'family': res.get('gene_family'),
                    'family_aa_multiplier_applied': res.get('family_aa_multiplier_applied', 1.0),
                    'final_score': res.get('final_score'),
                    'final_classification': res.get('final_classification'),
                })
            except Exception as e:
                results.append({'gene': gene, 'variant': prot, 'error': str(e)})

    return results


def main():
    if len(sys.argv) < 2:
        print("Usage: run_tests_from_tsv.py <path.tsv> [<path2.tsv> ...]")
        sys.exit(1)

    all_results = []
    for p in sys.argv[1:]:
        r = run_file(Path(p))
        all_results.extend(r)

    print(json.dumps(all_results, indent=2))


if __name__ == '__main__':
    main()

