#!/usr/bin/env python3
import csv
import json
import re
import sys
import os
import io
import argparse
from contextlib import redirect_stdout
from pathlib import Path
from datetime import datetime

from DNModeling.cascade.cascade_analyzer import CascadeAnalyzer

# Heuristics for column detection
GENE_CANDIDATES = [
    'Gene', 'Symbol', 'Gene Symbol', 'Gene_Symbol', 'HGNC', 'GeneName', 'Gene Name'
]
PROT_CANDIDATES = [
    'Protein', 'Protein change', 'Protein Change', 'AA Chg', 'AA Change',
    'HGVSp', 'HGVSp_Short', 'Protein Consequence', 'Amino acid change', 'Amino Acid Change'
]
VAR_NAME_CANDIDATES = ['Variation Name', 'Variant', 'Variant Name', 'Name']

GENE_IN_PARENS = re.compile(r"\(([^)]+)\)")  # like NM_000222.3(KIT):...
PROT_HGVS_3LETTER = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")
PROT_HGVS_1LETTER = re.compile(r"p\.([A-Z])(\d+)([A-Z])")
AA_COMPACT = re.compile(r"^([A-Z])(\d+)([A-Z])$")
AA3_TO_1 = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Glu':'E','Gln':'Q','Gly':'G','His':'H',
    'Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'
}


def pick_column(header, candidates):
    idx = None
    name = None
    for cand in candidates:
        if cand in header:
            idx = header.index(cand)
            name = cand
            break
    return idx, name


def parse_gene_and_protein(row, header):
    # Headerless TSV support: detect NM_/rs-style rows
    if header is None:
        gene, prot = None, None
        # Prefer NM_... column if present
        nm_idx = None
        for i, tok in enumerate(row):
            if tok and 'NM_' in tok and '(' in tok and ')' in tok:
                nm_idx = i; break
        if nm_idx is not None:
            var_name = row[nm_idx]
            m = GENE_IN_PARENS.search(var_name)
            if m:
                gene = m.group(1).strip()
            m1 = PROT_HGVS_1LETTER.search(var_name)
            if m1:
                ref, pos, alt = m1.groups(); prot = f"p.{ref}{pos}{alt}"
            else:
                m2 = PROT_HGVS_3LETTER.search(var_name)
                if m2:
                    ref3, pos, alt3 = m2.groups(); ref = AA3_TO_1.get(ref3); alt = AA3_TO_1.get(alt3)
                    if ref and alt: prot = f"p.{ref}{pos}{alt}"
        # If still no prot, try compact AA token like G368V
        if not prot:
            for tok in row:
                m3 = AA_COMPACT.search((tok or '').strip())
                if m3:
                    r, p, a = m3.groups(); prot = f"p.{r}{p}{a}"; break
        return gene, prot

    # Try explicit gene column
    gidx, _ = pick_column(header, GENE_CANDIDATES)
    gene = row[gidx].strip() if gidx is not None and gidx < len(row) else None

    # Try explicit protein column (prefer HGVS form)
    pidx, _ = pick_column(header, PROT_CANDIDATES)
    prot_src = row[pidx].strip() if pidx is not None and pidx < len(row) else ''

    prot = None
    if prot_src:
        m1 = PROT_HGVS_1LETTER.search(prot_src)
        if m1:
            ref, pos, alt = m1.groups()
            prot = f"p.{ref}{pos}{alt}"
        else:
            m2 = PROT_HGVS_3LETTER.search(prot_src)
            if m2:
                ref3, pos, alt3 = m2.groups()
                ref = AA3_TO_1.get(ref3); alt = AA3_TO_1.get(alt3)
                if ref and alt: prot = f"p.{ref}{pos}{alt}"
            else:
                m3 = AA_COMPACT.search(prot_src)
                if m3:
                    ref, pos, alt = m3.groups(); prot = f"p.{ref}{pos}{alt}"

    # Fallback: parse from Variation Name
    if not gene or not prot:
        vidx, _ = pick_column(header, VAR_NAME_CANDIDATES)
        var_name = row[vidx] if vidx is not None and vidx < len(row) else ''
        if not gene:
            m = GENE_IN_PARENS.search(var_name)
            if m:
                gene = m.group(1).strip()
        if not prot:
            m1 = PROT_HGVS_1LETTER.search(var_name)
            if m1:
                ref, pos, alt = m1.groups(); prot = f"p.{ref}{pos}{alt}"
            else:
                m2 = PROT_HGVS_3LETTER.search(var_name)
                if m2:
                    ref3, pos, alt3 = m2.groups()
                    ref = AA3_TO_1.get(ref3); alt = AA3_TO_1.get(alt3)
                    if ref and alt: prot = f"p.{ref}{pos}{alt}"

    return gene, prot


def process_file(path: Path, verbose: bool = False, emit_compact: bool = True):
    ca = CascadeAnalyzer()
    results = []

    # Sniff delimiter: try TSV first, then CSV
    with path.open('r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None)
        if header is None:
            return []
        # If header came as single token, try comma
        if len(header) == 1:
            f.seek(0)
            reader = csv.reader(f)
            header = next(reader, None)
            if header is None:
                return []
        # If header looks like data (contains NM_ or starts with rs or has parentheses), treat as headerless TSV
        if header and any((tok and ('NM_' in tok or tok.startswith('rs') or '(' in tok)) for tok in header):
            f.seek(0)
            reader = csv.reader(f, delimiter='\t')
            header = None

        for row in reader:
            if not row:
                continue
            try:
                gene, prot = parse_gene_and_protein(row, header)
            except Exception:
                gene, prot = None, None
            if not gene or not prot:
                continue
            try:
                if verbose:
                    res = ca.analyze_cascade_biological(gene, prot, gnomad_freq=0.0, variant_type='missense', sequence=None)
                else:
                    # Suppress verbose analyzer prints in compact mode
                    buf = io.StringIO()
                    with redirect_stdout(buf):
                        res = ca.analyze_cascade_biological(gene, prot, gnomad_freq=0.0, variant_type='missense', sequence=None)
                rec = {
                    'gene': gene,
                    'variant': prot,
                    'family': res.get('gene_family'),
                    'family_aa_multiplier_applied': res.get('family_aa_multiplier_applied', 1.0),
                    'final_score': res.get('final_score'),
                    'final_classification': res.get('final_classification'),
                }
                results.append(rec)
                if emit_compact and not verbose:
                    # Emoji-free, one-liner summary per variant
                    fs = rec.get('final_score')
                    fs_txt = f"{fs:.3f}" if isinstance(fs, (int, float)) else ""
                    print(f"{rec['gene']} {rec['variant']} | {rec.get('family') or ''} | {rec.get('family_aa_multiplier_applied',1.0):.3f} | {fs_txt} | {rec.get('final_classification') or ''}")
            except Exception as e:
                results.append({'gene': gene, 'variant': prot, 'error': str(e)})
    return results


def to_markdown(results, title):
    lines = [f"# {title}", "", f"Generated: {datetime.utcnow().isoformat()}Z", "", "| gene | variant | family | family_aa_mult | final_score | final_class |", "|---|---|---:|---:|---:|---|"]
    for r in results:
        if 'error' in r:
            lines.append(f"| {r.get('gene','')} | {r.get('variant','')} |  |  |  | ERROR: {r['error']} |")
        else:
            lines.append(
                f"| {r['gene']} | {r['variant']} | {r.get('family') or ''} | {r.get('family_aa_multiplier_applied',1.0):.3f} | {r.get('final_score') if r.get('final_score') is not None else ''} | {r.get('final_classification') or ''} |"
            )
    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Run CascadeAnalyzer on TSV/CSV and export JSON+Markdown.")
    parser.add_argument("files", nargs="+", help="Input TSV/CSV files")
    parser.add_argument("--verbose", action="store_true", help="Show full emoji-rich analyzer logs")
    parser.add_argument("--no-compact", dest="compact", action="store_false", help="Do not print compact per-variant lines")
    args = parser.parse_args()

    # Env override
    env_level = os.getenv("CASCADE_LOG_LEVEL", "").upper()
    verbose = args.verbose or env_level in ("DEBUG", "VERBOSE")
    emit_compact = args.compact and not verbose  # if verbose, skip compact lines

    out_dir = Path('DNModeling/cascade/reports')
    out_dir.mkdir(parents=True, exist_ok=True)

    for in_path in args.files:
        p = Path(in_path)
        results = process_file(p, verbose=verbose, emit_compact=emit_compact)
        base = p.stem.replace(' ', '_')
        json_path = out_dir / f"{base}_results.json"
        md_path = out_dir / f"{base}_results.md"

        with json_path.open('w', encoding='utf-8') as jf:
            json.dump(results, jf, indent=2)
        with md_path.open('w', encoding='utf-8') as mf:
            mf.write(to_markdown(results, title=base))

        print(f"Wrote {json_path} and {md_path}")


if __name__ == '__main__':
    main()

