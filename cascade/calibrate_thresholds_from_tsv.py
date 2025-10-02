#!/usr/bin/env python3
import csv
import io
import json
import re
import sys
from pathlib import Path
from contextlib import redirect_stdout

from DNModeling.cascade.cascade_analyzer import CascadeAnalyzer

# Heuristics similar to run_tsv_and_export
GENE_IN_PARENS = re.compile(r"\(([^)]+)\)")
PROT_HGVS_3LETTER = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")
PROT_HGVS_1LETTER = re.compile(r"p\.([A-Z])(\d+)([A-Z])")
AA_COMPACT = re.compile(r"^([A-Z])(\d+)([A-Z])$")
AA3_TO_1 = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Glu':'E','Gln':'Q','Gly':'G','His':'H',
    'Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'
}

CLIN_KEYS = (
    'Clinical significance and condition',
    'Clinical significance',
)


def parse_gene_and_protein_and_label(row, header):
    # Try detect clinical significance column
    clin_idx = None
    for k in CLIN_KEYS:
        if header and k in header:
            clin_idx = header.index(k)
            break
    # Heuristic when no header: find a column with Benign/Pathogenic keywords
    if clin_idx is None:
        for i, tok in enumerate(row):
            t = tok or ""
            if any(x in t for x in ("Pathogenic", "Benign", "Likely benign", "Likely pathogenic")):
                clin_idx = i
                break

    clin_text = (row[clin_idx] if clin_idx is not None and clin_idx < len(row) else "") or ""

    # Parse gene + protein
    # First try "Variation Name" if header has it
    var_idx = header.index('Variation Name') if header and 'Variation Name' in header else None
    if var_idx is not None:
        var_name = row[var_idx]
        m = GENE_IN_PARENS.search(var_name)
        gene = m.group(1).strip() if m else None
        prot = None
        m1 = PROT_HGVS_1LETTER.search(var_name)
        if m1:
            ref, pos, alt = m1.groups(); prot = f"p.{ref}{pos}{alt}"
        else:
            m2 = PROT_HGVS_3LETTER.search(var_name)
            if m2:
                ref3, pos, alt3 = m2.groups(); ref = AA3_TO_1.get(ref3); alt = AA3_TO_1.get(alt3)
                if ref and alt: prot = f"p.{ref}{pos}{alt}"
    else:
        # Headerless / different schema: heuristics
        gene, prot = None, None
        # Often the NM_...(...):... column is present; find it by pattern NM_
        nm_idx = None
        for i, tok in enumerate(row):
            if tok and 'NM_' in tok and '(' in tok and ')' in tok:
                nm_idx = i; break
        if nm_idx is not None:
            var_name = row[nm_idx]
            m = GENE_IN_PARENS.search(var_name)
            gene = m.group(1).strip() if m else None
            m1 = PROT_HGVS_1LETTER.search(var_name)
            if m1:
                ref, pos, alt = m1.groups(); prot = f"p.{ref}{pos}{alt}"
            else:
                m2 = PROT_HGVS_3LETTER.search(var_name)
                if m2:
                    ref3, pos, alt3 = m2.groups(); ref = AA3_TO_1.get(ref3); alt = AA3_TO_1.get(alt3)
                    if ref and alt: prot = f"p.{ref}{pos}{alt}"
        # If still no prot, try compact AA field
        if not prot:
            for tok in row:
                m3 = AA_COMPACT.match((tok or '').strip())
                if m3:
                    r, p, a = m3.groups(); prot = f"p.{r}{p}{a}"; break
        # Try grab gene from another column (rsID files may not have explicit gene column)
        # If not found, we will skip this row

    if not gene or not prot:
        return None, None, None

    # Label extraction: path if contains Pathogenic (without Benign), benign if contains Benign (without Pathogenic)
    low = clin_text.lower()
    is_path = (('pathogenic' in low) and ('benign' not in low))
    is_benign = (('benign' in low) and ('pathogenic' not in low))
    label = 'P' if is_path else ('B' if is_benign else None)

    return gene, prot, label


def calibrate(files, k=3, out_path=Path('DNModeling/cascade/resources/classification_thresholds.json')):
    ca = CascadeAnalyzer()
    fam_to_scores = {}

    for fp in files:
        p = Path(fp)
        with p.open('r', encoding='utf-8', newline='') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader, None)
            # If header looks wrong (single cell), try comma
            if header and len(header) == 1:
                f.seek(0)
                reader = csv.reader(f)
                header = next(reader, None)
            # If there is no header and the first row is data, rewind
            if header and any(s for s in header if 'NM_' in s or s.startswith('rs')):
                # Probably not a real header; treat it as data
                f.seek(0)
                reader = csv.reader(f, delimiter='\t')
                header = None

            for row in reader:
                if not row:
                    continue
                gene, prot, label = parse_gene_and_protein_and_label(row, header)
                if not gene or not prot or not label:
                    continue
                # Suppress verbose prints during calibration
                buf = io.StringIO()
                with redirect_stdout(buf):
                    res = ca.analyze_cascade_biological(gene, prot, gnomad_freq=0.0, variant_type='missense', sequence=None)
                fam = res.get('gene_family') or 'GENERAL'
                score = res.get('final_score')
                if score is None:
                    continue
                fam_to_scores.setdefault(fam, {'B': [], 'P': []})[label].append(float(score))

    # Compute thresholds per family
    thresholds = {}
    for fam, groups in fam_to_scores.items():
        ben = sorted(groups['B'], reverse=True)
        pat = sorted(groups['P'])
        if not ben or not pat:
            continue
        kk = max(1, min(k, len(ben), len(pat)))
        benign_upper = max(ben[:kk])
        path_lower = max(pat[:kk])  # the highest among the lowest kk pathogenic scores
        # Set LP just above hardest benign; P at or just above weak pathogenic cluster
        eps = 0.01
        lp = benign_upper + eps
        pthr = max(lp + eps, path_lower + eps)
        # Store only LP/P overrides; other thresholds use defaults
        thresholds[fam] = {"LP": round(lp, 3), "P": round(pthr, 3)}

    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Merge with existing (if present)
    merged = {}
    if out_path.exists():
        try:
            with out_path.open('r', encoding='utf-8') as fh:
                merged = json.load(fh) or {}
        except Exception:
            merged = {}
    merged.update(thresholds)
    with out_path.open('w', encoding='utf-8') as fh:
        json.dump(merged, fh, indent=2)
    return merged


def main():
    if len(sys.argv) < 2:
        print("Usage: calibrate_thresholds_from_tsv.py <tsv_or_csv> [more files...]")
        sys.exit(1)
    thresholds = calibrate(sys.argv[1:])
    print(json.dumps(thresholds, indent=2))


if __name__ == '__main__':
    main()

