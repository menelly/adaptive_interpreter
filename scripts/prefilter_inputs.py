#!/usr/bin/env python3
"""
Pre-filter cascade input files for speed and cleanliness.

For each *_input.tsv in inputs_missense_only/:
1. BACKFILL: For rows with empty hgvs_p, look up hgvs_g in the matching
   v2 cascade output. If v2 successfully converted it (has p.X variant),
   populate hgvs_p in the input.
2. DROP: After backfill, drop any remaining rows with empty hgvs_p.
   Those are either non-missense (5'UTR, intron, etc.) or unconvertable
   genomic coords; either way they can't be missense-cascaded.

Result: inputs that go straight to missense analysis without VEP API calls.
Original inputs are backed up to inputs_missense_only_PREFILTER_BACKUP/.
"""
import os
import re
import shutil
import csv

INPUT_DIR = '/home/Ace/AdaptiveInterpreter/inputs_missense_only'
BACKUP_DIR = '/home/Ace/AdaptiveInterpreter/inputs_missense_only_PREFILTER_BACKUP'
V2_DIR = '/home/Ace/AdaptiveInterpreter/outputs_missense_v2'


def build_hgvs_map_from_v2(v2_path: str) -> dict:
    """Read v2 cascade output, return {hgvs_g: variant (p.X)} for missense rows."""
    if not os.path.exists(v2_path):
        return {}
    out = {}
    with open(v2_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            variant = row.get('variant', '').strip()
            hgvs_g = row.get('hgvs', '').strip()
            mc = row.get('molecular_consequence', '').lower()
            # Only use rows where variant is a proper protein change AND it's missense
            if variant.startswith('p.') and hgvs_g and 'missense' in mc:
                out[hgvs_g] = variant
    return out


def prefilter_one(input_path: str, gene: str) -> tuple:
    """Returns (n_in, n_backfilled, n_dropped, n_out)."""
    v2_path = os.path.join(V2_DIR, f'{gene}_cascade.tsv')
    backfill_map = build_hgvs_map_from_v2(v2_path)

    rows_in = []
    with open(input_path) as f:
        header_line = f.readline().rstrip('\n')
        headers = header_line.split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            # pad to header width
            if len(parts) < len(headers):
                parts += [''] * (len(headers) - len(parts))
            rows_in.append(dict(zip(headers, parts)))

    rows_out = []
    n_backfilled = 0
    n_dropped = 0
    n_stop_variant = 0

    for r in rows_in:
        hgvs_p = r.get('hgvs_p', '').strip()
        hgvs_g = r.get('hgvs_g', '').strip()
        mc = r.get('molecular_consequence', '').lower()

        # Only keep missense (the folder is "missense_only" but has non-missense leakage)
        if mc and 'missense' not in mc:
            n_dropped += 1
            continue

        if not hgvs_p and hgvs_g and hgvs_g in backfill_map:
            r['hgvs_p'] = backfill_map[hgvs_g]
            n_backfilled += 1
            hgvs_p = r['hgvs_p']

        if not hgvs_p:
            n_dropped += 1
            continue

        # 🚫 Drop stop-gained (p.X*), stop-lost (p.*X), and stop-altered (p.*X*) variants.
        # These ARE NOT missense substitutions even though some made it into "missense" folders.
        # Examples: p.Y22* (stop-gained = nonsense), p.*16W (stop-lost = readthrough),
        #           p.*6*  (synonymous stop change). All have * in the AA change.
        if '*' in hgvs_p:
            n_stop_variant += 1
            n_dropped += 1
            continue

        rows_out.append(r)

    # Write back
    with open(input_path, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for r in rows_out:
            f.write('\t'.join(r.get(h, '') for h in headers) + '\n')

    return len(rows_in), n_backfilled, n_dropped, len(rows_out), n_stop_variant


def main():
    if not os.path.exists(BACKUP_DIR):
        print(f"Backing up inputs → {BACKUP_DIR}")
        shutil.copytree(INPUT_DIR, BACKUP_DIR)
    else:
        print(f"Backup already exists at {BACKUP_DIR} — skipping backup")

    files = sorted(f for f in os.listdir(INPUT_DIR) if f.endswith('_input.tsv'))
    print(f"\nProcessing {len(files)} input files...\n")
    print(f"{'Gene':<12} {'In':>7} {'Backfilled':>12} {'Dropped':>9} {'StopVar':>8} {'Out':>7}  v2-source")
    print('-' * 80)

    grand_in = grand_back = grand_drop = grand_out = grand_stop = 0
    no_v2 = []
    for fname in files:
        gene = fname.replace('_input.tsv', '')
        path = os.path.join(INPUT_DIR, fname)
        v2_path = os.path.join(V2_DIR, f'{gene}_cascade.tsv')
        v2_status = '✓' if os.path.exists(v2_path) else '✗ no v2'
        if not os.path.exists(v2_path):
            no_v2.append(gene)
        n_in, n_back, n_drop, n_out, n_stop = prefilter_one(path, gene)
        grand_in += n_in
        grand_back += n_back
        grand_drop += n_drop
        grand_out += n_out
        grand_stop += n_stop
        print(f"{gene:<12} {n_in:>7} {n_back:>12} {n_drop:>9} {n_stop:>8} {n_out:>7}  {v2_status}")

    print('-' * 80)
    print(f"{'TOTAL':<12} {grand_in:>7} {grand_back:>12} {grand_drop:>9} {grand_stop:>8} {grand_out:>7}")
    print()
    print(f"Backfill saved {grand_back:,} variants from re-running VEP")
    print(f"Dropped {grand_drop:,} unconvertable/non-missense rows ({grand_stop:,} were stop variants leaking into missense)")
    if no_v2:
        print(f"\n⚠️ {len(no_v2)} genes had no v2 cascade — only their hgvs_p-populated rows were kept:")
        print(f"   {', '.join(no_v2)}")


if __name__ == '__main__':
    main()
