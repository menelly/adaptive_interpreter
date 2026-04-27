#!/usr/bin/env python3
"""
Post-hoc clean v3 cascade outputs:
- Drop stop-gained / stop-lost / stop-altered variants (anything with * in p.X)
- Drop rows with empty `variant` field
- Keep only rows where variant looks like a real missense substitution

Writes filtered files to outputs_missense_v3_clean/ — leaves originals untouched.
This is for analysis-only; the proper fix is re-cascading with cleaned inputs.
"""
import os
import csv

SRC_DIR = '/home/Ace/AdaptiveInterpreter/outputs_missense_v3'
DST_DIR = '/home/Ace/AdaptiveInterpreter/outputs_missense_v3_clean'

os.makedirs(DST_DIR, exist_ok=True)


def is_real_missense(variant: str) -> bool:
    """Return True only for proper p.RefPosAlt missense (no stops, no asterisks)."""
    if not variant or not variant.startswith('p.'):
        return False
    if '*' in variant:
        return False
    # The form should be p.<letter><digits><letter> roughly
    return True


def clean_one(src_path: str, dst_path: str) -> tuple:
    n_in = n_out = n_stop = n_other = 0
    with open(src_path) as fin, open(dst_path, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        header = next(reader)
        writer.writerow(header)
        try:
            var_idx = header.index('variant')
        except ValueError:
            return 0, 0, 0, 0
        for row in reader:
            n_in += 1
            if len(row) <= var_idx:
                n_other += 1
                continue
            variant = row[var_idx]
            if '*' in (variant or ''):
                n_stop += 1
                continue
            if not is_real_missense(variant):
                n_other += 1
                continue
            writer.writerow(row)
            n_out += 1
    return n_in, n_out, n_stop, n_other


def main():
    files = sorted(f for f in os.listdir(SRC_DIR) if f.endswith('_cascade.tsv'))
    print(f"Cleaning {len(files)} cascade outputs → {DST_DIR}\n")
    print(f"{'Gene':<12} {'In':>7} {'Out':>7} {'Stop':>6} {'Other':>6}")
    print('-' * 50)
    grand_in = grand_out = grand_stop = grand_other = 0
    for fname in files:
        gene = fname.replace('_cascade.tsv', '')
        n_in, n_out, n_stop, n_other = clean_one(
            os.path.join(SRC_DIR, fname),
            os.path.join(DST_DIR, fname),
        )
        grand_in += n_in
        grand_out += n_out
        grand_stop += n_stop
        grand_other += n_other
        if n_stop > 0 or n_other > 0:
            print(f"{gene:<12} {n_in:>7} {n_out:>7} {n_stop:>6} {n_other:>6}")
    print('-' * 50)
    print(f"{'TOTAL':<12} {grand_in:>7} {grand_out:>7} {grand_stop:>6} {grand_other:>6}")


if __name__ == '__main__':
    main()
