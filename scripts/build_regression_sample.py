#!/usr/bin/env python3
"""
Build a frozen, stratified regression sample for validating scoring changes.

Pulls 50 Pathogenic + 50 VUS + 50 Benign missense variants (with protein notation)
from inputs_missense_only/, seeded for reproducibility. Freeze the output, run it
through the cascade to capture a baseline, then re-run after a scoring change and diff
— so we can prove a change (e.g. the GOF interpretation gate) didn't break the rest.

Usage:  python scripts/build_regression_sample.py [N_PER_BUCKET] [SEED]
Writes:  validation/regression_sample_<3N>.tsv  (gene, hgvs_p, expected_class, clinical_sig)
"""
import csv
import glob
import os
import random
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 50
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 42


def bucket_of(sig: str) -> str | None:
    s = sig.lower()
    if "conflict" in s:
        return None
    if "pathogenic" in s and "benign" not in s:
        return "P"
    if "benign" in s and "pathogenic" not in s:
        return "B"
    if "uncertain" in s:
        return "VUS"
    return None


def main():
    random.seed(SEED)
    rows = {"P": [], "VUS": [], "B": []}
    for f in sorted(glob.glob("inputs_missense_only/*.tsv")):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                hp = (row.get("hgvs_p") or "").strip()
                mc = (row.get("molecular_consequence") or "").lower()
                rs = (row.get("review_status") or "").lower()
                if "no_assertion" in rs or not hp or "missense" not in mc:
                    continue
                b = bucket_of(row.get("clinical_sig") or "")
                if b:
                    rows[b].append((row["gene"], hp, row.get("clinical_sig", "")))

    os.makedirs("validation", exist_ok=True)
    out = f"validation/regression_sample_{3 * N}.tsv"
    counts = {}
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "hgvs_p", "expected_class", "clinical_sig"])
        for b in ("P", "VUS", "B"):
            samp = random.sample(rows[b], min(N, len(rows[b])))
            counts[b] = len(samp)
            for gene, hp, sig in sorted(samp):
                w.writerow([gene, hp, b, sig])
    print(f"wrote {out}: " + ", ".join(f"{b}={counts[b]}" for b in ("P", "VUS", "B")))
    print(f"(seed={SEED}; reproducible — re-running gives the same 150)")


if __name__ == "__main__":
    main()
