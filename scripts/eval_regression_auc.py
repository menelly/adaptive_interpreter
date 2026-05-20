#!/usr/bin/env python3
"""
Evaluate a regression run (P-vs-B discrimination) — the measurement half of the harness.

Run after any scoring change and compare to the baseline to PROVE the change moved the
needle (or didn't break it). Reports:
  - pooled ROC-AUC for final_score + each mechanism (P=1, B=0)
  - gene overlap between P and B sets (low overlap => pooled AUC is confounded by
    cross-family score ranges, i.e. the normalization problem)
  - within-gene P-vs-B mean separation for genes that have both (the cleaner signal)

Usage:  python scripts/eval_regression_auc.py [validation/run_baseline.tsv]
"""
import csv
import statistics as st
import sys
from collections import defaultdict


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "validation/run_baseline.tsv"
    rows = [r for r in csv.DictReader(open(path), delimiter="\t")
            if r["final_classification"] != "ERROR"]

    def f(r, k):
        try:
            return float(r[k])
        except (ValueError, KeyError, TypeError):
            return None

    cols = {"final": "final_score", "DN": "DN", "LOF": "LOF", "GOF": "GOF"}
    data = {k: [] for k in cols}
    y = []
    byg = defaultdict(lambda: {"P": [], "B": []})
    Pg, Bg = set(), set()
    for r in rows:
        if r["expected_class"] not in ("P", "B"):
            continue
        y.append(1 if r["expected_class"] == "P" else 0)
        for k, c in cols.items():
            data[k].append(f(r, c))
        (Pg if r["expected_class"] == "P" else Bg).add(r["gene"])
        byg[r["gene"]][r["expected_class"]].append(f(r, "final_score"))

    print(f"file={path}  P={sum(y)}  B={len(y) - sum(y)}")

    # Mann-Whitney rank AUC (no sklearn dependency)
    def auc(x, yy):
        pairs = [(xi, yi) for xi, yi in zip(x, yy) if xi is not None]
        pos = [xi for xi, yi in pairs if yi == 1]
        neg = [xi for xi, yi in pairs if yi == 0]
        if not pos or not neg:
            return float("nan")
        c = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
        return c / (len(pos) * len(neg))

    print("\npooled ROC-AUC (P vs B):")
    for k in ("final", "LOF", "DN", "GOF"):
        print(f"  {k:6} {auc(data[k], y):.3f}")

    both = Pg & Bg
    print(f"\ngene overlap: {len(Pg)} P-genes, {len(Bg)} B-genes, {len(both)} with BOTH")
    if len(both) < 0.5 * min(len(Pg), len(Bg)):
        print("  ⚠️ low overlap → pooled AUC is confounded by cross-family score ranges")

    if both:
        print("\nwithin-gene separation (genes with both P and B):")
        deltas = []
        for g in sorted(both):
            pm, bm = st.mean(byg[g]["P"]), st.mean(byg[g]["B"])
            deltas.append(pm - bm)
            print(f"  {g:8} ΔP-B={pm - bm:+.3f}  {'✓' if pm > bm else '✗ INVERTED'}")
        print(f"  → P>B in {sum(d > 0 for d in deltas)}/{len(deltas)}, "
              f"mean Δ={st.mean(deltas):+.3f}")


if __name__ == "__main__":
    main()
