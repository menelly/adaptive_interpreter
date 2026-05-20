#!/usr/bin/env python3
"""
Run the frozen regression sample through the cascade and record results.

Use:
  1. Capture a BASELINE before a scoring change:
       python scripts/run_regression_sample.py validation/regression_sample_150.tsv baseline
  2. After the change, capture AGAIN:
       python scripts/run_regression_sample.py validation/regression_sample_150.tsv afterchange
  3. Diff the two output TSVs to see which calls moved (and check accuracy vs expected_class).

Runs the NEW run-all path (analyze_cascade) — the path the GOF interpretation gate lives in.
Writes results incrementally so a partial/interrupted run is still usable.
"""
import csv
import sys
import time

sys.path.insert(0, "/home/Ace")
from AdaptiveInterpreter.analyzers.cascade_analyzer import CascadeAnalyzer  # noqa: E402


def main():
    sample = sys.argv[1] if len(sys.argv) > 1 else "validation/regression_sample_150.tsv"
    tag = sys.argv[2] if len(sys.argv) > 2 else "baseline"
    out = f"validation/run_{tag}.tsv"

    with open(sample) as fh:
        variants = list(csv.DictReader(fh, delimiter="\t"))

    a = CascadeAnalyzer()
    t0 = time.time()
    with open(out, "w", newline="") as ofh:
        w = csv.writer(ofh, delimiter="\t")
        w.writerow(["gene", "hgvs_p", "expected_class", "final_score",
                    "final_classification", "DN", "LOF", "GOF", "status"])
        for i, row in enumerate(variants, 1):
            g, v, exp = row["gene"], row["hgvs_p"], row["expected_class"]
            try:
                r = a.analyze_cascade(g, v)
                s = r.get("scores", {})
                w.writerow([g, v, exp, f"{r.get('final_score', 0):.4f}",
                            r.get("final_classification", "?"),
                            f"{s.get('DN', 0):.4f}", f"{s.get('LOF', 0):.4f}",
                            f"{s.get('GOF', 0):.4f}", r.get("status", "?")])
            except Exception as e:
                w.writerow([g, v, exp, "", "ERROR", "", "", "", f"{type(e).__name__}: {e}"])
            ofh.flush()
            if i % 10 == 0:
                print(f"  {i}/{len(variants)}  ({time.time() - t0:.0f}s)", flush=True)
    print(f"done -> {out}  ({len(variants)} variants, {time.time() - t0:.0f}s total)", flush=True)


if __name__ == "__main__":
    main()
