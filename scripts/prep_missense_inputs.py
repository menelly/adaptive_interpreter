#!/usr/bin/env python3
"""
Convert inputs_missense_only/<GENE>_input.tsv files into the run-harness format
(gene / hgvs_p / expected_class / clinical_sig).

KEY DESIGN (Ren, 2026-05-21): we RUN every missense variant through the pipeline
— including VUS — because classifying VUS is the actual job. But VUS get
expected_class="VUS" so eval_regression_auc.py scores them and then EXCLUDES them
from the ROC/AUC math (it only counts P and B). Ground-truth-labeled P/B drive the
measurement; VUS ride along to get predictions.

Usage:
  # report per-gene label counts (pick sanity genes), no output file:
  python scripts/prep_missense_inputs.py --report
  # build a labeled run-file from specific genes:
  python scripts/prep_missense_inputs.py --genes TP53,FBN1,HNF1B --out validation/gen_sanity.tsv
  # build from ALL genes that have data:
  python scripts/prep_missense_inputs.py --all --out validation/gen_all.tsv
"""
import csv
import glob
import os
import sys

INPUT_DIR = "/home/Ace/AdaptiveInterpreter/inputs_missense_only"

# 15-gene validation bed — EXCLUDE from a generalization test (these are in-bed)
VALIDATION_BED = {
    "ATP7B", "CDH23", "CFTR", "COL5A2", "CTNNB1", "DES", "GAA", "HEXA",
    "MYH9", "PDE6B", "POLG", "PTPN11", "SLC2A1", "TUBB3",
}

P_LABELS = {"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"}
B_LABELS = {"benign", "likely_benign", "benign/likely_benign"}


def classify(clinical_sig: str) -> str:
    s = (clinical_sig or "").strip().lower()
    if s in P_LABELS:
        return "P"
    if s in B_LABELS:
        return "B"
    return "VUS"  # uncertain / conflicting / other / not_provided → ride along, not counted


def read_gene(gene: str):
    path = os.path.join(INPUT_DIR, f"{gene}_input.tsv")
    if not os.path.exists(path):
        return []
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("molecular_consequence") or "").strip().lower() != "missense":
                continue
            hgvs_p = (r.get("hgvs_p") or "").strip()
            if not hgvs_p:
                continue
            rows.append({
                "gene": gene,
                "hgvs_p": hgvs_p,
                "expected_class": classify(r.get("clinical_sig")),
                "clinical_sig": (r.get("clinical_sig") or "").strip(),
            })
    # dedup on (gene, hgvs_p), preferring a labeled call over VUS if duplicated
    best = {}
    for r in rows:
        k = (r["gene"], r["hgvs_p"])
        if k not in best or (best[k]["expected_class"] == "VUS" and r["expected_class"] != "VUS"):
            best[k] = r
    return list(best.values())


def all_genes_with_data():
    genes = []
    for path in sorted(glob.glob(os.path.join(INPUT_DIR, "*_input.tsv"))):
        gene = os.path.basename(path)[:-len("_input.tsv")]
        with open(path) as fh:
            if sum(1 for _ in fh) > 1:  # has data beyond header
                genes.append(gene)
    return genes


def report():
    print(f"{'gene':12} {'P':>4} {'B':>4} {'VUS':>5} {'tot':>5}  in_bed")
    cand = []
    for gene in all_genes_with_data():
        rows = read_gene(gene)
        if not rows:
            continue
        p = sum(r["expected_class"] == "P" for r in rows)
        b = sum(r["expected_class"] == "B" for r in rows)
        v = sum(r["expected_class"] == "VUS" for r in rows)
        in_bed = "BED" if gene in VALIDATION_BED else ""
        print(f"{gene:12} {p:>4} {b:>4} {v:>5} {len(rows):>5}  {in_bed}")
        # good sanity candidate: both labels present, not giant, out of bed
        if p >= 3 and b >= 3 and len(rows) <= 400 and gene not in VALIDATION_BED:
            cand.append((gene, p, b, len(rows)))
    print("\n=== sanity candidates (>=3 P, >=3 B, <=400 rows, out-of-bed) ===")
    for gene, p, b, tot in sorted(cand, key=lambda x: x[3]):
        print(f"  {gene:12} P={p} B={b} tot={tot}")


def build(genes, out):
    all_rows = []
    for gene in genes:
        all_rows.extend(read_gene(gene))
    with open(out, "w", newline="") as ofh:
        w = csv.writer(ofh, delimiter="\t")
        w.writerow(["gene", "hgvs_p", "expected_class", "clinical_sig"])
        for r in all_rows:
            w.writerow([r["gene"], r["hgvs_p"], r["expected_class"], r["clinical_sig"]])
    p = sum(r["expected_class"] == "P" for r in all_rows)
    b = sum(r["expected_class"] == "B" for r in all_rows)
    v = sum(r["expected_class"] == "VUS" for r in all_rows)
    print(f"wrote {out}: {len(all_rows)} variants ({p} P, {b} B, {v} VUS) across {len(genes)} genes")


def main():
    args = sys.argv[1:]
    if "--report" in args:
        report()
        return
    out = None
    if "--out" in args:
        out = args[args.index("--out") + 1]
    if "--all" in args:
        genes = all_genes_with_data()
    elif "--genes" in args:
        genes = args[args.index("--genes") + 1].split(",")
    else:
        print(__doc__)
        return
    if not out:
        print("need --out PATH")
        return
    build(genes, out)


if __name__ == "__main__":
    main()
