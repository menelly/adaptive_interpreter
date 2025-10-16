#!/usr/bin/env python3
import argparse
import csv
import os
from pathlib import Path

# Columns expected in discovery TSVs
REQUIRED_COLS = [
    "gene", "hgvs_p", "hgvs_g", "gnomad_freq", "clinical_sig",
    "molecular_consequence", "review_status"
]

NON_CODING = {
    "synonymous", "intron", "intronic", "utr", "5_prime_utr_variant", "3_prime_utr_variant",
    "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"
}

LOF = {"frameshift", "frameshift_variant", "stop_gained", "stop_gained_variant"}
SPLICE_KEYS = {"splice", "splice_acceptor", "splice_donor", "splice_region"}
TO_CONVERT = {"missense", "missense_variant", "stop_lost", "stop_lost_variant"}


def parse_freq(val: str) -> float:
    if val is None:
        return 0.0
    s = str(val).strip()
    if s == "" or s.lower() == "nan":
        return 0.0
    try:
        return float(s)
    except Exception:
        return 0.0


def norm_consequence(mc: str) -> str:
    if not mc:
        return "unknown"
    mc = mc.strip().lower()
    # ClinVar MC can be like "SO:0001583|missense_variant"; keep the trailing term
    if "|" in mc:
        mc = mc.split("|")[-1]
    return mc


def categorize(row: dict) -> str:
    mc = norm_consequence(row.get("molecular_consequence", ""))
    freq = parse_freq(row.get("gnomad_freq"))
    hgvs_p = (row.get("hgvs_p") or "").strip()

    if freq >= 0.05 or mc in NON_CODING:
        return "benign_prebucket"

    if mc in LOF:
        return "lof_precall"

    if any(k in mc for k in SPLICE_KEYS):
        return "splice_precall"

    if mc in TO_CONVERT:
        # Only send to conversion when we don't already have p.-notation
        if not hgvs_p:
            return "to_convert"
        else:
            return "protein_ready"  # already protein annotated; can be fed to cascade directly

    return "other"


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def write_rows(path: Path, header: list, rows: list):
    if not rows:
        return
    # Filter out None headers and unknown keys to avoid csv.DictWriter errors
    header = [h for h in header if h]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        for r in rows:
            safe_row = {k: r.get(k, "") for k in header}
            w.writerow(safe_row)


def process_tsv(tsv_path: Path, out_root: Path):
    # group by gene
    by_gene = {}

    with tsv_path.open() as f:
        # Discovery files are TSV by contract; avoid auto-detect because commas in values confuse sniffer
        reader = csv.DictReader(f, delimiter="\t")
        header = [h for h in (reader.fieldnames or []) if h]

        # sanity: if some required columns are missing, still operate best-effort
        for row in reader:
            gene = (row.get("gene") or "").strip()
            if not gene:
                gene = "_unknown"
            by_gene.setdefault(gene, {"header": header, "rows": []})
            by_gene[gene]["rows"].append(row)

    for gene, pack in by_gene.items():
        gene_dir = out_root / gene
        ensure_dir(gene_dir)

        header = pack["header"]
        rows = pack["rows"]

        benign, lof, splice, to_convert, protein_ready, other = [], [], [], [], [], []
        for r in rows:
            cat = categorize(r)
            if cat == "benign_prebucket":
                benign.append(r)
            elif cat == "lof_precall":
                lof.append(r)
            elif cat == "splice_precall":
                splice.append(r)
            elif cat == "to_convert":
                to_convert.append(r)
            elif cat == "protein_ready":
                protein_ready.append(r)
            else:
                other.append(r)

        base = gene
        write_rows(gene_dir / f"{base}.benign_prebucket.tsv", header, benign)
        write_rows(gene_dir / f"{base}.lof_precall.tsv", header, lof)
        write_rows(gene_dir / f"{base}.splice_precall.tsv", header, splice)
        write_rows(gene_dir / f"{base}.to_convert.discovery.tsv", header, to_convert)
        write_rows(gene_dir / f"{base}.protein_ready.discovery.tsv", header, protein_ready)
        # Keep the leftovers for audit if needed
        write_rows(gene_dir / f"{base}.other.tsv", header, other)


def main():
    ap = argparse.ArgumentParser(description="Split discovery TSVs into routing buckets")
    ap.add_argument("--input-dir", required=True, help="Folder containing *.discovery.tsv files")
    ap.add_argument("--output-dir", required=True, help="Folder to write per-gene splits")
    args = ap.parse_args()

    in_dir = Path(args.input_dir)
    out_dir = Path(args.output_dir)
    ensure_dir(out_dir)

    tsvs = [p for p in in_dir.glob("**/*.discovery.tsv") if p.is_file()]
    if not tsvs:
        print(f"No discovery TSVs found in {in_dir}")
        return

    for tsv in tsvs:
        print(f"Processing {tsv}")
        # output dir per source folder name (optional): we already group by gene inside
        process_tsv(tsv, out_dir)


if __name__ == "__main__":
    main()

