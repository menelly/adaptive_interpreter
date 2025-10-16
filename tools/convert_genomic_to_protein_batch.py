#!/usr/bin/env python3
import argparse
import csv
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

NC_TO_CHR = {
    # GRCh38 RefSeq accessions mapping
    "NC_000001.11": "1", "NC_000002.12": "2", "NC_000003.12": "3", "NC_000004.12": "4",
    "NC_000005.10": "5", "NC_000006.12": "6", "NC_000007.14": "7", "NC_000008.11": "8",
    "NC_000009.12": "9", "NC_000010.11": "10", "NC_000011.10": "11", "NC_000012.12": "12",
    "NC_000013.11": "13", "NC_000014.9": "14", "NC_000015.10": "15", "NC_000016.10": "16",
    "NC_000017.11": "17", "NC_000018.10": "18", "NC_000019.10": "19", "NC_000020.11": "20",
    "NC_000021.9": "21", "NC_000022.11": "22", "NC_000023.11": "X", "NC_000024.10": "Y",
    "NC_012920.1": "M",
}

HGVS_G_RE = re.compile(r"^(?P<acc>NC_\d+\.\d+):g\.(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$")


def parse_hgvs_g(hgvs_g: str) -> Tuple[str, int, str, str]:
    m = HGVS_G_RE.match(hgvs_g)
    if not m:
        raise ValueError(f"Unsupported hgvs_g format: {hgvs_g}")
    acc = m.group("acc")
    chrom = NC_TO_CHR.get(acc)
    if not chrom:
        raise ValueError(f"Unknown accession mapping for {acc}")
    pos = int(m.group("pos"))
    ref = m.group("ref")
    alt = m.group("alt")
    return chrom, pos, ref, alt


def write_vcf(records: List[Tuple[str, int, str, str, str]], path: Path):
    with path.open("w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for rid, (chrom, pos, ref, alt, info) in enumerate(records, start=1):
            vid = f"var{rid}"
            f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t{info}\n")


def run_snpeff(vcf_in: Path, vcf_out: Path, snpeff_config: Path, genome: str):
    cmd = [
        "snpEff", "-Xmx8g", "-hgvs", "-canon",
        "-c", str(snpeff_config), genome, str(vcf_in)
    ]
    with vcf_out.open("w") as out:
        subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.PIPE)


def parse_ann_field(info: str) -> str:
    # ANN fields like: ANN=ALT|consequence|impact|gene|gene_id|feature_type|transcript|biotype|rank/total|...|HGVS.c|HGVS.p|...
    # Return HGVS.p first available
    ann_prefix = "ANN="
    if "ANN=" not in info:
        return ""
    # info could be many semicolon-separated fields, find ANN
    ann_field = None
    for piece in info.split(";"):
        if piece.startswith("ANN="):
            ann_field = piece[len(ann_prefix):]
            break
    if not ann_field:
        return ""
    # multiple comma-separated annotations
    entries = ann_field.split(",")
    for entry in entries:
        cols = entry.split("|")
        if len(cols) >= 11:
            hgvs_p = cols[10] or ""
            if hgvs_p:
                return hgvs_p
    return ""


def annotate_rows(tsv_in: Path, tsv_out: Path, snpeff_config: Path, genome: str):
    rows: List[dict] = []
    header: List[str] = []
    with tsv_in.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames or []
        for r in reader:
            rows.append(r)

    # Build minimal VCF
    recs = []
    idx_map = []  # map VCF order back to row index
    for i, r in enumerate(rows):
        hgvs_g = (r.get("hgvs_g") or "").strip()
        if not hgvs_g:
            continue
        try:
            chrom, pos, ref, alt = parse_hgvs_g(hgvs_g)
            # Keep the original hgvs_g to carry through in INFO for fallback mapping
            recs.append((chrom, pos, ref, alt, f"ORIGHGVS={hgvs_g}"))
            idx_map.append(i)
        except Exception:
            # skip rows we cannot parse
            continue

    if not recs:
        # nothing to do; just copy
        with tsv_out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow(r)
        return

    with tempfile.TemporaryDirectory() as td:
        vcf_in = Path(td) / "input.vcf"
        vcf_out = Path(td) / "annotated.vcf"
        write_vcf(recs, vcf_in)
        run_snpeff(vcf_in, vcf_out, snpeff_config, genome)

        # Map back annotations
        ann_by_key: Dict[Tuple[str, int, str, str], str] = {}
        with vcf_out.open() as f:
            for line in f:
                if not line:
                    continue
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 8:
                    # skip malformed or empty lines
                    continue
                chrom, pos_s, _vid, ref, alt, _qual, _filter, info = parts[:8]
                try:
                    pos = int(pos_s)
                except ValueError:
                    continue
                hgvs_p = parse_ann_field(info)
                ann_by_key[(chrom, pos, ref, alt)] = hgvs_p

    # Update rows
    if "hgvs_p" not in header:
        header.append("hgvs_p")

    for i, r in enumerate(rows):
        hgvs_g = (r.get("hgvs_g") or "").strip()
        try:
            chrom, pos, ref, alt = parse_hgvs_g(hgvs_g)
            key = (chrom, pos, ref, alt)
            hgvs_p = ann_by_key.get(key, "")
        except Exception:
            hgvs_p = ""
        if hgvs_p:
            r["hgvs_p"] = hgvs_p

    with tsv_out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    ap = argparse.ArgumentParser(description="Batch convert genomic HGVS to protein using SnpEff")
    ap.add_argument("--input-dir", required=True, help="Root folder with per-gene .to_convert.discovery.tsv files")
    ap.add_argument("--snpeff-config", required=True, help="Path to snpEff.config")
    ap.add_argument("--genome", required=True, help="Genome key configured in snpEff.config (e.g., GRCh38.gencode.v46)")
    ap.add_argument("--max-workers", type=int, default=4, help="Parallel workers (not used in first cut)")
    args = ap.parse_args()

    root = Path(args.input_dir)
    conf = Path(args.snpeff_config)
    genome = args.genome

    genes = [d for d in root.iterdir() if d.is_dir()]
    for gdir in genes:
        gene = gdir.name
        in_tsv = gdir / f"{gene}.to_convert.discovery.tsv"
        out_tsv = gdir / f"{gene}.protein_ready.discovery.tsv"
        if not in_tsv.exists():
            continue
        print(f"[SnpEff] {gene}: annotating {in_tsv}")
        annotate_rows(in_tsv, out_tsv, conf, genome)


if __name__ == "__main__":
    main()

