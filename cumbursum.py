#!/usr/bin/env python3
"""CumBurSum — cumulative pathway burden analysis.

Cross-database pathway-burden analyzer for personal genomes. See
cumbursum/README.md for methodology and cumbursum/data/README.md for
required data files.

CLI:
    # CASCADE results (per-variant mechanism-scored)
    python3 cumbursum.py \\
        --cascade <subject>_cascade_results.tsv \\
        --person "Subject 001" \\
        --output report.md

    # Curated variant list (GENE/VARIANT/SYSTEMS/PRIORITY), optionally
    # with CASCADE overlay columns
    python3 cumbursum.py \\
        --curated my_variants.tsv \\
        --output report.md
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from cumbursum.burden_scorer import (
    MECHANISMS, compute_pathway_burdens, load_cascade_variants,
)
from cumbursum.consensus import build_consensus
from cumbursum.constraint_weights import ConstraintWeights
from cumbursum.curated_loader import load_curated_variants
from cumbursum.hallmarks_loader import HallmarksLoader
from cumbursum.null_distribution import (
    DEFAULT_N_PERM, benjamini_hochberg, permutation_significance,
)
from cumbursum.reactome_loader import ReactomeLoader
from cumbursum.report_generator import write_report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--cascade", help="CASCADE results TSV (per-variant mechanism scores)")
    src.add_argument("--curated", help="Curated 4-col TSV (GENE/VARIANT/SYSTEMS/PRIORITY); "
                                        "consequence-based burden applied for frameshift/nonsense/splice")
    ap.add_argument("--output", required=True, help="Output markdown report path")
    ap.add_argument("--person", default=None, help="Person name (defaults to input filename stem)")
    ap.add_argument("--min-score", type=float, default=0.3,
                    help="Drop variants with adj_score below this (default 0.3, curated ignores it)")
    ap.add_argument("--n-perm", type=int, default=DEFAULT_N_PERM,
                    help=f"Permutations per pathway (default {DEFAULT_N_PERM})")
    ap.add_argument("--top", type=int, default=20,
                    help="Top N pathways to show per mechanism ranking")
    ap.add_argument("--no-full-table", action="store_true",
                    help="Skip the full per-pathway table in the report")
    ap.add_argument("--include-mt", action="store_true",
                    help="Include MT-* genes (default: exclude; Dante alignment drift)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(levelname)s %(message)s",
    )

    input_path = Path(args.cascade or args.curated)
    input_kind = "cascade" if args.cascade else "curated"
    if not input_path.exists():
        print(f"ERROR: input file not found: {input_path}", file=sys.stderr)
        return 1

    person = args.person or input_path.stem.replace("_cascade_results", "")
    exclude_mt = not args.include_mt

    print(f"🧬 CumBurSum v0 running for {person} ({input_kind} input)")
    print(f"   Input: {input_path}")
    print(f"   Output: {args.output}")
    if input_kind == "cascade":
        print(f"   Perms: {args.n_perm}; min_score: {args.min_score}; exclude_mt: {exclude_mt}")
    else:
        print(f"   Perms: {args.n_perm}; exclude_mt: {exclude_mt}")
    print()

    # Load variants once.
    if input_kind == "cascade":
        variants = load_cascade_variants(
            input_path, exclude_mt=exclude_mt, min_score=args.min_score,
        )
    else:
        variants = load_curated_variants(input_path, exclude_mt=exclude_mt)

    if not variants:
        print("ERROR: no variants survived filtering.", file=sys.stderr)
        return 2

    # --- Reactome (with hierarchy rollup) ---
    print("📚 Loading Reactome (with hierarchy rollup)...")
    reactome = ReactomeLoader(exclude_mt=exclude_mt, rollup_hierarchy=True)
    reactome.load()
    r_burdens = compute_pathway_burdens(variants, reactome.pathway_to_genes)
    r_sig = permutation_significance(
        r_burdens, variants,
        gene_universe=sorted(reactome.gene_to_pathways.keys()),
        mechanisms=MECHANISMS, n_perm=args.n_perm, random_seed=args.seed,
    )
    benjamini_hochberg(r_sig, mechanisms=MECHANISMS)

    # --- MSigDB Hallmarks ---
    print("📚 Loading MSigDB Hallmarks...")
    hallmarks = HallmarksLoader(exclude_mt=exclude_mt)
    hallmarks.load()
    h_burdens = compute_pathway_burdens(
        variants, hallmarks.pathway_to_genes,
        min_pathway_size=5, max_pathway_size=5000,  # Hallmarks can be up to ~200; no upper limit issue
    )
    h_sig = permutation_significance(
        h_burdens, variants,
        gene_universe=sorted(hallmarks.gene_to_pathways.keys()),
        mechanisms=MECHANISMS, n_perm=args.n_perm, random_seed=args.seed,
    )
    benjamini_hochberg(h_sig, mechanisms=MECHANISMS)

    # --- Consensus ---
    print("🔗 Computing cross-database consensus...")
    themes = build_consensus(
        reactome_burdens=r_burdens, reactome_sig=r_sig,
        hallmark_burdens=h_burdens, hallmark_sig=h_sig,
        top_n_per_db=30,
        p_ceiling_reactome=0.01, p_ceiling_hallmarks=0.075,
        mechanism="total",
        min_gene_overlap=0.5,
    )

    # --- Gene constraint info (pLI, mis_oe) for display ---
    print("🧬 Loading gene constraint annotations...")
    constraint = ConstraintWeights()
    constraint_by_gene = {}
    for v in variants:
        if v.gene not in constraint_by_gene:
            info = constraint.info(v.gene)
            if info:
                constraint_by_gene[v.gene] = info

    # --- Report (leads with consensus, then per-DB breakdowns) ---
    out_path = write_report(
        Path(args.output),
        person_name=person,
        variants=variants,
        reactome_burdens=r_burdens, reactome_sig=r_sig,
        hallmark_burdens=h_burdens, hallmark_sig=h_sig,
        themes=themes,
        constraint_by_gene=constraint_by_gene,
        n_perm=args.n_perm,
        cascade_path=str(input_path),
        mechanisms=MECHANISMS,
        top_per_mechanism=args.top,
        include_full_table=not args.no_full_table,
    )

    # CLI summary.
    n_consensus = sum(1 for t in themes if t.is_consensus)
    print()
    print(f"✅ Report written: {out_path}")
    print(f"   Consensus themes (≥2 databases agree): {n_consensus}")
    print(f"   Single-DB themes: {len(themes) - n_consensus}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
