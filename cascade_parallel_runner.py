#!/usr/bin/env python3
"""
Parallel runner for CASCADE batch processor.
Runs multiple genes concurrently using subprocess isolation.
"""

import argparse
import os
import sys
import time
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# 🍬 Force this script's own stdout/stderr to be unbuffered so candy streams
# even when output is piped to tee. Without this, prints sit in a 4KB buffer.
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

# Run from /home/Ace with PYTHONPATH set
ROOT = Path("/home/Ace")
PYTHONPATH = str(ROOT)


def run_single_gene(gene: str, input_file: str, output_file: str) -> tuple:
    """Run cascade_batch_processor for a single gene. Returns (gene, success, elapsed, error)."""
    
    cmd = [
        "python3", "-u",  # 🍬 unbuffered: candy streams immediately, doesn't sit in 4KB buffer
        "-m", "AdaptiveInterpreter.analyzers.cascade_batch_processor",
        "--input", input_file,
        "--output", output_file
    ]

    env = os.environ.copy()
    env["PYTHONPATH"] = PYTHONPATH
    env["PYTHONUNBUFFERED"] = "1"  # belt-and-suspenders: also tell Python to be unbuffered
    
    print(f"[{gene}] Starting...", flush=True)
    start = time.time()

    # 🍬 BUGFIX 2026-04-26 (Ace): old version captured stdout with subprocess.PIPE
    # which silently swallowed all the cascade emoji candy. Now we let stdout/stderr
    # inherit the parent terminal so all the analyzer output streams live, exactly
    # like running cascade_batch_processor directly. With multiple workers, lines
    # from different genes will interleave — that's the price of parallel candy.
    try:
        proc = subprocess.run(
            cmd, cwd=str(ROOT), env=env,
            timeout=7200,  # 2 hours per gene (big ones like BRCA2/MSH2/RYR1/CFTR need it)
        )
        elapsed = time.time() - start

        if proc.returncode != 0:
            print(f"[{gene}] FAILED in {elapsed:.1f}s", flush=True)
            return (gene, False, elapsed, f"returncode={proc.returncode}")

        print(f"[{gene}] ✅ Done in {elapsed:.1f}s", flush=True)
        return (gene, True, elapsed, "")
        
    except subprocess.TimeoutExpired:
        elapsed = time.time() - start
        print(f"[{gene}] TIMEOUT after {elapsed:.1f}s")
        return (gene, False, elapsed, "TIMEOUT")


def main():
    ap = argparse.ArgumentParser(description="Parallel CASCADE batch processor")
    ap.add_argument("--input-dir", required=True, help="Directory containing *_input.tsv files")
    ap.add_argument("--output-dir", required=True, help="Directory for *_cascade.tsv outputs")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 1),
                    help="Number of parallel workers (default: CPU-1)")
    ap.add_argument("--genes", nargs="*", help="Specific genes to run (default: all in input-dir)")
    args = ap.parse_args()

    # CRITICAL: Convert to absolute paths!
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find genes to process
    if args.genes:
        genes = args.genes
    else:
        genes = sorted([f.stem.replace("_input", "") for f in input_dir.glob("*_input.tsv")])
    
    if not genes:
        print("No genes found!")
        sys.exit(1)

    print(f"CASCADE PARALLEL RUNNER")
    print(f"  Input:   {input_dir}")
    print(f"  Output:  {output_dir}")
    print(f"  Workers: {args.workers}")
    print(f"  Genes:   {len(genes)}")
    print()

    start_all = time.time()
    results = []
    
    # Build file paths as strings (absolute paths!)
    jobs = []
    for g in genes:
        input_file = str(input_dir / f"{g}_input.tsv")
        output_file = str(output_dir / f"{g}_cascade.tsv")
        jobs.append((g, input_file, output_file))
    
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {
            ex.submit(run_single_gene, g, inf, outf): g 
            for g, inf, outf in jobs
        }
        for fut in as_completed(futs):
            gene = futs[fut]
            try:
                result = fut.result()
                results.append(result)
            except Exception as e:
                print(f"[{gene}] Exception: {e}")
                results.append((gene, False, 0, str(e)))

    elapsed_all = time.time() - start_all
    
    # Summary
    succeeded = [r for r in results if r[1]]
    failed = [r for r in results if not r[1]]
    
    print(f"\n{'='*50}")
    print(f"COMPLETE: {len(succeeded)}/{len(genes)} succeeded in {elapsed_all:.1f}s")
    
    if failed:
        print(f"\nFailed genes:")
        for gene, _, elapsed, err in failed:
            print(f"  {gene}: {err[:100] if err else 'unknown error'}")


if __name__ == "__main__":
    main()
