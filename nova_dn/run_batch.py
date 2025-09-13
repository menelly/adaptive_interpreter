"""
Minimal batch runner for Nova DN Analyzer.
Input: JSON list of items with fields:
  - seq_file: path to FASTA/plain AA
  - variant: e.g., "p.R273H"
  - protein: e.g., "TP53"
  - annotations_json: optional path to annotations JSON
Outputs a compact markdown table and JSONL to stdout.
"""
from __future__ import annotations
import sys, json, subprocess
from typing import List, Dict

PY = sys.executable or "python3"

USAGE = "usage: python3 -m nova_dn.run_batch path/to/batch.json"


def main(argv: List[str] | None = None) -> int:
    argv = argv or sys.argv[1:]
    if not argv:
        print(USAGE, file=sys.stderr)
        return 2
    path = argv[0]
    with open(path, "r", encoding="utf-8") as f:
        items: List[Dict] = json.load(f)

    print("| variant | protein | top | interface | active_site | lattice | trafficking | because |\n|---|---|---|---:|---:|---:|---:|---|")
    # If a tuned weights file exists and none specified in items, use it by default
    default_weights = "resources/tuned_weights.json"
    for it in items:
        cmd = [PY, "-m", "nova_dn.analyzer",
               "--seq-file", it["seq_file"],
               "--variant", it["variant"]]
        if it.get("annotations_json") and it.get("protein"):
            cmd += ["--annotations-json", it["annotations_json"], "--protein", it["protein"]]
        # Allow per-item weights_json override; else use default if present
        wj = it.get("weights_json")
        if not wj and __import__('os').path.exists(default_weights):
            wj = default_weights
        if wj:
            cmd += ["--weights-json", wj]
        if it.get("gene"):
            cmd += ["--gene", it["gene"]]
        if it.get("uniprot"):
            cmd += ["--uniprot", it["uniprot"]]
        cmd += ["--json"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print(f"| {it['variant']} | {it.get('protein','')} | ERROR | - | - | - | - | {proc.stderr.strip()} |")
            continue
        data = json.loads(proc.stdout)
        ms = data["mechanism_scores"]
        print(f"| {data['variant']} | {it.get('protein','')} | {data['top_mechanism']} | {ms['interface_poisoning']:.2f} | {ms['active_site_jamming']:.2f} | {ms['lattice_disruption']:.2f} | {ms['trafficking_maturation']:.2f} | {data['explanation']} |")
        # Also emit JSONL for programmatic use
        print(json.dumps(data), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

