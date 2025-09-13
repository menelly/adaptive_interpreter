"""
Weight/threshold tuner for Nova DN analyzer using expected_labels.json

Usage:
  python -m nova_dn.tune [--labels resources/expected_labels.json] \
      [--annotations resources/protein_annotations.json] \
      [--out-weights resources/tuned_weights.json] \
      [--out-threshold resources/tuned_threshold.json] \
      [--report resources/tuning_report.md]
"""
from __future__ import annotations
import json
import os
import itertools
from typing import Dict, List, Tuple

from .analyzer import NovaDNAnalyzer, parse_variant
from .context import load_annotations_json, build_position_context

# Minimal FASTA mapping for this repo
PROT_FASTA = {
    "TP53": "resources/tp53.fasta",
    "FGFR3": "resources/fgfr3.fasta",
    "COL1A1": "resources/col1a1.fasta",
    "VWF": "resources/vwf.fasta",
    "FBN1": "resources/fbn1.fasta",
}


def _read_fasta(path: str) -> str:
    seq_lines: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line.strip())
    return "".join(seq_lines)


def load_labels(path: str) -> List[Dict]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def predict_top(scores: Dict[str, float], threshold: float) -> str:
    mech, val = max(scores.items(), key=lambda kv: kv[1])
    return mech if val >= threshold else "none"


def evaluate(an: NovaDNAnalyzer, anns: Dict, labels: List[Dict], weights: Dict[str, float], threshold: float) -> Tuple[float, float, List[Dict]]:
    correct_mech = 0
    correct_path = 0
    rows = []
    for item in labels:
        var = item["variant"]
        prot = item["protein"]
        expect_mech = item["expected_top"]
        expect_path = bool(item.get("pathogenic", True))
        seq = _read_fasta(PROT_FASTA[prot])
        # Build context for this position
        _, pos1, _ = parse_variant(var)
        ctx = build_position_context(anns, prot, pos1)
        if weights:
            ctx.setdefault("_weights", {}).update(weights)
        res = an.analyze(seq, var, ctx)
        pred_mech_raw = res["top_mechanism"]
        pred_mech = pred_mech_raw if res["mechanism_scores"][pred_mech_raw] >= threshold else "none"
        is_path = pred_mech != "none"
        if pred_mech == expect_mech:
            correct_mech += 1
        if is_path == expect_path:
            correct_path += 1
        rows.append({
            "variant": var,
            "protein": prot,
            "expected_top": expect_mech,
            "pred_top": pred_mech,
            "pred_raw": pred_mech_raw,
            "scores": res["mechanism_scores"],
            "threshold": threshold,
            "pathogenic_expected": expect_path,
            "pathogenic_pred": is_path,
        })
    return correct_mech / len(labels), correct_path / len(labels), rows


def grid_search(an: NovaDNAnalyzer, anns: Dict, labels: List[Dict]):
    # Define a small search grid for key weights + global threshold
    grid = {
        "active_site_jamming.active_site_proximity": [0.25, 0.3, 0.35, 0.4],
        "active_site_jamming.flexible_loop": [-0.1, -0.15, -0.2],
        "active_site_jamming.aromatic_swap": [0.1, 0.15],
        "interface_poisoning.interface_likelihood": [0.35, 0.4, 0.45],
        "trafficking_maturation.disulfide_network_change": [0.5, 0.6, 0.7, 0.8, 0.9],
        "trafficking_maturation.NXS/T_gained": [0.3, 0.35, 0.4],
        "trafficking_maturation.NXS/T_lost": [0.25, 0.3, 0.35],
    }
    thresholds = [0.35, 0.4, 0.45, 0.5]

    keys = list(grid.keys())
    best = None
    best_detail = None
    # Iterate combinations (kept small)
    for values in itertools.product(*[grid[k] for k in keys]):
        weights = {k: v for k, v in zip(keys, values)}
        for thr in thresholds:
            mech_acc, path_acc, rows = evaluate(an, anns, labels, weights, thr)
            score = (mech_acc, path_acc)
            if best is None or score > best:
                best = score
                best_detail = (weights, thr, rows)
    return best_detail, best


def write_report(path: str, best_weights: Dict[str, float], threshold: float, mech_acc: float, path_acc: float, rows: List[Dict]):
    lines = []
    lines.append("# Tuning Report\n")
    lines.append(f"Mechanism accuracy: {mech_acc*100:.1f}%\n")
    lines.append(f"Pathogenicity accuracy: {path_acc*100:.1f}%\n\n")
    lines.append("## Weights\n")
    for k, v in sorted(best_weights.items()):
        lines.append(f"- {k}: {v}\n")
    lines.append(f"\n## Threshold\n- global: {threshold}\n\n")
    lines.append("## Per-variant\n")
    lines.append("| variant | protein | expected | pred | raw_top | interface | active_site | lattice | trafficking |\n")
    lines.append("|---|---|---|---|---|---:|---:|---:|---:|\n")
    for r in rows:
        sc = r["scores"]
        lines.append(
            f"| {r['variant']} | {r['protein']} | {r['expected_top']} | {r['pred_top']} | {r['pred_raw']} | "
            f"{sc['interface_poisoning']:.2f} | {sc['active_site_jamming']:.2f} | {sc['lattice_disruption']:.2f} | {sc['trafficking_maturation']:.2f} |\n"
        )
    with open(path, "w", encoding="utf-8") as f:
        f.write("".join(lines))


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Tune feature weights and threshold against expected labels")
    ap.add_argument("--labels", default="resources/expected_labels.json")
    ap.add_argument("--annotations", default="resources/protein_annotations.json")
    ap.add_argument("--out-weights", default="resources/tuned_weights.json")
    ap.add_argument("--out-threshold", default="resources/tuned_threshold.json")
    ap.add_argument("--report", default="resources/tuning_report.md")
    args = ap.parse_args()

    labels = load_labels(args.labels)
    anns = load_annotations_json(args.annotations)
    an = NovaDNAnalyzer()

    (best_weights, best_thr, rows), (mech_acc, path_acc) = grid_search(an, anns, labels)

    # Save artifacts
    with open(args.out_weights, "w", encoding="utf-8") as f:
        json.dump(best_weights, f, indent=2)
    with open(args.out_threshold, "w", encoding="utf-8") as f:
        json.dump({"global": best_thr}, f, indent=2)
    write_report(args.report, best_weights, best_thr, mech_acc, path_acc, rows)

    print(json.dumps({
        "mechanism_accuracy": mech_acc,
        "pathogenicity_accuracy": path_acc,
        "threshold": best_thr,
        "weights": best_weights,
        "report": args.report,
    }, indent=2))


if __name__ == "__main__":
    main()

