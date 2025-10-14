#!/usr/bin/env python3
"""
Compute per-gene classification metrics from discovery_results_*.tsv files and aggregate.

- Reads one or more TSVs (each produced by lumen_discovery_pipeline.py for a single gene)
- Computes per-gene confusion matrix (bucketed) and rates using ground-truth agreement buckets
- Aggregates across all provided files via count-summing, then recomputes rates
- Optionally writes a per-row DISAGREE export with gene labels for audit

Usage:
  python3 AdaptiveInterpreter/metrics/per_gene_metrics.py discovery_results_*.tsv \
      --per_gene_out per_gene_metrics.tsv \
      --aggregate_out aggregate_metrics.tsv \
      --disagree_out per_gene_disagrees.tsv

Notes:
- Uses agreement_analysis.compare / normalize_* if available for perfect parity with the paper
- Falls back to internal normalizers if import fails
"""
from __future__ import annotations
import argparse
import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple


import pandas as pd

# Try to use the ground-truth logic if available
try:
    from agreement_analysis import compare as agree_compare
    from agreement_analysis import normalize_clinvar_sig, normalize_ai_class
except Exception:
    agree_compare = None

    def normalize_clinvar_sig(x: str) -> str:
        if x is None:
            return "VUS"
        s = str(x).strip().upper()
        # Common ClinVar values
        if any(k in s for k in ["PATHOGENIC", "LIKELY PATHOGENIC"]):
            # Treat LP/P as separate labels but both in P-bucket
            return "P" if "PATHOGENIC" in s and "LIKELY" not in s else "LP"
        if "BENIGN" in s:
            return "B" if "LIKELY" not in s else "LB"
        if "UNCERTAIN" in s or s in {"VUS", "UNKNOWN"}:
            return "VUS"
        return "VUS"

    def normalize_ai_class(x: str) -> str:
        if x is None:
            return "VUS"
        s = str(x).strip().upper()
        # Expect from our pipeline one of {P, LP, VUS-P, VUS, LB, B}
        # Be forgiving
        mapping = {
            "PATHOGENIC": "P",
            "LIKELY PATHOGENIC": "LP",
            "VUS-P": "VUS-P",
            "VUS": "VUS",
            "LIKELY BENIGN": "LB",
            "BENIGN": "B",
        }
        if s in mapping:
            return mapping[s]
        if s in {"P", "LP", "VUS-P", "VUS", "LB", "B"}:
            return s
        return "VUS"

# Buckets
POS_BUCKET = {"P", "LP", "VUS-P"}
NEG_BUCKET = {"VUS", "LB", "B"}


def pick_column(cols: List[str], candidates: List[str]) -> str | None:
    """Pick first column whose lowercase contains any candidate substrings in order."""
    low = [c.lower() for c in cols]
    for cand in candidates:
        cand_l = cand.lower()
        for i, name in enumerate(low):
            if cand_l in name:
                return cols[i]
    return None


def detect_columns(df: pd.DataFrame) -> Tuple[str, str]:
    """Return (clinvar_col, ai_col) best guesses."""
    cols = list(df.columns)
    clin = pick_column(cols, [
        "clinvar_sig", "clinvar_significance", "clinvar_clinical_significance",
        "clin_sig", "clinvar", "clinical_significance"
    ])
    ai = None
    # Prefer 'cascade_class' then 'final_classification'; fallback to variants and fuzzy
    for cand in ["cascade_class", "final_classification", "ai_class", "classification", "final", "ai"]:
        if cand in df.columns:
            ai = cand
            break
    if ai is None:
        # fuzzy pick by substring if exact not found
        ai = pick_column(cols, ["cascade_class", "final_classification", "cascade", "ai_class", "classification"])
    if clin is None:
        raise ValueError("Could not detect ClinVar significance column. Known examples: 'clinvar_sig', 'clinvar_clinical_significance'.")
    if ai is None:
        raise ValueError("Could not detect AI classification column. Known examples: 'cascade_class', 'final_classification'.")
    return clin, ai


def is_synonymous_hgvs_p(hgvs_p: str) -> bool:
    """Detect synonymous from protein HGVS (p.XnX, p.XxxnXxx, or equals form)."""
    if not isinstance(hgvs_p, str):
        return False
    try:
        m1 = re.search(r"p\.([A-Z\*])(\d+)([A-Z\*])", hgvs_p)
        if m1:
            ref, _, alt = m1.groups()
            if ref == alt and alt != '*':
                return True
        m2 = re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=)", hgvs_p)
        if m2:
            ref_long, _, alt_long = m2.groups()
            if alt_long == '=' or ref_long == alt_long:
                return True
    except Exception:
        return False
    return False


def summarize_file(path: Path, prefer_cv: str | None = None, prefer_ai: str | None = None) -> Tuple[pd.DataFrame, str, str]:
    df = pd.read_csv(path, sep='\t', dtype=str)
    if prefer_cv and prefer_cv in df.columns:
        clin_col = prefer_cv
    else:
        clin_col = None
    if prefer_ai and prefer_ai in df.columns:
        ai_col = prefer_ai
    else:
        ai_col = None
    if clin_col is None or ai_col is None:
        det_clin, det_ai = detect_columns(df)
        clin_col = clin_col or det_clin
        ai_col = ai_col or det_ai
    return df, clin_col, ai_col


def to_bucket(label: str, positive: bool) -> bool:
    return label in POS_BUCKET if positive else label in NEG_BUCKET


def compare_bucketed(cv: str, ai: str) -> str:
    """Return AGREE / DISAGREE / BETTER_DATA using ground-truth if available; else emulate."""
    cv_n = normalize_clinvar_sig(cv)
    ai_n = normalize_ai_class(ai)
    if agree_compare is not None:
        return agree_compare(cv_n, ai_n)
    # Emulate: VUS<->VUS-P counts as AGREE; otherwise same-bucket agree
    if (cv_n == "VUS" and ai_n == "VUS-P") or (cv_n == "VUS-P" and ai_n == "VUS"):
        return "AGREE"
    cv_pos = cv_n in POS_BUCKET
    ai_pos = ai_n in POS_BUCKET
    if cv_n == "VUS" and ai_n in {"B", "LB", "P", "LP"}:
        return "BETTER_DATA"
    return "AGREE" if cv_pos == ai_pos else "DISAGREE"


def agreement_counts(cv_series: pd.Series, ai_series: pd.Series) -> Dict[str, int]:
    AG = BD = DG = 0
    n = 0
    for cv, ai in zip(cv_series, ai_series):
        n += 1
        kind = compare_bucketed(cv, ai)
        if kind == "AGREE":
            AG += 1
        elif kind == "BETTER_DATA":
            BD += 1
        elif kind == "DISAGREE":
            DG += 1
    return {"AGREE": AG, "BETTER_DATA": BD, "DISAGREE": DG, "N": n}


def confusion_counts_strict(cv_series: pd.Series, ai_series: pd.Series) -> Dict[str, int]:
    """Strict clinical: binary only on {P,LP} vs {B,LB}. Drop any rows with VUS/VUS-P on either side."""
    TP = FP = TN = FN = 0
    for cv, ai in zip(cv_series, ai_series):
        cv_n = normalize_clinvar_sig(cv)
        ai_n = normalize_ai_class(ai)
        if cv_n in {"VUS", "VUS-P"} or ai_n in {"VUS", "VUS-P"}:
            continue  # exclude from strict binary
        cv_pos = cv_n in {"P", "LP"}
        ai_pos = ai_n in {"P", "LP"}
        if cv_pos and ai_pos:
            TP += 1
        elif (not cv_pos) and (not ai_pos):
            TN += 1
        elif (not cv_pos) and ai_pos:
            FP += 1
        else:
            FN += 1
    return {"TP": TP, "FP": FP, "TN": TN, "FN": FN}


def confusion_counts_lenient(cv_series: pd.Series, ai_series: pd.Series) -> Dict[str, int]:
    """Lenient directional: POS={P,LP,VUS-P}, NEG={VUS,LB,B}, but drop VUS↔VUS-P pairs."""
    TP = FP = TN = FN = 0
    for cv, ai in zip(cv_series, ai_series):
        cv_n = normalize_clinvar_sig(cv)
        ai_n = normalize_ai_class(ai)
        # drop special equivalence pairs from binary
        if (cv_n == "VUS" and ai_n == "VUS-P") or (cv_n == "VUS-P" and ai_n == "VUS"):
            continue
        cv_pos = cv_n in POS_BUCKET
        ai_pos = ai_n in POS_BUCKET
        if cv_pos and ai_pos:
            TP += 1
        elif (not cv_pos) and (not ai_pos):
            TN += 1
        elif (not cv_pos) and ai_pos:
            FP += 1
        else:
            FN += 1
    return {"TP": TP, "FP": FP, "TN": TN, "FN": FN}


def rates_binary(TP: int, FP: int, TN: int, FN: int) -> Dict[str, float]:
    def safe(num, den):
        return float(num) / den if den else 0.0
    return {
        "sensitivity": safe(TP, TP + FN),
        "specificity": safe(TN, TN + FP),
        "PPV": safe(TP, TP + FP),
        "NPV": safe(TN, TN + FN),
    }




def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+", help="discovery_results_*.tsv files")
    ap.add_argument("--per_gene_out", default="per_gene_metrics.tsv")
    ap.add_argument("--aggregate_out", default="aggregate_metrics.tsv")
    ap.add_argument("--disagree_out", default="per_gene_disagrees.tsv")
    ap.add_argument("--clinvar_col", default=None, help="Override ClinVar significance column name")
    ap.add_argument("--ai_col", default=None, help="Override AI classification column name (e.g., cascade_class)")
    ap.add_argument("--exclude_synonymous", action="store_true", help="Exclude p.= and p.XnX synonyms from metrics (export separately)")
    ap.add_argument("--syn_out", default="synonymous_excluded.tsv", help="Where to write excluded synonymous rows")
    args = ap.parse_args()

    rows = []
    agg_counts = {k: 0 for k in ["TP","FP","TN","FN","AGREE","BETTER_DATA","DISAGREE","N"]}
    disagree_frames = []

    for inp in args.inputs:
        p = Path(inp)
        if not p.exists():
            print(f"[skip] missing file: {p}")
            continue
        m = re.match(r"discovery_results_(.+)\.tsv$", p.name)
        gene = m.group(1) if m else p.stem
        print(f"[metrics] {gene} <- {p}")
        df, clin_col, ai_col = summarize_file(p, args.clinvar_col, args.ai_col)

        # Optional: exclude synonymous rows from metrics (still export separately)
        if args.exclude_synonymous and 'hgvs_p' in df.columns:
            syn_mask = df['hgvs_p'].map(is_synonymous_hgvs_p)
            if syn_mask.any():
                syn_df = df.loc[syn_mask].copy()
                syn_df['gene'] = gene
                cols = ['gene'] + [c for c in syn_df.columns if c != 'gene']
                syn_df = syn_df[cols]
                # append to global list via special marker in disagree_frames (reuse, we concatenate later)
                disagree_frames.append(syn_df.assign(__is_synonym__=True))
            df = df.loc[~syn_mask].copy()
            # re-detect columns if structure changed
            clin_col, ai_col = detect_columns(df)

        # Agreement/buckets (for AGREE/BETTER_DATA/DISAGREE + N)
        agree_ct = agreement_counts(df[clin_col], df[ai_col])

        # Strict and lenient binary confusion
        strict_ct = confusion_counts_strict(df[clin_col], df[ai_col])
        lenient_ct = confusion_counts_lenient(df[clin_col], df[ai_col])

        strict_rates = rates_binary(**strict_ct)
        lenient_rates = rates_binary(**lenient_ct)

        # accumulate per-gene
        row = {
            "gene": gene,
            **agree_ct,
            "strict_TP": strict_ct["TP"], "strict_FP": strict_ct["FP"], "strict_TN": strict_ct["TN"], "strict_FN": strict_ct["FN"],
            "strict_sensitivity": strict_rates["sensitivity"], "strict_specificity": strict_rates["specificity"],
            "strict_PPV": strict_rates["PPV"], "strict_NPV": strict_rates["NPV"],
            "lenient_TP": lenient_ct["TP"], "lenient_FP": lenient_ct["FP"], "lenient_TN": lenient_ct["TN"], "lenient_FN": lenient_ct["FN"],
            "lenient_sensitivity": lenient_rates["sensitivity"], "lenient_specificity": lenient_rates["specificity"],
            "lenient_PPV": lenient_rates["PPV"], "lenient_NPV": lenient_rates["NPV"],
        }
        rows.append(row)

        # sum aggregates
        for k in ["AGREE","BETTER_DATA","DISAGREE","N"]:
            agg_counts[k] += agree_ct[k]
        # add binary sums tracked separately
        agg_counts.setdefault("strict_TP",0); agg_counts["strict_TP"] += strict_ct["TP"]
        agg_counts.setdefault("strict_FP",0); agg_counts["strict_FP"] += strict_ct["FP"]
        agg_counts.setdefault("strict_TN",0); agg_counts["strict_TN"] += strict_ct["TN"]
        agg_counts.setdefault("strict_FN",0); agg_counts["strict_FN"] += strict_ct["FN"]
        agg_counts.setdefault("lenient_TP",0); agg_counts["lenient_TP"] += lenient_ct["TP"]
        agg_counts.setdefault("lenient_FP",0); agg_counts["lenient_FP"] += lenient_ct["FP"]
        agg_counts.setdefault("lenient_TN",0); agg_counts["lenient_TN"] += lenient_ct["TN"]
        agg_counts.setdefault("lenient_FN",0); agg_counts["lenient_FN"] += lenient_ct["FN"]

        # collect disagrees (on the df used for metrics)
        if args.clinvar_col and args.clinvar_col in df.columns:
            clin_col = args.clinvar_col
        else:
            clin_col = None
        if args.ai_col and args.ai_col in df.columns:
            ai_col = args.ai_col
        else:
            ai_col = None
        if clin_col is None or ai_col is None:
            det_clin, det_ai = detect_columns(df)
            clin_col = clin_col or det_clin
            ai_col = ai_col or det_ai
        mask = [compare_bucketed(cv, ai) == "DISAGREE" for cv, ai in zip(df[clin_col], df[ai_col])]
        if any(mask):
            dd = df.loc[mask].copy()
            # Ensure a 'gene' column exists (overwrite to be safe) and move it to front
            dd['gene'] = gene
            cols = ['gene'] + [c for c in dd.columns if c != 'gene']
            dd = dd[cols]
            disagree_frames.append(dd)

    if not rows:
        print("No inputs found; exiting.")
        sys.exit(1)

    per_gene_df = pd.DataFrame(rows)
    per_gene_df.sort_values(["DISAGREE", "BETTER_DATA"], ascending=[False, False], inplace=True)
    per_gene_df.to_csv(args.per_gene_out, sep='\t', index=False)
    print(f"[write] per-gene metrics -> {args.per_gene_out}")

    # Aggregate: compute strict and lenient rates from summed counts, plus agreement rates
    strict_rates_agg = rates_binary(
        agg_counts.get("strict_TP",0), agg_counts.get("strict_FP",0), agg_counts.get("strict_TN",0), agg_counts.get("strict_FN",0)
    )
    lenient_rates_agg = rates_binary(
        agg_counts.get("lenient_TP",0), agg_counts.get("lenient_FP",0), agg_counts.get("lenient_TN",0), agg_counts.get("lenient_FN",0)
    )

    agg_row = {
        "gene": "ALL",
        "AGREE": agg_counts["AGREE"], "BETTER_DATA": agg_counts["BETTER_DATA"], "DISAGREE": agg_counts["DISAGREE"], "N": agg_counts["N"],
        "strict_TP": agg_counts.get("strict_TP",0), "strict_FP": agg_counts.get("strict_FP",0), "strict_TN": agg_counts.get("strict_TN",0), "strict_FN": agg_counts.get("strict_FN",0),
        "strict_sensitivity": strict_rates_agg["sensitivity"], "strict_specificity": strict_rates_agg["specificity"],
        "strict_PPV": strict_rates_agg["PPV"], "strict_NPV": strict_rates_agg["NPV"],
        "lenient_TP": agg_counts.get("lenient_TP",0), "lenient_FP": agg_counts.get("lenient_FP",0), "lenient_TN": agg_counts.get("lenient_TN",0), "lenient_FN": agg_counts.get("lenient_FN",0),
        "lenient_sensitivity": lenient_rates_agg["sensitivity"], "lenient_specificity": lenient_rates_agg["specificity"],
        "lenient_PPV": lenient_rates_agg["PPV"], "lenient_NPV": lenient_rates_agg["NPV"],
    }
    pd.DataFrame([agg_row]).to_csv(args.aggregate_out, sep='\t', index=False)
    print(f"[write] aggregate metrics -> {args.aggregate_out}")

    if disagree_frames:
        # Split out synonymous exports (if present) vs actual DISAGREE rows
        dfs = []
        syn_list = []
        for d in disagree_frames:
            if isinstance(d, pd.DataFrame) and '__is_synonym__' in d.columns:
                syn_list.append(d.drop(columns=['__is_synonym__']))
            else:
                dfs.append(d)
        if dfs:
            all_dis = pd.concat(dfs, ignore_index=True)
            all_dis.to_csv(args.disagree_out, sep='\t', index=False)
            print(f"[write] per-row DISAGREE export -> {args.disagree_out}  (rows={len(all_dis)})")
        if syn_list:
            syn_all = pd.concat(syn_list, ignore_index=True)
            syn_all.to_csv(args.syn_out, sep='\t', index=False)
            print(f"[write] synonymous excluded export -> {args.syn_out}  (rows={len(syn_all)})")

if __name__ == "__main__":
    main()

