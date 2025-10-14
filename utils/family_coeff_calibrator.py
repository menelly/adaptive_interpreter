#!/usr/bin/env python3
"""
Family Coefficient Calibrator

Turns per-family training data into interpretable coefficients that replace
magic numbers in analyzers.

Outputs JSON per family at:
  AdaptiveInterpreter/cascade/resources/family_models/{family}_coefficients.json

Schema (example):
{
  "family": "collagen_fibrillar",
  "n_variants": 12345,
  "aa_effects": {
    "G": {"ref_loss_multiplier": 1.60, "gain_multiplier": 1.25, "confidence": 0.85, "n": 980},
    "R": {"ref_loss_multiplier": 1.30, "gain_multiplier": 1.40, "confidence": 0.72, "n": 640}
  },
  "domain_modifiers": {
    "gly_xy_repeat_boost": 1.30,  # optional, for collagens only when in Gly-X-Y region
    "notes": "Derived from conditional means with Laplace smoothing"
  },
  "conservation_mapping": {
    "method": "empirical_breakpoints_v1",
    "breakpoints": [1.0, 3.0, 5.0, 7.0],
    "multipliers": [1.0, 1.2, 1.5, 2.0, 2.5]
  }
}

Notes
- This module is intentionally decoupled from the trainers. Provide a dataframe
  with the required columns. A thin bridge in the training pipeline can call
  these functions after training/aggregation runs and write artifacts.
- Gentle, conservative defaults: we blend multipliers toward 1.0 by a
  confidence factor, and clamp outputs to avoid runaway effects.
"""
from __future__ import annotations
import json
from pathlib import Path
from typing import Dict, Any, Optional

import math
import numpy as np
import pandas as pd

AA_LIST = list("ARNDCEQGHILKMFPSTWYV")


def _laplace_rate(success: float, total: float, alpha: float = 1.0, beta: float = 1.0) -> float:
    """Laplace/Beta(1,1) smoothing for rates; returns (success+alpha)/(total+alpha+beta)."""
    return (success + alpha) / (total + alpha + beta) if total >= 0 else 0.5


def _safe_ratio(num: float, den: float, default: float = 1.0) -> float:
    if den <= 0:
        return default
    return num / den


def _blend_toward_one(mult: float, confidence: float) -> float:
    confidence = max(0.0, min(1.0, float(confidence)))
    return 1.0 + (mult - 1.0) * confidence


def _gentle_clamp(mult: float, lo: float = 0.7, hi: float = 1.8) -> float:
    return max(lo, min(hi, mult))


def calibrate_family_coefficients(df: pd.DataFrame, family: str, min_n: int = 30) -> Dict[str, Any]:
    """
    Compute per-AA coefficients for a given family from a labeled variant dataframe.

    Required columns in df (per variant row):
      - family (str)
      - ref_aa (one-letter)
      - alt_aa (one-letter)
      - one of: is_pathogenic ({0,1}) or pathogenicity_score (float in [0,1])
    Optional columns:
      - in_gly_x_y_repeat (bool/int)
      - one of: phylop (float) or phylop_score (float) for empirical conservation mapping

    Returns a JSON-serializable dict per schema above. Missing parts are omitted.
    """
    fam_df = df[df["family"].str.lower() == str(family).lower()].copy()
    out: Dict[str, Any] = {
        "family": family,
        "n_variants": int(len(fam_df)),
        "aa_effects": {},
    }

    if len(fam_df) == 0:
        return out

    # Normalize pathogenic label
    if "is_pathogenic" in fam_df.columns:
        fam_df["is_pathogenic"] = fam_df["is_pathogenic"].astype(float)
    elif "pathogenicity_score" in fam_df.columns:
        fam_df["is_pathogenic"] = (fam_df["pathogenicity_score"].astype(float) >= 0.5).astype(float)
    else:
        raise ValueError("DataFrame must contain is_pathogenic or pathogenicity_score column")

    # Baseline pathogenic rate for the family
    base_rate = _laplace_rate(fam_df["is_pathogenic"].sum(), len(fam_df))

    # --- Per-AA effects ---
    aa_effects: Dict[str, Any] = {}
    for aa in AA_LIST:
        sub_ref = fam_df[fam_df["ref_aa"] == aa]
        sub_alt = fam_df[fam_df["alt_aa"] == aa]
        n_ref = int(len(sub_ref))
        n_alt = int(len(sub_alt))

        if n_ref >= min_n:
            ref_rate = _laplace_rate(sub_ref["is_pathogenic"].sum(), n_ref)
            ref_mult = _safe_ratio(ref_rate, base_rate, default=1.0)
        else:
            ref_mult = 1.0

        if n_alt >= min_n:
            alt_rate = _laplace_rate(sub_alt["is_pathogenic"].sum(), n_alt)
            alt_mult = _safe_ratio(alt_rate, base_rate, default=1.0)
        else:
            alt_mult = 1.0

        # Confidence grows with sample size (log-shaped, capped at ~1 with 500+ samples)
        n = n_ref + n_alt
        conf = min(1.0, math.log10(max(10.0, n)) / 3.0)  # ~0.33 at n=100, ~0.66 at n=1000

        # Gentle clamp and record
        aa_effects[aa] = {
            "ref_loss_multiplier": round(_gentle_clamp(ref_mult), 3),
            "gain_multiplier": round(_gentle_clamp(alt_mult), 3),
            "confidence": round(conf, 3),
            "n": int(n),
        }

    out["aa_effects"] = aa_effects

    # --- Domain modifier: collagen Gly-X-Y repeat boost ---
    if "in_gly_x_y_repeat" in fam_df.columns:
        try:
            # Compare pathogenic rates for Gly ref-loss inside Gly-X-Y vs outside
            gly_ref = fam_df[(fam_df["ref_aa"] == "G")]
            in_rep = gly_ref[gly_ref["in_gly_x_y_repeat"] == 1]
            out_rep = gly_ref[gly_ref["in_gly_x_y_repeat"] == 0]
            if len(in_rep) >= min_n and len(out_rep) >= min_n:
                rate_in = _laplace_rate(in_rep["is_pathogenic"].sum(), len(in_rep))
                rate_out = _laplace_rate(out_rep["is_pathogenic"].sum(), len(out_rep))
                boost = _gentle_clamp(_safe_ratio(rate_in, rate_out, default=1.0), lo=1.0, hi=1.5)
                out["domain_modifiers"] = {
                    "gly_xy_repeat_boost": round(boost, 3),
                    "notes": "Derived from conditional means with Laplace smoothing",
                }
        except Exception:
            pass

    # --- Empirical conservation mapping (optional) ---
    # Normalize phylop column name if needed
    if "phylop" not in fam_df.columns and "phylop_score" in fam_df.columns:
        fam_df["phylop"] = fam_df["phylop_score"]

    if "phylop" in fam_df.columns:
        try:
            # Bin phylop by typical breakpoints and compute per-bin relative risk
            bins = [-21, 1.0, 3.0, 5.0, 7.0, 21]
            fam_df["ph_bin"] = pd.cut(fam_df["phylop"], bins=bins, labels=False, include_lowest=True)
            rr = []
            for b in sorted(fam_df["ph_bin"].dropna().unique()):
                sub = fam_df[fam_df["ph_bin"] == b]
                rate = _laplace_rate(sub["is_pathogenic"].sum(), len(sub))
                rr.append(_gentle_clamp(_safe_ratio(rate, base_rate, default=1.0), lo=0.8, hi=2.5))
            # Map to upper bin edges to align with cascade logic
            out["conservation_mapping"] = {
                "method": "empirical_breakpoints_v1",
                "breakpoints": [1.0, 3.0, 5.0, 7.0],
                "multipliers": [1.0] + [round(x, 2) for x in rr[1:]],
            }
        except Exception:
            pass

    return out


def save_family_coefficients(artifact: Dict[str, Any], out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    family = artifact.get("family", "unknown")
    path = out_dir / f"{family}_coefficients.json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(artifact, f, indent=2)
    return path


def calibrate_and_save(df: pd.DataFrame, family: str, resources_root: Optional[Path] = None, min_n: int = 30) -> Optional[Path]:
    """
    Convenience: calibrate a single family and write JSON next to cascade resources.
    """
    if resources_root is None:
        # Default relative to cascade package
        resources_root = Path(__file__).resolve().parent.parent / "cascade" / "resources" / "family_models"
    artifact = calibrate_family_coefficients(df, family, min_n=min_n)
    return save_family_coefficients(artifact, resources_root)

