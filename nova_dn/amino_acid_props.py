"""Amino acid property tables and simple delta helpers (no external deps).
Properties include Kyte-Doolittle hydropathy, Zamyatnin-like volume (A^3),
net charge at ~physiological pH (approx), polarity class, and aromatic flag.
"""
from __future__ import annotations

AA_PROPS = {
    # hydropathy (Kyte-Doolittle), volume (A^3 approx), charge (-1,0,+1), polarity (0 nonpolar, 1 polar), aromatic(0/1)
    "A": {"hyd": 1.8,  "vol":  88, "chg": 0, "pol": 0, "aro": 0, "name": "Ala"},
    "R": {"hyd": -4.5, "vol": 173, "chg": +1, "pol": 1, "aro": 0, "name": "Arg"},
    "N": {"hyd": -3.5, "vol": 114, "chg": 0, "pol": 1, "aro": 0, "name": "Asn"},
    "D": {"hyd": -3.5, "vol": 111, "chg": -1, "pol": 1, "aro": 0, "name": "Asp"},
    "C": {"hyd": 2.5,  "vol": 108, "chg": 0, "pol": 0, "aro": 0, "name": "Cys"},
    "Q": {"hyd": -3.5, "vol": 143, "chg": 0, "pol": 1, "aro": 0, "name": "Gln"},
    "E": {"hyd": -3.5, "vol": 138, "chg": -1, "pol": 1, "aro": 0, "name": "Glu"},
    "G": {"hyd": -0.4, "vol":  60, "chg": 0, "pol": 0, "aro": 0, "name": "Gly"},
    "H": {"hyd": -3.2, "vol": 153, "chg": +1, "pol": 1, "aro": 1, "name": "His"},
    "I": {"hyd": 4.5,  "vol": 166, "chg": 0, "pol": 0, "aro": 0, "name": "Ile"},
    "L": {"hyd": 3.8,  "vol": 166, "chg": 0, "pol": 0, "aro": 0, "name": "Leu"},
    "K": {"hyd": -3.9, "vol": 168, "chg": +1, "pol": 1, "aro": 0, "name": "Lys"},
    "M": {"hyd": 1.9,  "vol": 162, "chg": 0, "pol": 0, "aro": 0, "name": "Met"},
    "F": {"hyd": 2.8,  "vol": 189, "chg": 0, "pol": 0, "aro": 1, "name": "Phe"},
    "P": {"hyd": -1.6, "vol": 112, "chg": 0, "pol": 0, "aro": 0, "name": "Pro"},
    "S": {"hyd": -0.8, "vol":  89, "chg": 0, "pol": 1, "aro": 0, "name": "Ser"},
    "T": {"hyd": -0.7, "vol": 116, "chg": 0, "pol": 1, "aro": 0, "name": "Thr"},
    "W": {"hyd": -0.9, "vol": 227, "chg": 0, "pol": 0, "aro": 1, "name": "Trp"},
    "Y": {"hyd": -1.3, "vol": 193, "chg": 0, "pol": 1, "aro": 1, "name": "Tyr"},
    "V": {"hyd": 4.2,  "vol": 140, "chg": 0, "pol": 0, "aro": 0, "name": "Val"},
}

# Rough ranges for normalization
RANGES = {
    "hyd": (-4.5, 4.5),
    "vol": (60, 230),
    "chg": (-1, 1),
}

AROMATICS = {"F", "W", "Y", "H"}


def get_props(aa: str) -> dict:
    aa = aa.upper()
    if aa not in AA_PROPS:
        raise ValueError(f"Unknown amino acid: {aa}")
    return AA_PROPS[aa]


def norm(value: float, key: str) -> float:
    lo, hi = RANGES[key]
    if hi == lo:
        return 0.0
    v = (value - lo) / (hi - lo)
    if v < 0:
        v = 0.0
    if v > 1:
        v = 1.0
    return v


def delta(ref: str, alt: str) -> dict:
    """Compute property deltas ref→alt and convenience flags.
    Returns dict with absolute deltas and identity flags.
    """
    pr = get_props(ref)
    pa = get_props(alt)
    d_hyd = pa["hyd"] - pr["hyd"]
    d_vol = pa["vol"] - pr["vol"]
    d_chg = pa["chg"] - pr["chg"]
    return {
        "d_hyd": d_hyd,
        "d_vol": d_vol,
        "d_chg": d_chg,
        "abs_hyd": abs(d_hyd),
        "abs_vol": abs(d_vol),
        "abs_chg": abs(d_chg),
        "norm_abs_hyd": norm(abs(d_hyd), "hyd"),
        "norm_abs_vol": norm(abs(d_vol), "vol"),
        "norm_abs_chg": norm(abs(d_chg), "chg"),
        "proline_introduced": alt.upper() == "P",
        "glycine_removed": ref.upper() == "G" and alt.upper() != "G",
        "cysteine_gain": alt.upper() == "C" and ref.upper() != "C",
        "cysteine_loss": ref.upper() == "C" and alt.upper() != "C",
        "aromatic_gain": (alt.upper() in AROMATICS) and (ref.upper() not in AROMATICS),
        "aromatic_loss": (ref.upper() in AROMATICS) and (alt.upper() not in AROMATICS),
    }

