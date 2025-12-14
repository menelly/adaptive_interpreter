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
    Returns dict with absolute deltas, identity flags, and Grantham distance.
    """
    pr = get_props(ref)
    pa = get_props(alt)
    d_hyd = pa["hyd"] - pr["hyd"]
    d_vol = pa["vol"] - pr["vol"]
    d_chg = pa["chg"] - pr["chg"]

    # Get Grantham distance (lazy import to avoid circular dependency)
    grantham = _get_grantham_distance_internal(ref, alt)

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
        # 🧬 NEW: Grantham distance for substitution severity
        "grantham": grantham,
        "grantham_conservative": grantham < 50,  # Very similar AAs
        "grantham_moderate": 50 <= grantham < 100,
        "grantham_radical": grantham >= 100,  # Very different AAs
    }


# Internal helper to avoid forward reference issues
def _get_grantham_distance_internal(ref: str, alt: str) -> int:
    """Internal Grantham lookup (matrix defined later in file)."""
    ref, alt = ref.upper(), alt.upper()
    if ref == alt:
        return 0
    # Matrix defined at end of file, use globals
    if 'GRANTHAM_MATRIX' not in globals():
        return 100  # Fallback before matrix is defined
    dist = GRANTHAM_MATRIX.get((ref, alt))
    if dist is None:
        dist = GRANTHAM_MATRIX.get((alt, ref))
    return dist if dist is not None else 100


# Helper functions for Nova's Universal Lattice Disruption Framework

def get_amino_acid_volume_change(ref: str, alt: str) -> float:
    """Get volume change in Å³ (positive = expansion, negative = contraction)."""
    ref_props = get_props(ref)
    alt_props = get_props(alt)
    return alt_props["vol"] - ref_props["vol"]


def get_charge_change(ref: str, alt: str) -> int:
    """Get charge change (positive = more positive, negative = more negative)."""
    ref_props = get_props(ref)
    alt_props = get_props(alt)
    return alt_props["chg"] - ref_props["chg"]


def get_hydrophobicity_change(ref: str, alt: str) -> float:
    """Get hydrophobicity change (positive = more hydrophobic)."""
    ref_props = get_props(ref)
    alt_props = get_props(alt)
    return alt_props["hyd"] - ref_props["hyd"]


# 🧬 GRANTHAM DISTANCE MATRIX
# Grantham R. (1974) Amino acid difference formula to help explain protein evolution.
# Higher values = more radical substitution
GRANTHAM_MATRIX = {
    ('A', 'R'): 112, ('A', 'N'): 111, ('A', 'D'): 126, ('A', 'C'): 195,
    ('A', 'Q'): 91, ('A', 'E'): 107, ('A', 'G'): 60, ('A', 'H'): 86, ('A', 'I'): 94,
    ('A', 'L'): 96, ('A', 'K'): 106, ('A', 'M'): 84, ('A', 'F'): 113, ('A', 'P'): 27,
    ('A', 'S'): 99, ('A', 'T'): 58, ('A', 'W'): 148, ('A', 'Y'): 112, ('A', 'V'): 64,
    ('R', 'N'): 86, ('R', 'D'): 96, ('R', 'C'): 180, ('R', 'Q'): 43, ('R', 'E'): 54,
    ('R', 'G'): 125, ('R', 'H'): 29, ('R', 'I'): 97, ('R', 'L'): 102, ('R', 'K'): 26,
    ('R', 'M'): 91, ('R', 'F'): 97, ('R', 'P'): 103, ('R', 'S'): 110, ('R', 'T'): 71,
    ('R', 'W'): 101, ('R', 'Y'): 77, ('R', 'V'): 96,
    ('N', 'D'): 23, ('N', 'C'): 139, ('N', 'Q'): 46, ('N', 'E'): 42, ('N', 'G'): 80,
    ('N', 'H'): 68, ('N', 'I'): 149, ('N', 'L'): 153, ('N', 'K'): 94, ('N', 'M'): 142,
    ('N', 'F'): 158, ('N', 'P'): 91, ('N', 'S'): 46, ('N', 'T'): 65, ('N', 'W'): 174,
    ('N', 'Y'): 143, ('N', 'V'): 133,
    ('D', 'C'): 154, ('D', 'Q'): 61, ('D', 'E'): 45, ('D', 'G'): 94, ('D', 'H'): 81,
    ('D', 'I'): 168, ('D', 'L'): 172, ('D', 'K'): 101, ('D', 'M'): 160, ('D', 'F'): 177,
    ('D', 'P'): 108, ('D', 'S'): 65, ('D', 'T'): 85, ('D', 'W'): 181, ('D', 'Y'): 160,
    ('D', 'V'): 152,
    ('C', 'Q'): 154, ('C', 'E'): 170, ('C', 'G'): 159, ('C', 'H'): 174, ('C', 'I'): 198,
    ('C', 'L'): 198, ('C', 'K'): 202, ('C', 'M'): 196, ('C', 'F'): 205, ('C', 'P'): 169,
    ('C', 'S'): 112, ('C', 'T'): 149, ('C', 'W'): 215, ('C', 'Y'): 194, ('C', 'V'): 192,
    ('Q', 'E'): 29, ('Q', 'G'): 87, ('Q', 'H'): 24, ('Q', 'I'): 109, ('Q', 'L'): 113,
    ('Q', 'K'): 53, ('Q', 'M'): 101, ('Q', 'F'): 116, ('Q', 'P'): 76, ('Q', 'S'): 68,
    ('Q', 'T'): 42, ('Q', 'W'): 130, ('Q', 'Y'): 99, ('Q', 'V'): 96,
    ('E', 'G'): 98, ('E', 'H'): 40, ('E', 'I'): 134, ('E', 'L'): 138, ('E', 'K'): 56,
    ('E', 'M'): 126, ('E', 'F'): 140, ('E', 'P'): 93, ('E', 'S'): 80, ('E', 'T'): 65,
    ('E', 'W'): 152, ('E', 'Y'): 122, ('E', 'V'): 121,
    ('G', 'H'): 98, ('G', 'I'): 135, ('G', 'L'): 138, ('G', 'K'): 127, ('G', 'M'): 127,
    ('G', 'F'): 153, ('G', 'P'): 42, ('G', 'S'): 56, ('G', 'T'): 59, ('G', 'W'): 184,
    ('G', 'Y'): 147, ('G', 'V'): 109,
    ('H', 'I'): 94, ('H', 'L'): 99, ('H', 'K'): 32, ('H', 'M'): 87, ('H', 'F'): 100,
    ('H', 'P'): 77, ('H', 'S'): 89, ('H', 'T'): 47, ('H', 'W'): 115, ('H', 'Y'): 83,
    ('H', 'V'): 84,
    ('I', 'L'): 5, ('I', 'K'): 102, ('I', 'M'): 10, ('I', 'F'): 21, ('I', 'P'): 95,
    ('I', 'S'): 142, ('I', 'T'): 89, ('I', 'W'): 61, ('I', 'Y'): 33, ('I', 'V'): 29,
    ('L', 'K'): 107, ('L', 'M'): 15, ('L', 'F'): 22, ('L', 'P'): 98, ('L', 'S'): 145,
    ('L', 'T'): 92, ('L', 'W'): 61, ('L', 'Y'): 36, ('L', 'V'): 32,
    ('K', 'M'): 95, ('K', 'F'): 102, ('K', 'P'): 103, ('K', 'S'): 121, ('K', 'T'): 78,
    ('K', 'W'): 110, ('K', 'Y'): 85, ('K', 'V'): 97,
    ('M', 'F'): 28, ('M', 'P'): 87, ('M', 'S'): 135, ('M', 'T'): 81, ('M', 'W'): 67,
    ('M', 'Y'): 36, ('M', 'V'): 21,
    ('F', 'P'): 114, ('F', 'S'): 155, ('F', 'T'): 103, ('F', 'W'): 40, ('F', 'Y'): 22,
    ('F', 'V'): 50,
    ('P', 'S'): 74, ('P', 'T'): 38, ('P', 'W'): 147, ('P', 'Y'): 110, ('P', 'V'): 68,
    ('S', 'T'): 58, ('S', 'W'): 177, ('S', 'Y'): 144, ('S', 'V'): 124,
    ('T', 'W'): 128, ('T', 'Y'): 92, ('T', 'V'): 69,
    ('W', 'Y'): 37, ('W', 'V'): 88,
    ('Y', 'V'): 55,
}


def get_grantham_distance(ref: str, alt: str) -> int:
    """Get Grantham distance between two amino acids.

    Returns:
        int: Grantham distance (0-215 scale)
             0 = identical, <50 = conservative, 50-100 = moderate, >100 = radical
    """
    ref, alt = ref.upper(), alt.upper()
    if ref == alt:
        return 0

    # Try both orientations
    dist = GRANTHAM_MATRIX.get((ref, alt))
    if dist is None:
        dist = GRANTHAM_MATRIX.get((alt, ref))

    # Fallback for unknown pairs
    return dist if dist is not None else 100


def get_grantham_severity(ref: str, alt: str) -> str:
    """Classify substitution severity based on Grantham distance.

    Returns:
        str: 'conservative', 'moderate', or 'radical'
    """
    dist = get_grantham_distance(ref, alt)
    if dist < 50:
        return 'conservative'
    elif dist < 100:
        return 'moderate'
    else:
        return 'radical'

