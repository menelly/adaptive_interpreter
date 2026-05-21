#!/usr/bin/env python3
"""
GOF mechanism bank — mechanism-first, structure-driven, gene-agnostic.

Philosophy (Ren's hard lines, 2026-05-21):
  * NO hardcoded genes. Every functional label comes per-protein from UniProt
    annotation; every geometric fact comes from the AlphaFold structure.
  * Frequency and conservation are NOT mechanics. They do not appear here at
    all — they only apply a bounded nudge to the FINAL score, outside this
    module.
  * GOF is about POSITION x MECHANISM, not the magnitude of chemical change.
    Big disruptive changes tend to BREAK a protein (LOF); the subtle ones at a
    functional switch are what gain function. So we never gate GOF on Grantham
    and never zero a variant for being "chemically mild".

Each detector answers one question: "does THIS change, at THIS functional site,
trigger a known way to gain function?" A variant only needs ONE mechanism to
fire, so the overall GOF signal is a soft-OR (max) across detectors.

Mechanisms implemented:
  1. degradation_resistance   — phosphodegron / phosphosite loss (S/T/Y -> non)
  2. phosphomimetic_activation— S/T/Y -> D/E at a regulatory phosphosite
  3. dimerization_cysteine    — new surface cysteine -> constitutive disulfide
  4. autoinhibition_release   — break an interdomain salt-bridge "latch"
  5. activation_segment       — charge/structure change in a kinase catalytic
                                domain near the active machinery

A negative check (buried destabilization == LOF, not GOF) keeps the structural
genes from false positives.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

from AdaptiveInterpreter.analyzers.structure_features import get_structure_features

logger = logging.getLogger(__name__)

PHOSPHO_ACCEPTORS = set("STY")
ACIDIC = set("DE")
BASIC = set("KR")           # His treated as weak/conditional, handled separately
CHARGED = ACIDIC | BASIC

# Kyte-Doolittle hydropathy and side-chain volume (A^3) for change-effect calls.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0,
}


def _charge(aa: str) -> int:
    if aa in BASIC:
        return 1
    if aa in ACIDIC:
        return -1
    return 0


def _annotation(uniprot_id: str) -> Dict:
    """Load cached UniProt features; tolerate any failure (returns {})."""
    try:
        from AdaptiveInterpreter.data_processing.universal_protein_annotator import (
            UniversalProteinAnnotator,
        )
        return UniversalProteinAnnotator().get_uniprot_features(uniprot_id) or {}
    except Exception as e:
        logger.debug(f"annotation unavailable for {uniprot_id}: {e}")
        return {}


def _in_regions(pos: int, regions: List[Dict], needles: List[str]) -> Optional[str]:
    """Return the description of the first region covering pos whose text matches."""
    for r in regions or []:
        try:
            if r["start"] <= pos <= r["end"]:
                desc = (r.get("description") or "").lower()
                if not needles or any(n in desc for n in needles):
                    return r.get("description") or r.get("type")
        except Exception:
            continue
    return None


def _near_site(pos: int, sf, sites: List, contact_window: int = 0) -> bool:
    """True if pos IS an annotated site, or is in 3D contact with one."""
    site_positions = set()
    for s in sites or []:
        if isinstance(s, int):
            site_positions.add(s)
        elif isinstance(s, dict):
            for k in ("position", "start"):
                if k in s:
                    site_positions.add(int(s[k]))
    if pos in site_positions:
        return True
    if sf and sf.has(pos):
        nbrs = set(sf.contacts(pos))
        if site_positions & nbrs:
            return True
    return False


# ----------------------------------------------------------------------------
# Mechanism detectors. Each returns {mechanism, score, explanation} or None.
# ----------------------------------------------------------------------------

def detect_degradation_resistance(ref, alt, pos, sf, ann) -> Optional[Dict]:
    """Loss of a phosphosite/degron -> protein can't be tagged for destruction.

    A Ser/Thr phosphodegron (GSK3/CK1 etc.) is destroyed by any substitution
    that isn't another Ser/Thr — including S->Y, since tyrosine kinases don't
    rescue a Ser/Thr site (this is exactly CTNNB1 S33Y). A Tyr phosphosite is
    destroyed by any substitution away from Tyr.
    """
    if ref in ("S", "T"):
        if alt in ("S", "T"):
            return None  # Ser/Thr phosphorylatability preserved
    elif ref == "Y":
        if alt == "Y":
            return None
    else:
        return None
    phospho = ann.get("phosphorylation") or []
    annotated = pos in {int(p) if isinstance(p, int) else int(p.get("position", p.get("start", -1)))
                         for p in phospho if isinstance(p, (int, dict))}
    degron_region = _in_regions(pos, ann.get("regions"), ["degron", "destruction", "degrad"])
    if annotated:
        return {"mechanism": "degradation_resistance", "score": 0.85,
                "explanation": f"{ref}{pos}{alt} removes annotated phosphosite -> degradation resistance"}
    if degron_region:
        return {"mechanism": "degradation_resistance", "score": 0.8,
                "explanation": f"{ref}{pos}{alt} in degron region '{degron_region}' -> degradation resistance"}
    # structural proxy: phospho-acceptor in a disordered, surface, regulatory locus
    if sf and sf.has(pos) and sf.is_disordered(pos) and sf.burial(pos) == "surface":
        return {"mechanism": "degradation_resistance", "score": 0.55,
                "explanation": f"{ref}{pos}{alt} lost phospho-acceptor in disordered surface region (likely regulatory)"}
    return None


def detect_phosphomimetic(ref, alt, pos, sf, ann) -> Optional[Dict]:
    """S/T/Y -> D/E mimics constitutive phosphorylation -> constitutive activation."""
    if ref not in PHOSPHO_ACCEPTORS or alt not in ACIDIC:
        return None
    phospho = ann.get("phosphorylation") or []
    annotated = pos in {int(p) if isinstance(p, int) else int(p.get("position", p.get("start", -1)))
                        for p in phospho if isinstance(p, (int, dict))}
    if annotated:
        return {"mechanism": "phosphomimetic_activation", "score": 0.8,
                "explanation": f"{ref}{pos}{alt} phosphomimetic at annotated phosphosite -> constitutive activation"}
    if sf and sf.has(pos) and (sf.is_disordered(pos) or sf.burial(pos) == "surface"):
        return {"mechanism": "phosphomimetic_activation", "score": 0.5,
                "explanation": f"{ref}{pos}{alt} phosphomimetic in flexible/surface regulatory locus"}
    return None


def detect_dimerization_cysteine(ref, alt, pos, sf, ann) -> Optional[Dict]:
    """New surface cysteine -> unpaired thiol -> intermolecular disulfide -> constitutive dimer."""
    if alt != "C" or ref == "C":
        return None
    if not (sf and sf.has(pos)):
        return None
    if sf.burial(pos) == "buried":
        return None  # buried Cys won't reach a partner chain
    ecto = bool(ann.get("transmembrane")) or bool(ann.get("signal_peptide")) or bool(ann.get("disulfide_bonds"))
    score = 0.65 if ecto else 0.5
    where = "ecto/disulfide-rich" if ecto else "surface"
    return {"mechanism": "dimerization_cysteine", "score": score,
            "explanation": f"{ref}{pos}{alt} introduces unpaired {where} cysteine -> constitutive disulfide dimer"}


def detect_autoinhibition_release(ref, alt, pos, sf, ann) -> Optional[Dict]:
    """Break a long-range salt-bridge 'latch' that holds an autoinhibited fold."""
    if _charge(ref) == 0:
        return None
    if _charge(alt) == _charge(ref):
        return None  # charge preserved -> latch intact
    if not (sf and sf.has(pos)):
        return None
    partners = sf.salt_bridge_partners(pos)
    long_range = [j for j in partners if abs(j - pos) > 20]  # tertiary / interdomain
    if not long_range:
        return None
    score = min(0.55 + 0.1 * len(long_range), 0.75)
    return {"mechanism": "autoinhibition_release", "score": score,
            "explanation": f"{ref}{pos}{alt} breaks long-range salt-bridge latch to {long_range} -> autoinhibition release"}


def detect_activation_segment(ref, alt, pos, sf, ann) -> Optional[Dict]:
    """Charge/structure change in a kinase catalytic domain near the active machinery."""
    domains = ann.get("domains") or []
    kinase_dom = None
    for d in domains:
        desc = (d.get("description") or "").lower()
        if any(k in desc for k in ("kinase", "pi3k", "pi4k", "p110")):
            try:
                if d["start"] <= pos <= d["end"]:
                    kinase_dom = d.get("description")
                    break
            except Exception:
                continue
    if not kinase_dom:
        return None
    # the change must actually do something: introduce/flip charge, or proline, or big size jump
    dcharge = _charge(alt) - _charge(ref)
    dvol = abs(VOLUME.get(alt, 120) - VOLUME.get(ref, 120))
    meaningful = (dcharge != 0) or (alt == "P") or (ref == "P") or (dvol > 40)
    if not meaningful:
        return None
    # functional proximity: at/near an active or binding site, or a named functional region
    near = _near_site(pos, sf, ann.get("active_sites")) or _near_site(pos, sf, ann.get("binding_sites"))
    region = _in_regions(pos, ann.get("regions"),
                         ["activation", "catalytic", "g-loop", "loop", "helical", "kinase"])
    confident = sf.has(pos) and sf.plddt.get(pos, 0) >= 70 if sf else False
    if near:
        score = 0.75
        why = "adjacent to catalytic/binding machinery"
    elif region:
        score = 0.65
        why = f"in functional region '{region}'"
    elif confident:
        score = 0.55
        why = "structured catalytic-domain residue"
    else:
        return None
    chg = "charge change" if dcharge else ("proline" if "P" in (ref, alt) else "bulk change")
    return {"mechanism": "activation_segment", "score": score,
            "explanation": f"{ref}{pos}{alt} ({chg}) in {kinase_dom} domain, {why} -> constitutive activation"}


def _looks_like_lof(ref, alt, pos, sf) -> bool:
    """Buried destabilization == loss of function, not gain. Used only to veto weak GOF."""
    if not (sf and sf.has(pos)) or sf.burial(pos) != "buried":
        return False
    # proline into a helix/sheet, or large hydropathy/volume swing in the core
    if (alt == "P" and sf.ss.get(pos) in ("H", "E")):
        return True
    if abs(HYDROPATHY.get(alt, 0) - HYDROPATHY.get(ref, 0)) > 5 or \
       abs(VOLUME.get(alt, 120) - VOLUME.get(ref, 120)) > 60:
        return True
    return False


DETECTORS = [
    detect_degradation_resistance,
    detect_phosphomimetic,
    detect_dimerization_cysteine,
    detect_autoinhibition_release,
    detect_activation_segment,
]


def score_gof_mechanisms(ref: str, alt: str, pos: int, uniprot_id: str,
                         sequence: Optional[str] = None) -> Dict:
    """Run the mechanism bank. Returns gof_score (soft-OR) + per-mechanism detail.

    conservation/frequency are intentionally absent — they nudge the final
    score elsewhere, never here.
    """
    ref, alt = ref.upper(), alt.upper()
    sf = get_structure_features(uniprot_id) if uniprot_id else None
    ann = _annotation(uniprot_id) if uniprot_id else {}

    fired: List[Dict] = []
    for det in DETECTORS:
        try:
            hit = det(ref, alt, pos, sf, ann)
        except Exception as e:
            logger.debug(f"{det.__name__} failed: {e}")
            hit = None
        if hit:
            fired.append(hit)

    fired.sort(key=lambda h: h["score"], reverse=True)
    gof = fired[0]["score"] if fired else 0.0

    # LOF veto: a buried destabilizing change with only a weak GOF signal is LOF.
    if gof < 0.6 and _looks_like_lof(ref, alt, pos, sf):
        gof = min(gof, 0.1)

    return {
        "gof_score": round(gof, 3),
        "mechanisms": fired,
        "top_mechanism": fired[0]["mechanism"] if fired else None,
        "explanation": fired[0]["explanation"] if fired else "no GOF mechanism detected",
        "structure_context": sf.summary(pos) if (sf and sf.has(pos)) else None,
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)
    tests = [
        ("H", "R", 1047, "P42336", "PIK3CA H1047R (oncogenic activation, expect GOF)"),
        ("S", "Y", 33,   "P35222", "CTNNB1 S33Y (phosphodegron loss, expect GOF)"),
        ("S", "F", 37,   "P35222", "CTNNB1 S37F (phosphodegron loss, expect GOF)"),
        ("R", "H", 175,  "P04637", "TP53 R175H (buried structural, expect LOW GOF)"),
        ("G", "D", 1000, "P05997", "COL5A2 G->D triple helix (structural, expect LOW GOF)"),
    ]
    for ref, alt, pos, up, label in tests:
        r = score_gof_mechanisms(ref, alt, pos, up)
        print(f"\n{label}\n   GOF={r['gof_score']}  top={r['top_mechanism']}")
        print(f"   {r['explanation']}")
