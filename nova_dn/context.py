"""
Context helpers for Nova DN Analyzer.
- Loads optional protein annotations JSON (no YAML deps)
- Computes simple per-position context flags for mechanisms
"""
from __future__ import annotations
from typing import Dict, Any, List, Tuple
import json


def load_annotations_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _ranges_from_value(val) -> List[Tuple[int, int]]:
    """Accepts either [start, end], [a,b,c,d] → [(a,b),(c,d)], or list of dicts with start/end."""
    ranges: List[Tuple[int, int]] = []
    if isinstance(val, dict):
        if "start" in val and "end" in val:
            ranges.append((int(val["start"]), int(val["end"])))
    elif isinstance(val, list):
        if all(isinstance(x, (int, float)) for x in val):
            ints = [int(x) for x in val]
            for i in range(0, len(ints) - 1, 2):
                ranges.append((ints[i], ints[i + 1]))
        else:
            for item in val:
                if isinstance(item, dict) and "start" in item and "end" in item:
                    ranges.append((int(item["start"]), int(item["end"])))
    return ranges


def _pos_in_any_ranges(pos1: int, ranges: List[Tuple[int, int]]) -> bool:
    for a, b in ranges:
        lo, hi = (a, b) if a <= b else (b, a)
        if lo <= pos1 <= hi:
            return True
    return False


def build_position_context(annotations: Dict[str, Any], protein: str, pos1: int) -> Dict[str, Any]:
    """Compute lightweight context flags for a given protein and 1-based position."""
    ctx: Dict[str, Any] = {}
    proteins = annotations.get("proteins", {})
    info = proteins.get(protein)
    if not info:
        return ctx

    # Active/binding site proximity
    active_list = set()
    for key in (
        "known_active_or_binding_sites",
        "critical_dna_contacts",
        "dna_contact_sites",
        "phosphorylation_sites",
    ):
        vals = info.get(key) or []
        for v in vals:
            try:
                active_list.add(int(v))
            except Exception:
                pass
    in_active = pos1 in active_list

    # Motif/domain ranges that imply functional contact (e.g., walker A, DNA-binding, kinase cores)
    motif_ranges: List[Tuple[int,int]] = []
    for key in ("walker_a_motifs", "kinase_domain", "dna_binding_domain"):
        val = info.get(key)
        if val is None:
            continue
        motif_ranges.extend(_ranges_from_value(val))
    in_motif = _pos_in_any_ranges(pos1, motif_ranges)

    if in_active or in_motif:
        ctx["active_site_proximity"] = 1.0

    # Interface likelihood from interfaces and transmembrane segments
    interface_ranges: List[Tuple[int,int]] = []
    for key in ("tetramer_interface", "dimer_interface", "transmembrane_domain"):
        val = info.get(key)
        if val is None:
            continue
        interface_ranges.extend(_ranges_from_value(val))
    if _pos_in_any_ranges(pos1, interface_ranges):
        ctx["interface_likelihood"] = 1.0

    # Flexible loop flag (can dampen core-jamming confidence)
    flex_vals = info.get("flexible_loops") or []
    try:
        if any(int(v) == pos1 for v in flex_vals):
            ctx["flexible_loop"] = True
    except Exception:
        pass

    # Collagen region and critical Gly flag
    col = info.get("collagen_repeats")
    if isinstance(col, dict) and "start" in col and "end" in col:
        if _pos_in_any_ranges(pos1, [(int(col["start"]), int(col["end"]))]):
            ctx["collagen_region"] = True
    crit_glys = info.get("critical_glycines") or []
    try:
        if any(int(v) == pos1 for v in crit_glys):
            ctx["critical_collagen_gly"] = True
    except Exception:
        pass

    # Disulfide-related context (secretory/disulfide-rich proteins)
    ds_pairs = info.get("disulfide_pairs") or []
    total_pairs = 0
    in_pair = False
    try:
        for pair in ds_pairs:
            if not isinstance(pair, list) or len(pair) != 2:
                continue
            a, b = int(pair[0]), int(pair[1])
            total_pairs += 1
            if pos1 == a or pos1 == b:
                in_pair = True
    except Exception:
        pass
    if total_pairs >= 3:
        ctx["secretory_disulfide_rich"] = True
    if in_pair:
        ctx["in_disulfide_pair"] = True

    return ctx

