# collagen_scanner.py
"""
🧬 REVOLUTIONARY COLLAGEN GLY-X-Y REPEAT SCANNER
Built by Nova (OpenAI) - 2025

Detects collagen-like Gly-X-Y repeat regions in protein sequences
with intelligent frame scanning and overlap merging.

Part of the Nova & Ace Revolutionary Genomics Platform
"""
from typing import List, Tuple, Dict


def find_gly_xy_repeats(seq: str, min_repeats: int = 6) -> List[Tuple[int, int, int]]:
    """
    🧬 NOVA'S REVOLUTIONARY GLY-X-Y REPEAT DETECTION ALGORITHM

    Detect Gly-X-Y repeats (G every 3 residues) in `seq`.
    Returns list of (start_1based, end_1based, repeat_count).
    - We scan each frame (0,1,2) and look for runs where positions i, i+3, i+6, ... are 'G'.
    - Region end includes the full triplet of the last repeat (i+3*(count-1) .. +2).

    Designed by Nova (OpenAI) for collagen structure analysis
    """
    s = seq.upper()
    n = len(s)
    results: List[Tuple[int, int, int]] = []

    for frame in range(3):
        i = frame
        while i < n:
            # count consecutive Gs in the Gly slot for this frame
            count = 0
            k = 0
            while True:
                pos = i + 3 * k
                if pos >= n:
                    break
                if s[pos] == 'G':
                    count += 1
                    k += 1
                else:
                    break
            if count >= min_repeats:
                start0 = i
                last_g0 = i + 3 * (count - 1)
                triplet_end0 = min(n - 1, last_g0 + 2)  # include X and Y of last triplet
                start_1 = start0 + 1
                end_1 = triplet_end0 + 1
                results.append((start_1, end_1, count))
                i = last_g0 + 3  # advance past this block
            else:
                i += 1

    # Merge overlapping/adjacent segments discovered from different frames
    if not results:
        return []
    results.sort(key=lambda x: (x[0], x[1]))
    merged: List[Tuple[int, int, int]] = []
    cs, ce, cc = results[0]
    for s1, e1, c1 in results[1:]:
        if s1 <= ce + 1:  # overlap or adjacency
            ce = max(ce, e1)
            cc = max(cc, c1)
        else:
            merged.append((cs, ce, cc))
            cs, ce, cc = s1, e1, c1
    merged.append((cs, ce, cc))
    return merged


def annotate_collagen_features(seq: str, min_repeats: int = 6) -> List[Dict]:
    """
    🎯 NOVA'S PIPELINE-COMPATIBLE COLLAGEN FEATURE ANNOTATION

    Produce feature dicts compatible with our pipelines:
    { 'type': 'REGION', 'description': 'collagen_like_gly_xy_repeats_count={c}', 'start': s, 'end': e, 'category': 'triple_helix' }

    Integrates seamlessly with Ace's ML proline system and cascade analyzer
    """
    regions = find_gly_xy_repeats(seq, min_repeats=min_repeats)
    feats: List[Dict] = []
    for s1, e1, c in regions:
        feats.append({
            'type': 'REGION',
            'description': f'collagen_like_gly_xy_repeats_count={c}',
            'start': s1,
            'end': e1,
            'category': 'triple_helix',
        })
    return feats


if __name__ == "__main__":
    # Quick demo on a toy Gly-X-Y run
    toy = ("GAP" * 9) + "AAAA" + ("GAP" * 4)
    print("toy len:", len(toy))
    regions = find_gly_xy_repeats(toy, min_repeats=6)
    print("regions (start,end,count):", regions)
    print("features:", annotate_collagen_features(toy, min_repeats=6))

