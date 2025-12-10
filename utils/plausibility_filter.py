"""
plausibility_filter.py

Biological Pathogenic Mechanism Plausibility Framework
------------------------------------------------------
Implements domain-aware filtering of mechanism scores
to reflect biological plausibility across gene families.

Pipeline steps:
1. Collect raw analyzer scores
2. Classify gene into a family (enzymes, channels, etc.)
3. Apply plausibility weights based on family rules
4. Preserve provenance + rationales for explainability
5. Pass filtered scores downstream to synergy scoring
6. Final classification incorporates plausibility + synergy

Author: Nova (with Ren + Ace)
"""

from typing import Dict, Any, List
import json
import logging
from pathlib import Path
from collections import defaultdict


# ======================================================
# STEP 3: GENE FAMILY RULES (tiered plausibility weights)
# ======================================================

PATHOGENICITY_RULES: Dict[str, Dict[str, float]] = {
    "ENZYME": {"LOF": 1.0, "DN": 1.0, "GOF": 0.0},
    "ION_CHANNEL": {"LOF": 1.0, "DN": 1.0, "GOF": 1.0},
    "STRUCTURAL": {"LOF": 1.0, "DN": 1.0, "GOF": 0.25},
    # --- Collagen subclasses (BALANCED for AD inheritance) ---
    "COLLAGEN_FIBRILLAR": {"LOF": 0.85, "DN": 1.25, "GOF": 0.02},  # Balanced: DN favored but LOF viable
    "COLLAGEN_NETWORK": {"LOF": 0.85, "DN": 1.25, "GOF": 0.02},    # Same balance
    "COLLAGEN_ANCHORING": {"LOF": 0.9, "DN": 1.3, "GOF": 0.0},     # Slightly more LOF tolerance
    "COLLAGEN_FACIT": {"LOF": 0.9, "DN": 1.15, "GOF": 0.05},       # Balanced
    # --- Newly added families ---
    "INTERMEDIATE_FILAMENT": {"LOF": 0.8, "DN": 1.35, "GOF": 0.0},
    "FIBRILLIN": {"LOF": 0.85, "DN": 1.2, "GOF": 0.05},  # Balanced for Marfan (AD with both LOF+DN)
    "ELASTIN": {"LOF": 1.2, "DN": 0.9, "GOF": 0.0},
    "CYTOSKELETON_POLYMER": {"LOF": 0.9, "DN": 1.25, "GOF": 0.1},
    "LAMIN": {"LOF": 1.0, "DN": 1.1, "GOF": 0.1},
    "RTK_MAPK": {"LOF": 0.6, "DN": 0.8, "GOF": 1.3},
    "NEGATIVE_REGULATOR": {"LOF": 1.2, "DN": 0.9, "GOF": 0.0},
    # ----------------------------
    "TRANSCRIPTION_FACTOR": {"LOF": 1.0, "DN": 1.0, "GOF": 1.0},
    "TUMOR_SUPPRESSOR": {"LOF": 1.0, "DN": 1.0, "GOF": 0.05},  # GOF is NOT pathogenic for tumor suppressors!
    "ONCOGENE": {"LOF": 1.0, "DN": 1.0, "GOF": 1.0},  # All mechanisms possible for oncogenes
    "MUSCULAR_DYSTROPHY": {"LOF": 1.0, "DN": 0.5, "GOF": 0.0},  # Muscular dystrophy genes
    "RIBOSOMAL_PROTEIN": {"LOF": 1.0, "DN": 0.5, "GOF": 0.0},  # More ribosomes = good for growth
    "MOTOR_PROTEIN": {"LOF": 1.0, "DN": 1.0, "GOF": 0.25},  # Myosins: GOF rarely pathogenic
    "DNA_REPAIR": {"LOF": 1.0, "DN": 0.5, "GOF": 0.0},  # DNA repair: more repair = good
    "GENERAL": {"LOF": 1.0, "DN": 0.5, "GOF": 0.5},  # fallback
}


# ======================================================
# STEP 5: WEIGHTED GENE FAMILY CLASSIFICATION
# ======================================================

# Load category keywords from JSON configuration
def load_category_config():
    """Load category keywords and metadata from JSON file"""
    config_path = Path(__file__).parent.parent / "config" / "category_keywords.json"
    try:
        with open(config_path, "r") as f:
            config = json.load(f)
            return config.get("categories", {}), config.get("meta", {})
    except FileNotFoundError:
        logging.warning(f"Category keywords file not found at {config_path}, using fallback")
        return {}, {}

# Global category keywords and metadata (loaded once)
CATEGORY_KEYWORDS, CATEGORY_META = load_category_config()
PRIORITY_ORDER = CATEGORY_META.get("priority_order", [])
DEFAULT_THRESHOLD = CATEGORY_META.get("threshold_margin", 0.2)

def classify_protein_weighted(function_text: str, threshold_margin: float = None):
    """
    Weighted keyword-based classification of protein families.
    Returns top category + full score breakdown.

    Based on Nova's brilliant weighted scoring approach!
    """
    if threshold_margin is None:
        threshold_margin = DEFAULT_THRESHOLD

    function_lower = function_text.lower()
    scores = defaultdict(float)

    # Score each category by keyword hits
    for category, keywords in CATEGORY_KEYWORDS.items():
        for kw, weight in keywords.items():
            if kw in function_lower:
                scores[category] += weight

    if not scores:
        return "UNCLASSIFIED", {}

    # Find the best category
    sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    top_category, top_score = sorted_scores[0]

    # Check for close competitors (MULTIROLE detection)
    if len(sorted_scores) > 1:
        second_category, second_score = sorted_scores[1]
        if (top_score - second_score) / top_score < threshold_margin:
            # Use priority order for tie-breaking if available
            if PRIORITY_ORDER:
                try:
                    top_priority = PRIORITY_ORDER.index(top_category)
                    second_priority = PRIORITY_ORDER.index(second_category)
                    if top_priority < second_priority:  # Lower index = higher priority
                        return top_category, dict(sorted_scores)
                    else:
                        return second_category, dict(sorted_scores)
                except ValueError:
                    # If categories not in priority list, use MULTIROLE
                    return f"MULTIROLE ({top_category}, {second_category})", dict(sorted_scores)
            else:
                return f"MULTIROLE ({top_category}, {second_category})", dict(sorted_scores)

    # Otherwise return single top category
    return top_category, dict(sorted_scores)

def classify_gene_family(gene_symbol: str, uniprot_function: str, go_terms: List[str],
                        override_family: str = None) -> str:
    """
    Classify gene into pathogenicity-relevant families.
    Uses Nova's weighted scoring approach with biological intelligence!

    🧬 REVOLUTIONARY WEIGHTED CLASSIFICATION SYSTEM!

    Args:
        gene_symbol: Gene symbol (e.g., 'MYO7A')
        uniprot_function: UniProt function description
        go_terms: List of GO terms
        override_family: Optional manual override (e.g., 'ONCOGENE' for KRAS)
    """

    # 🎯 GENETICIST OVERRIDE - Trust the expert!
    if override_family:
        logging.info(f"🔧 Manual override: {gene_symbol} → {override_family} (geneticist specified)")
        return override_family.upper()
    terms = [t.lower() for t in go_terms]
    function_lower = uniprot_function.lower()

    # Combine UniProt function and GO terms for comprehensive analysis
    combined_text = f"{uniprot_function} {' '.join(go_terms)}"

    # Use Nova's weighted scoring system
    classification, score_breakdown = classify_protein_weighted(combined_text)

    # Handle MULTIROLE classifications
    if classification.startswith("MULTIROLE"):
        # For now, extract the primary classification from MULTIROLE
        # Format: "MULTIROLE (PRIMARY, SECONDARY)"
        primary = classification.split("(")[1].split(",")[0].strip()
        logging.info(f"🔥 MULTIROLE protein detected: {classification}")
        logging.info(f"📊 Score breakdown: {score_breakdown}")
        return primary

    # Handle UNCLASSIFIED - fall back to legacy patterns for special cases
    if classification == "UNCLASSIFIED":
        # 🧬 CHECK SPECIFIC DISEASE FAMILIES FIRST (before generic patterns)

        # Muscular dystrophy genes
        muscular_dystrophy_keywords = [
            "muscular dystrophy", "muscle dystrophy", "limb-girdle", "dysferlin",
            "dystrophin", "sarcoglycan", "muscle membrane", "dystroglycan", "fkrp"
        ]
        if any(keyword in function_lower for keyword in muscular_dystrophy_keywords) or gene_symbol == "FKRP":
            return "MUSCULAR_DYSTROPHY"

        # 🧬 SMART AD STRUCTURAL DETECTION
        ad_structural_patterns = [
            "marfan", "fibrillin", "connective tissue",
            "osteogenesis imperfecta", "brittle bone", "collagen",
            "ehlers-danlos", "eds", "hypermobility",
            "aortic dilatation", "aortic aneurysm",
            "dominant negative", "poison subunit", "haploinsufficiency",
            "autosomal dominant", "dominant disorder"
        ]

        if any(pattern in function_lower for pattern in ad_structural_patterns):
            # Use gene symbol for collagen subtype detection
            if "collagen" in function_lower or gene_symbol.startswith("COL"):
                if any(x in function_lower for x in ["type i", "type ii", "type iii", "type v", "type xi"]):
                    return "COLLAGEN_FIBRILLAR"
                elif "type iv" in function_lower:
                    return "COLLAGEN_NETWORK"
                elif "type vii" in function_lower or gene_symbol == "COL7A1":
                    return "COLLAGEN_ANCHORING"
                elif any(x in function_lower for x in ["type ix", "type xii", "type xiv"]):
                    return "COLLAGEN_FACIT"
                else:
                    return "STRUCTURAL"
            elif "fibrillin" in function_lower:
                return "FIBRILLIN"
            else:
                return "STRUCTURAL"

        # Additional legacy patterns for special cases

        ribosomal_keywords = ["ribosomal protein", "ribosome", "rpl", "rps", "60s ribosomal", "40s ribosomal"]
        if any(keyword in function_lower for keyword in ribosomal_keywords):
            return "RIBOSOMAL_PROTEIN"

        motor_keywords = ["myosin", "kinesin", "dynein", "motor protein", "actin-based motor", "microtubule motor"]
        if any(keyword in function_lower for keyword in motor_keywords):
            return "MOTOR_PROTEIN"

        # Fall back to GO term analysis
        if any("channel" in term or "transport" in term for term in terms):
            return "ION_CHANNEL"
        elif any(term in terms for term in ["structural", "cytoskeleton", "basement membrane", "extracellular matrix"]):
            return "STRUCTURAL"
        elif any(term in terms for term in ["transcription", "methyltransferase", "dna-binding"]):
            return "TRANSCRIPTION_FACTOR"
        else:
            return "GENERAL"

    # Log the successful weighted classification
    logging.info(f"🎯 Weighted classification: {classification}")
    logging.info(f"📊 Score breakdown: {score_breakdown}")

    return classification


# ======================================================
# STEP 2 + 4: PLAUSIBILITY FILTER
# ======================================================

def _dn_evidence_score(gene_symbol: str, uniprot_function: str = "", go_terms: List[str] | None = None) -> tuple[float, list[str]]:
    """Heuristic dominant-negative evidence score from local text signals.
    - Searches UniProt function text for phrases like 'dominant negative'.
    - Optionally looks at GO terms (rare direct hits expected).
    Returns (score, hits).
    """
    if go_terms is None:
        go_terms = []
    score = 0.0
    hits: list[str] = []
    func = (uniprot_function or "").lower()
    # Core phrases
    patterns = [
        "dominant negative",
        "dominant-negative",
        "exerts a dominant negative",
        "acts in a dominant negative",
    ]
    for p in patterns:
        if p in func:
            score += 1.0
            hits.append(f"uniprot:{p}")
            break
    # GO rarely encodes this explicitly; keep very conservative
    for term in go_terms:
        t = (term or "").lower()
        if "dominant-negative" in t or "dominant negative" in t:
            score += 0.5
            hits.append("go:dominant-negative")
            break
    return score, hits


def _adjust_dn_weight(base_weight: float, evidence_score: float) -> float:
    """Map evidence score → adjusted DN weight conservatively.
    - strong (>=1.0): at least 1.0
    - moderate (>=0.5): at least 0.9
    - else: unchanged
    """
    if evidence_score >= 1.0:
        return max(base_weight, 1.0)
    if evidence_score >= 0.5:
        return max(base_weight, 0.9)
    return base_weight


def apply_pathogenicity_filter(
    raw_scores: Dict[str, float],
    gene_family: str,
    gene_symbol: str = "",
    uniprot_function: str = "",
    go_terms: List[str] | None = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Apply gene-family-specific plausibility rules to raw mechanism scores.

    🧬 REN'S ELEGANT REFACTOR (December 2025):
    Now uses ADVISORY mode instead of hard filtering!
    - Raw scores are preserved for classification (trust the analyzer's work)
    - Atypical mechanisms are FLAGGED, not squashed
    - "The math says this is a problem. THAT SAID, this family isn't usually <mechanism>"

    Args:
        raw_scores: dict of {mechanism: score}
        gene_family: str, e.g. "ENZYME", "ION_CHANNEL"
        gene_symbol/uniprot_function/go_terms: used for DN evidence-based adjustment

    Returns:
        dict with mechanism → {raw_score, weighted_score, status, rationale, atypical_flag, ...}
    """
    if go_terms is None:
        go_terms = []
    rules = PATHOGENICITY_RULES.get(gene_family, PATHOGENICITY_RULES["GENERAL"])
    filtered: Dict[str, Dict[str, Any]] = {}

    # Compute DN evidence once per gene
    dn_ev_score, dn_ev_hits = _dn_evidence_score(gene_symbol, uniprot_function, go_terms)

    for mech, score in raw_scores.items():
        base_weight = rules.get(mech, 0.0)
        weight = base_weight
        evidence_applied = False
        if mech == "DN":
            new_weight = _adjust_dn_weight(base_weight, dn_ev_score)
            if new_weight != base_weight:
                weight = new_weight
                evidence_applied = True
        weighted = score * weight

        if weight == 1.0:
            status = "kept"
        elif weight in (0.5, 0.25, 0.8, 0.9):
            status = "downweighted"
        elif weight == 0.0:
            status = "zeroed"
        else:
            status = "weighted"

        rationale = f"{mech} {status} for {gene_family} (weight {weight})"

        # 🧬 ATYPICAL MECHANISM FLAG (Ren's advisory mode)
        # If the weight < 1.0 AND the raw score is significant (>0.4), flag for review
        atypical_flag = None
        if weight < 1.0 and score >= 0.4:
            atypical_flag = f"⚠️ {mech} is atypical for {gene_family} (family weight: {weight}) - verify mechanism"

        entry = {
            "raw_score": score,
            "weighted_score": weighted,
            "status": status,
            "rationale": rationale,
            "atypical_flag": atypical_flag,  # 🧬 NEW: Advisory flag instead of hard filter
            "family_weight": weight,  # 🧬 NEW: Expose the weight for transparency
        }
        if mech == "DN":
            entry.update({
                "dn_weight_base": base_weight,
                "dn_weight_final": weight,
                "dn_evidence_score": dn_ev_score,
                "dn_evidence_hits": dn_ev_hits,
                "dn_evidence_applied": evidence_applied,
            })

        filtered[mech] = entry

    return filtered


# ======================================================
# STEP 6: PIPELINE FUNCTION
# ======================================================

def plausibility_pipeline(
    gene_symbol: str,
    raw_scores: Dict[str, float],
    uniprot_function: str = "",
    go_terms: List[str] = None,
    override_family: str = None,
    use_advisory_mode: bool = True,  # 🧬 NEW: Ren's elegant refactor!
) -> Dict[str, Any]:
    """
    Full plausibility pipeline: classify → filter → return explainable scores.

    🧬 REN'S ELEGANT REFACTOR (December 2025):
    Now supports ADVISORY MODE (default: True)!
    - advisory_mode=True: Use RAW scores for classification, flag atypical mechanisms
    - advisory_mode=False: Use FILTERED scores (legacy behavior)

    Args:
        gene_symbol: str, e.g. "DLD"
        raw_scores: dict of {mechanism: score}
        uniprot_function: UniProt function description
        go_terms: list of GO terms (optional)
        override_family: Optional manual family override (e.g., 'ONCOGENE' for KRAS)
        use_advisory_mode: If True (default), use raw scores + advisory flags instead of hard filtering

    Returns:
        dict with:
            - gene_family
            - filtered_scores (per mechanism with rationale + atypical_flag)
            - final_scores (RAW values if advisory_mode, else weighted)
            - atypical_mechanisms (list of mechanisms flagged as atypical for this family)
            - advisory_mode (whether advisory mode was used)
    """
    if go_terms is None:
        go_terms = []

    # Step 1: raw_scores are provided
    # Step 2: classify gene family (with optional override)
    gene_family = classify_gene_family(gene_symbol, uniprot_function, go_terms, override_family)

    # Step 3-4: apply plausibility filter (now includes atypical_flag)
    filtered = apply_pathogenicity_filter(
        raw_scores,
        gene_family,
        gene_symbol=gene_symbol,
        uniprot_function=uniprot_function,
        go_terms=go_terms,
    )

    # 🧬 REN'S ADVISORY MODE: Trust the analyzer, flag atypical mechanisms
    if use_advisory_mode:
        # Use RAW scores for classification (trust the analyzer's work)
        final_scores = {mech: vals["raw_score"] for mech, vals in filtered.items()}
    else:
        # Legacy: use weighted scores
        final_scores = {mech: vals["weighted_score"] for mech, vals in filtered.items()}

    # Collect atypical mechanism flags
    atypical_mechanisms = []
    for mech, vals in filtered.items():
        if vals.get("atypical_flag"):
            atypical_mechanisms.append({
                "mechanism": mech,
                "flag": vals["atypical_flag"],
                "raw_score": vals["raw_score"],
                "family_weight": vals.get("family_weight", 1.0),
            })

    return {
        "gene_symbol": gene_symbol,
        "gene_family": gene_family,
        "filtered_scores": filtered,
        "final_scores": final_scores,
        "atypical_mechanisms": atypical_mechanisms,  # 🧬 NEW: Advisory flags
        "advisory_mode": use_advisory_mode,  # 🧬 NEW: Track which mode was used
    }


# ======================================================
# Example usage
# ======================================================
if __name__ == "__main__":
    raw = {"DN": 0.45, "LOF": 0.62, "GOF": 0.84}
    result = plausibility_pipeline("DLD", raw, "pyruvate dehydrogenase enzyme", ["enzyme", "catalytic"])
    from pprint import pprint
    pprint(result)
