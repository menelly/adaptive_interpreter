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
    "TUMOR_SUPPRESSOR": {"LOF": 1.0, "DN": 1.0, "GOF": 0.25},
    "ONCOGENE": {"LOF": 1.0, "DN": 1.0, "GOF": 1.0},  # All mechanisms possible for oncogenes
    "AUTOSOMAL_RECESSIVE": {"LOF": 1.3, "DN": 0.4, "GOF": 0.0},  # AR diseases: STRONG LOF boost - designed to lose function!
    "MUSCULAR_DYSTROPHY": {"LOF": 1.0, "DN": 0.5, "GOF": 0.0},  # MD: Let GO terms handle AD vs AR distinction
    "RIBOSOMAL_PROTEIN": {"LOF": 1.0, "DN": 0.5, "GOF": 0.0},  # More ribosomes = good for growth
    "MOTOR_PROTEIN": {"LOF": 1.0, "DN": 1.0, "GOF": 0.25},  # Myosins: GOF rarely pathogenic
    "DNA_REPAIR": {"LOF": 1.0, "DN": 0.5, "GOF": 0.0},  # DNA repair: more repair = good
    "GENERAL": {"LOF": 1.0, "DN": 0.5, "GOF": 0.5},  # fallback
}


# ======================================================
# STEP 5: CLASSIFY GENE FAMILY
# ======================================================

def classify_gene_family(gene_symbol: str, uniprot_function: str, go_terms: List[str]) -> str:
    """
    Classify gene into pathogenicity-relevant families.
    Uses GO terms + UniProt function + keyword parsing.

    🧬 SMART INHERITANCE PATTERN DETECTION!
    """
    terms = [t.lower() for t in go_terms]
    function_lower = uniprot_function.lower()

    # 🧬 SMART AR DETECTION (based on biological clues)
    ar_disease_patterns = [
        # Classic AR diseases
        "cystic fibrosis", "cf", "cftr",
        "stargardt", "abca4", "retinal dystrophy",
        "muscular dystrophy", "limb-girdle", "fkrp", "dysferlin",
        "sickle cell", "thalassemia", "hemoglobin",
        "phenylketonuria", "pku", "phenylalanine hydroxylase",
        "galactosemia", "galt", "galactose-1-phosphate",
        "glycogen storage", "pompe", "gaucher",
        # AR functional patterns
        "enzyme deficiency", "metabolic disorder", "inborn error",
        "recessive disorder", "autosomal recessive", "compound heterozygous"
    ]

    # Check if this looks like an AR gene
    if any(pattern in function_lower for pattern in ar_disease_patterns):
        return "AUTOSOMAL_RECESSIVE"

    # 🧬 SMART AD STRUCTURAL DETECTION
    ad_structural_patterns = [
        # Classic AD structural diseases
        "marfan", "fibrillin", "connective tissue",
        "osteogenesis imperfecta", "brittle bone", "collagen",
        "ehlers-danlos", "eds", "hypermobility",
        "aortic dilatation", "aortic aneurysm",
        # AD functional patterns
        "dominant negative", "poison subunit", "haploinsufficiency",
        "autosomal dominant", "dominant disorder"
    ]

    if any(pattern in function_lower for pattern in ad_structural_patterns):
        # Further classify structural proteins
        if "collagen" in function_lower:
            return "COLLAGEN_FIBRILLAR"
        elif "fibrillin" in function_lower:
            return "FIBRILLIN"
        else:
            return "STRUCTURAL"

    # Check for oncogene keywords first (before general kinase check)
    # 🔥 FIXED: More specific keywords to avoid false positives like FBN1
    oncogene_keywords = [
        "cell growth", "cell proliferation", "mitogenic signals",
        "growth factor receptor", "tyrosine kinase", "tyrosine-protein kinase",
        "receptor tyrosine kinase", "proto-oncogene", "oncogene",
        "fibroblast growth factor receptor", "epidermal growth factor receptor",
        "phosphoinositide-3-kinase", "pi3k", "ras protein",
        "growth factor signaling", "mitogen-activated protein kinase"
    ]
    # 🚫 REMOVED: "growth factor" (too broad - matches FBN1 which BINDS growth factors)
    # 🚫 REMOVED: "fibroblast growth factor" (too broad)
    # 🚫 REMOVED: "epidermal growth factor" (too broad)
    if any(keyword in function_lower for keyword in oncogene_keywords):
        return "ONCOGENE"

    # Check for DNA repair keywords FIRST (before tumor suppressor to avoid conflict)
    dna_repair_keywords = ["dna repair", "dna damage", "homologous recombination", "double-strand break", "mismatch repair"]
    if any(keyword in function_lower for keyword in dna_repair_keywords):
        return "DNA_REPAIR"

    # Check for tumor suppressor keywords (after DNA repair check)
    tumor_suppressor_keywords = [
        "tumor suppressor", "cell cycle checkpoint",
        "apoptosis regulation", "growth inhibition"
    ]
    if any(keyword in function_lower for keyword in tumor_suppressor_keywords):
        return "TUMOR_SUPPRESSOR"

    # Check for muscular dystrophy keywords (🔥 FIXED: removed "muscle fiber" - too broad!)
    muscular_dystrophy_keywords = [
        "muscular dystrophy", "muscle dystrophy", "limb-girdle", "dysferlin",
        "dystrophin", "sarcoglycan", "muscle membrane"
    ]
    if any(keyword in function_lower for keyword in muscular_dystrophy_keywords):
        return "MUSCULAR_DYSTROPHY"

    # Check for autosomal recessive disease patterns
    ar_keywords = [
        "autosomal recessive", "recessive disorder", "recessive disease",
        "homozygous", "compound heterozygous"
    ]
    if any(keyword in function_lower for keyword in ar_keywords):
        return "AUTOSOMAL_RECESSIVE"

    # Check for ribosomal protein keywords
    ribosomal_keywords = ["ribosomal protein", "ribosome", "rpl", "rps", "60s ribosomal", "40s ribosomal"]
    if any(keyword in function_lower for keyword in ribosomal_keywords):
        return "RIBOSOMAL_PROTEIN"

    # Check for motor protein keywords
    motor_keywords = ["myosin", "kinesin", "dynein", "motor protein", "actin-based motor", "microtubule motor"]
    if any(keyword in function_lower for keyword in motor_keywords):
        return "MOTOR_PROTEIN"

    # Check for enzyme keywords (after oncogene check to avoid kinase conflicts)
    enzyme_keywords = ["dehydrogenase", "synthase", "synthetase", "transferase", "hydrolase", "isomerase", "ligase", "catalytic", "phosphatase", "reductase", "oxidase"]
    if any(keyword in function_lower for keyword in enzyme_keywords):
        return "ENZYME"

    # --- Collagen overrides ---
    is_collagen = ("collagen" in function_lower or any("collagen" in t for t in terms) or
                   gene_symbol.startswith("COL") or "anchoring fibril" in function_lower)

    if is_collagen:
        if any(x in function_lower for x in ["type i", "type ii", "type iii", "type v", "type xi"]):
            return "COLLAGEN_FIBRILLAR"
        elif "type iv" in function_lower:
            return "COLLAGEN_NETWORK"
        elif "type vii" in function_lower or gene_symbol == "COL7A1" or "anchoring fibril" in function_lower:
            return "COLLAGEN_ANCHORING"
        elif any(x in function_lower for x in ["type ix", "type xii", "type xiv"]):
            return "COLLAGEN_FACIT"
        else:
            return "STRUCTURAL"

    # --- Intermediate filaments ---
    if any(keyword in function_lower for keyword in ["keratin", "desmin", "vimentin", "intermediate filament"]):
        return "INTERMEDIATE_FILAMENT"

    # --- Fibrillin / microfibrils ---
    if "fibrillin" in function_lower or "microfibril" in function_lower:
        return "FIBRILLIN"

    # --- Elastin ---
    if "elastin" in function_lower:
        return "ELASTIN"

    # --- Actins / tubulins ---
    if any(keyword in function_lower for keyword in ["actin", "tubulin", "microtubule"]):
        return "CYTOSKELETON_POLYMER"

    # --- Lamins ---
    if "lamin" in function_lower:
        return "LAMIN"

    # --- RTKs / MAPK pathway ---
    if any(keyword in function_lower for keyword in ["receptor tyrosine kinase", "rtk", "ras", "raf", "mapk", "egfr", "fgfr"]):
        return "RTK_MAPK"

    # --- Negative regulators ---
    if any(keyword in function_lower for keyword in ["negative regulator", "repressor", "patched", "neurofibromin"]):
        return "NEGATIVE_REGULATOR"

    # Check for ion channel keywords (🔥 ADDED: calcium channel detection for RYR1!)
    ion_channel_keywords = [
        "calcium channel", "sodium channel", "potassium channel", "chloride channel",
        "ion channel", "calcium-activated", "voltage-gated", "ligand-gated"
    ]
    if any(keyword in function_lower for keyword in ion_channel_keywords):
        return "ION_CHANNEL"

    # Check GO terms and function for other families
    if any("channel" in term or "transport" in term for term in terms):
        return "ION_CHANNEL"
    elif any(term in terms for term in ["structural", "cytoskeleton", "basement membrane", "extracellular matrix"]):
        return "STRUCTURAL"
    elif any(term in terms for term in ["transcription", "methyltransferase", "dna-binding"]):
        return "TRANSCRIPTION_FACTOR"
    else:
        return "GENERAL"


# ======================================================
# STEP 2 + 4: PLAUSIBILITY FILTER
# ======================================================

def apply_pathogenicity_filter(
    raw_scores: Dict[str, float], gene_family: str
) -> Dict[str, Dict[str, Any]]:
    """
    Apply gene-family-specific plausibility rules to raw mechanism scores.

    Args:
        raw_scores: dict of {mechanism: score}
        gene_family: str, e.g. "ENZYME", "ION_CHANNEL"

    Returns:
        dict with mechanism → {raw_score, weighted_score, status, rationale}
    """
    rules = PATHOGENICITY_RULES.get(gene_family, PATHOGENICITY_RULES["GENERAL"])
    filtered: Dict[str, Dict[str, Any]] = {}

    for mech, score in raw_scores.items():
        weight = rules.get(mech, 0.0)
        weighted = score * weight

        if weight == 1.0:
            status = "kept"
        elif weight == 0.5 or weight == 0.25:
            status = "downweighted"
        else:
            status = "zeroed"

        rationale = f"{mech} {status} for {gene_family} (weight {weight})"

        filtered[mech] = {
            "raw_score": score,
            "weighted_score": weighted,
            "status": status,
            "rationale": rationale,
        }

    return filtered


# ======================================================
# STEP 6: PIPELINE FUNCTION
# ======================================================

def plausibility_pipeline(
    gene_symbol: str,
    raw_scores: Dict[str, float],
    uniprot_function: str = "",
    go_terms: List[str] = None,
) -> Dict[str, Any]:
    """
    Full plausibility pipeline: classify → filter → return explainable scores.

    Args:
        gene_symbol: str, e.g. "DLD"
        raw_scores: dict of {mechanism: score}
        uniprot_function: UniProt function description
        go_terms: list of GO terms (optional)

    Returns:
        dict with:
            - gene_family
            - filtered_scores (per mechanism with rationale)
            - final_scores (simplified weighted values)
    """
    if go_terms is None:
        go_terms = []

    # Step 1: raw_scores are provided
    # Step 2: classify gene family
    gene_family = classify_gene_family(gene_symbol, uniprot_function, go_terms)

    # Step 3-4: apply plausibility filter
    filtered = apply_pathogenicity_filter(raw_scores, gene_family)

    # Step 5-6: extract final weighted scores for downstream use
    final_scores = {mech: vals["weighted_score"] for mech, vals in filtered.items()}

    return {
        "gene_symbol": gene_symbol,
        "gene_family": gene_family,
        "filtered_scores": filtered,
        "final_scores": final_scores,
    }


# ======================================================
# Example usage
# ======================================================
if __name__ == "__main__":
    raw = {"DN": 0.45, "LOF": 0.62, "GOF": 0.84}
    result = plausibility_pipeline("DLD", raw, "pyruvate dehydrogenase enzyme", ["enzyme", "catalytic"])
    from pprint import pprint
    pprint(result)
