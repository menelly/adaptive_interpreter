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
    # --- BUGFIX 2026-04-26 (Ace): These labels are produced by classify_protein_weighted
    # via category_keywords.json but had NO entries in the rules table, silently falling
    # through to GENERAL (DN=0.5). That suppressed legitimate DN signals in POLG, ATP5F1A,
    # CDH1, NF2, BMPR2, ATP7B, CFTR and others. Per-family weights below reflect biology:
    "METABOLIC_ENZYME": {"LOF": 1.0, "DN": 1.0, "GOF": 0.0},     # mirror ENZYME (POLG, ATP5F1A, SDHB/C, DLD)
    "SIGNALING_REGULATOR": {"LOF": 1.0, "DN": 0.9, "GOF": 0.1},  # NF2/BMPR2/PCSK9 — mostly LOF, DN possible
    "SCAFFOLD_ADAPTOR": {"LOF": 1.0, "DN": 1.1, "GOF": 0.1},     # CDH1 — complex assembly disrupted by DN
    "TRANSPORTER": {"LOF": 1.0, "DN": 0.9, "GOF": 0.5},          # ATP7B/CFTR — channel/transport
    "AUTOSOMAL_RECESSIVE": {"LOF": 1.0, "DN": 0.7, "GOF": 0.3},  # AR-flagged: DN possible but LOF dominant
    # ----------------------------
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

    # 🛡️ CURATED OVERRIDES — for genes the keyword classifier consistently misclassifies.
    # Added 2026-04-27 (mechanism-refactor branch) after diagnosing that:
    #   ATP7A/ATP7B were getting METABOLIC_ENZYME (catalytic-keyword overlap) instead of TRANSPORTER
    #   MYH9 was getting CYTOSKELETON_POLYMER via 'actin' keyword instead of MOTOR_PROTEIN
    #   GAA was getting ELASTIN via priority-tie-breaking instead of METABOLIC_ENZYME
    #   CDH23/CDH1 were getting ION_CHANNEL via legacy GO-term fallback instead of SCAFFOLD_ADAPTOR
    # The keyword classifier and priority order need broader review, but for these
    # well-characterized genes, hardcoding the right answer is cleaner than tuning weights.
    CURATED_OVERRIDES = {
        # Transporters / channels (P-type ATPases, solute carriers)
        'ATP7A': 'TRANSPORTER',     # Menkes — copper transporter, NOT metabolic enzyme
        'ATP7B': 'TRANSPORTER',     # Wilson — copper transporter
        'SLC2A1': 'TRANSPORTER',    # GLUT1 — glucose transporter
        'CFTR':   'ION_CHANNEL',    # CFTR is genuinely a chloride channel
        # Motor proteins
        'MYH7':   'MOTOR_PROTEIN',  # cardiac myosin
        'MYH9':   'MOTOR_PROTEIN',  # non-muscle myosin IIA
        'MYO7A':  'MOTOR_PROTEIN',  # Usher 1B — unconventional myosin
        'MYO5B':  'MOTOR_PROTEIN',
        # Scaffold/adapter proteins (cadherins, catenins, junction proteins)
        'CDH1':   'SCAFFOLD_ADAPTOR',  # E-cadherin — cell adhesion
        'CDH23':  'SCAFFOLD_ADAPTOR',  # Usher 1D — cadherin (NOT ion channel)
        'CTNNB1': 'SCAFFOLD_ADAPTOR',  # β-catenin
        'CTNNA1': 'SCAFFOLD_ADAPTOR',  # α-catenin
        # Metabolic enzymes (lysosomal, mitochondrial, cytosolic)
        'GAA':    'METABOLIC_ENZYME',  # Pompe — alpha-glucosidase (NOT elastin)
        'HEXA':   'METABOLIC_ENZYME',  # Tay-Sachs — beta-hexosaminidase
        'GBA':    'METABOLIC_ENZYME',  # Gaucher — glucocerebrosidase
        'GALC':   'METABOLIC_ENZYME',  # Krabbe
        'PMM2':   'METABOLIC_ENZYME',  # CDG-Ia
        'MUT':    'METABOLIC_ENZYME',  # methylmalonic acidemia
        'POLG':   'METABOLIC_ENZYME',  # mitochondrial DNA polymerase
        'PAH':    'METABOLIC_ENZYME',  # PKU — phenylalanine hydroxylase
        'DLD':    'METABOLIC_ENZYME',  # dihydrolipoamide dehydrogenase
        'PYGL':   'METABOLIC_ENZYME',  # glycogen phosphorylase
        'ATP5F1A': 'METABOLIC_ENZYME', # mito ATP synthase alpha
        'PDE6B':  'METABOLIC_ENZYME',  # phosphodiesterase
        # Tumor suppressors (LOF haploinsufficient)
        'BRCA1':  'TUMOR_SUPPRESSOR',
        'BRCA2':  'TUMOR_SUPPRESSOR',
        'TP53':   'TUMOR_SUPPRESSOR',
        'RB1':    'TUMOR_SUPPRESSOR',
        'PTEN':   'TUMOR_SUPPRESSOR',
        'NF1':    'TUMOR_SUPPRESSOR',
        'NF2':    'TUMOR_SUPPRESSOR',
        'TSC1':   'TUMOR_SUPPRESSOR',
        'TSC2':   'TUMOR_SUPPRESSOR',
        'VHL':    'TUMOR_SUPPRESSOR',
        'APC':    'TUMOR_SUPPRESSOR',
        'CDKN2A': 'TUMOR_SUPPRESSOR',
        'STK11':  'TUMOR_SUPPRESSOR',
        # Mismatch repair (DNA_REPAIR family — currently has no JSON keywords; using TUMOR_SUPPRESSOR)
        # When DNA_REPAIR keywords get added to category_keywords.json, swap these.
        'MSH2':   'TUMOR_SUPPRESSOR',
        'MSH6':   'TUMOR_SUPPRESSOR',
        'MLH1':   'TUMOR_SUPPRESSOR',
        'MLH3':   'TUMOR_SUPPRESSOR',
        'PMS2':   'TUMOR_SUPPRESSOR',
        'MUTYH':  'TUMOR_SUPPRESSOR',
        # RAS pathway / RASopathies (MAPK signaling)
        'PTPN11': 'RTK_MAPK',
        'KRAS':   'RTK_MAPK',
        'NRAS':   'RTK_MAPK',
        'HRAS':   'RTK_MAPK',
        'BRAF':   'RTK_MAPK',
        'PIK3CA': 'ONCOGENE',
        # Ion channels — let weighted classifier handle most;
        # only override if it's clearly being misclassified
    }
    if gene_symbol in CURATED_OVERRIDES:
        chosen = CURATED_OVERRIDES[gene_symbol]
        logging.info(f"🛡️ Curated override: {gene_symbol} → {chosen}")
        return chosen
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

def _dn_evidence_score(gene_symbol: str, uniprot_function: str = "", go_terms: List[str] | None = None,
                       inheritance_patterns: List[str] | None = None) -> tuple[float, list[str]]:
    """Dominant-negative evidence score from UniProt data.
    Sources:
    - UniProt function text ("dominant negative" phrases)
    - GO terms (rare)
    - UniProt disease inheritance patterns (AD = DN evidence, AR-only = counter-evidence)
    Returns (score, hits). Positive = DN evidence, negative = AR-only counter-evidence.
    """
    if go_terms is None:
        go_terms = []
    if inheritance_patterns is None:
        inheritance_patterns = []
    score = 0.0
    hits: list[str] = []
    func = (uniprot_function or "").lower()

    # Core phrases in function text
    patterns = [
        "dominant negative",
        "dominant-negative",
        "exerts a dominant negative",
        "acts in a dominant negative",
    ]
    for p in patterns:
        if p in func:
            score += 1.0
            hits.append(f"uniprot_function:{p}")
            break

    # GO terms (rare direct hits)
    for term in go_terms:
        t = (term or "").lower()
        if "dominant-negative" in t or "dominant negative" in t:
            score += 0.5
            hits.append("go:dominant-negative")
            break

    # UniProt disease inheritance — the strongest signal
    if "AD" in inheritance_patterns:
        score += 1.5  # Strong: gene has known AD disease
        hits.append(f"uniprot_disease:AD_inheritance")
    if "AR" in inheritance_patterns and "AD" not in inheritance_patterns:
        score -= 0.5  # Counter-evidence: AR-only gene, DN less likely
        hits.append(f"uniprot_disease:AR_only_inheritance")

    return score, hits


def _adjust_dn_weight(base_weight: float, evidence_score: float) -> float:
    """Map evidence score → adjusted DN weight.
    - strong positive (>=1.0): at least 1.0 (AD disease known)
    - moderate positive (>=0.5): at least 0.9
    - negative (<0): halve the weight (AR-only, DN unlikely)
    - else: unchanged
    """
    if evidence_score >= 1.0:
        return max(base_weight, 1.0)
    if evidence_score >= 0.5:
        return max(base_weight, 0.9)
    if evidence_score < 0:
        return base_weight * 0.5  # AR-only counter-evidence
    return base_weight


def _self_association_score(uniprot_function: str = "", go_terms: List[str] | None = None) -> tuple[float, list[str]]:
    """Self-association / oligomerization evidence — the DN PRECONDITION.

    A dominant-negative variant works by co-assembling with its own wild-type copies and
    poisoning the shared assembly. That is IMPOSSIBLE for an obligate monomer — so DN must
    be gated on evidence that the protein self-associates (homo-oligomer, structural polymer,
    coiled-coil). We read this from GO terms (the cleanest available signal): homo-binding
    molecular functions and structural-unit COMPONENTS (collagen trimer, intermediate
    filament, microtubule). We deliberately exclude PROCESS terms ("regulation of...",
    "...organization") — those mean the protein ACTS ON a polymer (e.g. PIK3CA regulating
    actin), not that it IS one. Gene-agnostic; no hardcoding.

    Returns (score, hits): score > 0 == self-association evidence (DN mechanistically possible).
    """
    if go_terms is None:
        go_terms = []
    strong = ("identical protein binding", "homodimer", "homotetramer", "homotrimer",
              "homo-oligomer", "homooligomer", "homodimerization", "homotetramerization",
              "homo-oligomerization", "protein tetramerization", "tetramerization",
              "oligomerization", "intermediate filament", "microtubule", "collagen type",
              "collagen trimer", "myosin filament", "homophilic")
    disq = ("regulation", "organization", "depolymeriz", "based movement", "signaling",
            "biosynth", "catabolic", "transport", "binding to")
    score = 0.0
    hits: list[str] = []
    for t in go_terms:
        tl = (t or "").lower()
        if any(k in tl for k in strong):
            if tl == "identical protein binding" or not any(d in tl for d in disq):
                score += 0.5
                hits.append(f"go:{t}")
    fl = (uniprot_function or "").lower()
    if any(k in fl for k in ("homodimer", "homotrimer", "homotetramer", "homo-oligomer",
                             "self-assembl", "forms a coiled coil", "triple helix")):
        score += 0.5
        hits.append("function:self_association_language")
    return min(score, 1.5), hits


def _adjust_gof_weight(base_weight: float, evidence_score: float) -> float:
    """Map GOF licensing-evidence score → adjusted GOF weight (mirror of _adjust_dn_weight).

    GOF is the rare mechanism, so the coarse family weight stands UNLESS the protein's own
    molecular evidence licenses it. This is what lets a protein whose disease mechanism is
    atypical for its structural family (e.g. CTNNB1 — a SCAFFOLD_ADAPTOR that GOFs via
    degron-loss stabilization) express a real GOF signal that the family bucket would mute.
    Scale matches _gof_evidence_score: gain-language=+1.5, capable-class=+0.8, AR-only=-0.5.
    - strong gain-language (>=1.5): at least 1.0
    - moderate capable-class (>=0.5): at least 0.9
    - negative (veto family / AR-only, no gain): halve the weight
    - else: unchanged
    """
    if evidence_score >= 1.5:
        return max(base_weight, 1.0)
    if evidence_score >= 0.5:
        return max(base_weight, 0.9)
    if evidence_score < 0:
        return base_weight * 0.5
    return base_weight


# Families where GOF is biologically nonsensical (loss/structural). A bare capable-class
# keyword must NOT license GOF for these — only explicit gain-language can (rare but real).
_GOF_VETO_FAMILIES = {
    "METABOLIC_ENZYME", "ENZYME", "DNA_REPAIR", "RIBOSOMAL_PROTEIN", "TUMOR_SUPPRESSOR",
    "STRUCTURAL", "INTERMEDIATE_FILAMENT", "FIBRILLIN", "ELASTIN", "MUSCULAR_DYSTROPHY",
}


def _gof_evidence_score(gene_symbol: str, uniprot_function: str = "", go_terms: List[str] | None = None,
                        inheritance_patterns: List[str] | None = None,
                        disease_text: str = "", gene_family: str = "") -> tuple[float, list[str]]:
    """GOF licensing-evidence score — the GOF analog of _dn_evidence_score.

    Design (Ren + Nova + Ace, 2026-05-20): GOF is the RARE mechanism. A random missense is
    far more likely to break a protein (LOF) than to give it new/enhanced function. So GOF
    must be *licensed* by a positive signal, not assumed. This returns:
      score > 0  → GOF biologically plausible (licensed)
      score <= 0 → GOF implausible (no signal, or recessive counter-evidence)
    Sources searched: UniProt function text, GO terms, AND disease descriptions (the gain-
    language for genes like MEFV — "autoinflammatory" — lives in the disease text, not the
    function text). Returns (score, hits) for transparency in the output explanation.

    NOTE: pure predicate. Does not change any score on its own — the interpretation gate
    (built + calibration-validated separately) consumes this.
    """
    if go_terms is None:
        go_terms = []
    if inheritance_patterns is None:
        inheritance_patterns = []
    score = 0.0
    hits: list[str] = []
    blob = " ".join([(uniprot_function or ""), " ".join(go_terms), (disease_text or "")]).lower()
    veto_family = (gene_family or "").upper() in _GOF_VETO_FAMILIES or (gene_family or "").upper().startswith("COLLAGEN")

    # Strong: explicit gain-language anywhere in function/disease text. Trusted even for
    # veto families (e.g. a metabolic gene with a documented gain-of-function is rare but real).
    gain_phrases = [
        "gain-of-function", "gain of function", "gain-of-function mutation",
        "constitutive", "constitutively active", "constitutive activation",
        "ligand-independent", "ligand independent", "autoinflammatory",
        "increased activity", "enhanced activity", "hyperactiv",
    ]
    has_gain_language = False
    for p in gain_phrases:
        if p in blob:
            score += 1.5
            hits.append(f"gain_language:{p}")
            has_gain_language = True
            break  # one strong hit is enough; don't stack

    # Moderate: molecular class that CAN gain function (receptor / channel / GTPase / kinase / TF).
    # Suppressed for veto families — a metabolic enzyme with a "channel" domain (e.g. ATP
    # synthase) must not be licensed by the keyword alone.
    capable_classes = [
        "receptor", "channel", "gtpase", "guanine nucleotide", "kinase",
        "transcription factor", "signal transduction", "g protein",
    ]
    if not veto_family:
        for c in capable_classes:
            if c in blob:
                score += 0.8
                hits.append(f"gof_capable_class:{c}")
                break

    # Counter-evidence: AR-only inheritance — recessive disease is almost always loss.
    # (Gain-language above can still license it — the escape hatch for MEFV-type genes.)
    if "AR" in inheritance_patterns and "AD" not in inheritance_patterns:
        score -= 0.5
        hits.append("inheritance:AR_only_gof_counter")

    # Veto family with no explicit gain-language → force implausible.
    if veto_family and not has_gain_language:
        hits.append(f"family_veto:{(gene_family or '').upper()}")
        score = min(score, -0.5)

    return score, hits


def _adjust_gof_weight(base_weight: float, evidence_score: float) -> float:
    """Map GOF evidence score → adjusted GOF weight (mirror of _adjust_dn_weight).
    - strong positive (>=1.0, licensed): at least 1.0 — OVERRIDES a low family weight,
      so a mislabeled gene (MEFV → CYTOSKELETON_POLYMER, base ~0.1) can still gain.
    - moderate positive (>=0.5): at least 0.8
    - negative (implausible / veto): suppress to <=0.1
    - else (no signal): unchanged (defer to family weight)
    """
    if evidence_score >= 1.0:
        return max(base_weight, 1.0)
    if evidence_score >= 0.5:
        return max(base_weight, 0.8)
    if evidence_score < 0:
        return min(base_weight, 0.1)
    return base_weight


def _annotation_desert(go_terms: List[str] | None, inheritance_patterns: List[str] | None,
                       disease_text: str = "") -> bool:
    """True when we have NO biological context to reason from — no GO terms, no disease,
    no inheritance. The interpretation gate is then "tossing paint at a wall" and confidence
    must be floored (~0.25) with an honest explanation. (Ren, 2026-05-20.)"""
    return (not go_terms) and (not inheritance_patterns) and (not (disease_text or "").strip())


def _is_negative_regulator(uniprot_function: str = "", go_terms: List[str] | None = None) -> bool:
    """True when the protein's job is to inhibit/repress/negatively-regulate. This is the
    'brakes' case: a dominant-negative hit on a brake RELEASES it = functional gain, so DN
    synergy is licensed in the GOF branch only for these genes (not for plain receptors)."""
    if go_terms is None:
        go_terms = []
    blob = " ".join([(uniprot_function or ""), " ".join(go_terms)]).lower()
    neg_phrases = [
        "negative regulation", "negative regulator", "negatively regulat",
        "inhibitor of", "inhibits", "repressor", "represses", "suppressor of",
        "down-regulat", "downregulat", "antagonist",
    ]
    return any(p in blob for p in neg_phrases)


def apply_pathogenicity_filter(
    raw_scores: Dict[str, float],
    gene_family: str,
    gene_symbol: str = "",
    uniprot_function: str = "",
    go_terms: List[str] | None = None,
    inheritance_patterns: List[str] | None = None,
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

    # Compute DN + GOF licensing-evidence once per gene (includes inheritance patterns).
    # These let a mechanism that is atypical for the structural family still express when the
    # protein's OWN molecular evidence supports it — gene-agnostic, no hardcoding.
    dn_ev_score, dn_ev_hits = _dn_evidence_score(gene_symbol, uniprot_function, go_terms, inheritance_patterns)
    gof_ev_score, gof_ev_hits = _gof_evidence_score(gene_symbol, uniprot_function, go_terms,
                                                    inheritance_patterns, gene_family=gene_family)
    # DN precondition: can this protein even self-assemble? (DN poisons a shared assembly;
    # an obligate monomer cannot.) Soft-gates DN below, overridable by explicit DN text.
    self_assoc_score, self_assoc_hits = _self_association_score(uniprot_function, go_terms)

    for mech, score in raw_scores.items():
        # 🎯 NORMALIZE: family weights are capped at 1.0. Some rules historically used
        # >1.0 weights as multiplicative BOOSTS (e.g. COLLAGEN_FIBRILLAR DN=1.25), which
        # pushed scores past [0,1] and made families incomparable (collagen DN floated 25%
        # above everyone). The boost is applied equally to P and B within a family, so it
        # added zero discrimination — verified 2026-05-20: capping leaves per-gene AUC
        # byte-identical (0.764) while pooled AUC rose 0.687->0.721 and DN returned to [0,1].
        # A family that "favors DN" should show it via a high RAW score, not a post-hoc boost.
        base_weight = min(rules.get(mech, 0.0), 1.0)
        weight = base_weight
        evidence_applied = False
        if mech == "DN":
            new_weight = _adjust_dn_weight(base_weight, dn_ev_score)
            # Oligomerization precondition (SOFT, overridable): DN needs self-assembly to
            # poison. If there's no self-association evidence AND no explicit "dominant
            # negative" text licensing it, halve the DN weight — don't zero it, because
            # absent annotation != confirmed monomer (Ren's soft-branch rule).
            explicit_dn = any(h.startswith("uniprot_function:") for h in dn_ev_hits)
            if self_assoc_score <= 0 and not explicit_dn:
                new_weight = new_weight * 0.5
            if new_weight != base_weight:
                weight = new_weight
                evidence_applied = True
        elif mech == "GOF":
            new_weight = _adjust_gof_weight(base_weight, gof_ev_score)
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
                "self_association_score": self_assoc_score,
                "self_association_hits": self_assoc_hits,
            })
        elif mech == "GOF":
            entry.update({
                "gof_weight_base": base_weight,
                "gof_weight_final": weight,
                "gof_evidence_score": gof_ev_score,
                "gof_evidence_hits": gof_ev_hits,
                "gof_evidence_applied": evidence_applied,
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
    inheritance_patterns: List[str] = None,
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

    # Step 3-4: apply plausibility filter (now includes atypical_flag + inheritance)
    filtered = apply_pathogenicity_filter(
        raw_scores,
        gene_family,
        gene_symbol=gene_symbol,
        uniprot_function=uniprot_function,
        go_terms=go_terms,
        inheritance_patterns=inheritance_patterns,
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
