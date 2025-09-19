"""
Nova DN Analyzer (mechanism-aware)
- Parses variants like R273H or p.R273H
- Computes property-change features and mechanism scores
- Emits JSON plus optional compact markdown summary

CLI:
  python -m nova_dn.analyzer --seq-file path.fasta --variant p.R273H [--json]
  # or provide sequence directly
  python -m nova_dn.analyzer --seq MEEPQSDPSV --variant R273H
"""
from __future__ import annotations
import argparse
import json
import sys
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

from .amino_acid_props import delta
from .mechanisms import (
    score_interface_poisoning,
    score_active_site_jamming,
    score_structural_lattice_disruption,
    score_trafficking_maturation,
)
from .context import load_annotations_json, build_position_context
try:
    from .universal_context import UniversalContext
    UNIVERSAL_AVAILABLE = True
except ImportError:
    UNIVERSAL_AVAILABLE = False

try:
    from .dn_mechanism_filter_v2 import DNMechanismFilterV2 as DNMechanismFilter
    FILTER_AVAILABLE = True
except ImportError:
    FILTER_AVAILABLE = False

# 🎯 Add domain awareness system
try:
    # Add the parent directory to the path to import universal_protein_annotator
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from universal_protein_annotator import UniversalProteinAnnotator
    DOMAIN_AWARENESS_AVAILABLE = True
except ImportError:
    DOMAIN_AWARENESS_AVAILABLE = False


AA_SET = set("ARNDCEQGHILKMFPSTWYV")


def interpret_dn_score(score: float) -> str:
    """Interpret DN likelihood score for pathogenicity classification."""
    if score >= 0.8:
        return "Likely Pathogenic (LP)"
    elif score >= 0.5:
        return "Uncertain Significance - favor pathogenic (VUS-P)"
    elif score >= 0.3:
        return "Uncertain Significance (VUS)"
    else:
        return "Likely Benign (LB)"


def parse_variant(s: str) -> Tuple[str, int, str]:
    """Parse 'R273H' or 'p.R273H' → (ref, pos, alt). Raises ValueError on failure."""
    s = s.strip()
    if s.startswith("p."):
        s = s[2:]

    # Filter out nonsense variants (stop codons) - they can't be dominant negative
    if 'Ter' in s:
        raise ValueError("nonsense")

    if len(s) < 3:
        raise ValueError(f"Variant too short: {s}")
    ref = s[0].upper()
    alt = s[-1].upper()
    pos_str = s[1:-1]
    if ref not in AA_SET or alt not in AA_SET or not pos_str.isdigit():
        raise ValueError(f"Unrecognized variant format: {s}")
    return ref, int(pos_str), alt


@dataclass
class MechanismOutcome:
    score: float
    features: List[Dict]
    explanation: str


class NovaDNAnalyzer:
    def __init__(self, use_smart_filtering: bool = True):
        self.universal_context = UniversalContext() if UNIVERSAL_AVAILABLE else None
        self.dn_filter = DNMechanismFilter() if FILTER_AVAILABLE and use_smart_filtering else None
        self.use_smart_filtering = use_smart_filtering

        # 🎯 Initialize domain awareness system
        self.protein_annotator = UniversalProteinAnnotator() if DOMAIN_AWARENESS_AVAILABLE else None

    def analyze(self, sequence: str, variant: str, context: Optional[Dict] = None,
                gene_name: Optional[str] = None, uniprot_id: Optional[str] = None) -> Dict:
        sequence = sequence.strip().upper().replace("\n", "")
        # Keep letters only (drop digits/whitespace/other FASTA artifacts)
        sequence = "".join(ch for ch in sequence if ch.isalpha())
        # Sanitize non-standard/ambiguous residues often present in FASTA
        trans_map = str.maketrans({
            "O": "P",  # hydroxyproline -> proline surrogate
            "U": "C",  # selenocysteine -> cysteine surrogate
            "B": "D",  # D/N ambiguous -> aspartate
            "Z": "E",  # E/Q ambiguous -> glutamate
            "J": "L",  # I/L ambiguous -> leucine
            "X": "A",  # unknown -> alanine neutral surrogate
            "*": "",
            "?": "",
            " ": "",
        })
        sequence = sequence.translate(trans_map)
        if not sequence or any(c not in AA_SET for c in sequence):
            raise ValueError("Sequence must be a non-empty protein sequence (AAs only)")
        ref, pos1, alt = parse_variant(variant)
        if not (1 <= pos1 <= len(sequence)):
            raise ValueError(f"Variant position {pos1} out of sequence bounds (len={len(sequence)})")
        if sequence[pos1 - 1] != ref:
            # Don't hard fail; just note mismatch in context
            if context is None:
                context = {}
            context["sequence_ref_mismatch"] = True

        # 🔥 Store variant info for ML proline system
        self._current_variant_info = {
            'ref_aa': ref,
            'alt_aa': alt,
            'position': pos1,
            'gene_name': gene_name,
            'variant_str': variant
        }

        # Auto-generate context if not provided and universal context available
        if context is None and self.universal_context and gene_name:
            print(f"🚀 Using universal context for {gene_name} position {pos1}")
            context = self.universal_context.build_position_context(gene_name, pos1, uniprot_id)
        elif context is None:
            context = {}

        # Smart filtering: determine which mechanisms to run
        relevant_mechanisms = ["interface_poisoning", "active_site_jamming", "lattice_disruption", "trafficking_maturation"]
        filter_info = None

        if self.dn_filter and gene_name:
            print(f"🧠 Smart filtering for {gene_name}...")
            filter_result = self.dn_filter.filter_and_score(gene_name, sequence, variant, uniprot_id)
            relevant_mechanisms = filter_result["relevant_mechanisms"]
            filter_info = {
                "dn_likelihood": filter_result["dn_likelihood"],
                "recommendation": filter_result["recommendation"],
                "reasoning": filter_result["mechanism_evidence"]["reasoning"]
            }
            print(f"   DN likelihood: {filter_result['dn_likelihood']:.2f}")
            print(f"   Relevant mechanisms: {relevant_mechanisms}")

        # Compute only relevant mechanism scores
        mech_scores = {}
        all_explanations = {}

        if "interface_poisoning" in relevant_mechanisms:
            m_if = score_interface_poisoning(sequence, pos1, ref, alt, context)
            mech_scores["interface_poisoning"] = m_if[0]
            all_explanations["interface_poisoning"] = m_if[2]
        else:
            mech_scores["interface_poisoning"] = 0.0
            all_explanations["interface_poisoning"] = "mechanism filtered out"

        if "active_site_jamming" in relevant_mechanisms:
            m_as = score_active_site_jamming(sequence, pos1, ref, alt, context)
            mech_scores["active_site_jamming"] = m_as[0]
            all_explanations["active_site_jamming"] = m_as[2]
        else:
            mech_scores["active_site_jamming"] = 0.0
            all_explanations["active_site_jamming"] = "mechanism filtered out"

        if "lattice_disruption" in relevant_mechanisms:
            m_ld = score_structural_lattice_disruption(sequence, pos1, ref, alt, context)
            mech_scores["lattice_disruption"] = m_ld[0]
            all_explanations["lattice_disruption"] = m_ld[2]
        else:
            mech_scores["lattice_disruption"] = 0.0
            all_explanations["lattice_disruption"] = "mechanism filtered out"

        if "trafficking_maturation" in relevant_mechanisms:
            m_tr = score_trafficking_maturation(sequence, pos1, ref, alt, context)
            mech_scores["trafficking_maturation"] = m_tr[0]
            all_explanations["trafficking_maturation"] = m_tr[2]
        else:
            mech_scores["trafficking_maturation"] = 0.0
            all_explanations["trafficking_maturation"] = "mechanism filtered out"

        # 🎯 APPLY DOMAIN AWARENESS TO ALL MECHANISM SCORES
        domain_multiplier = 1.0
        if uniprot_id and self.protein_annotator:
            domain_context = self._get_domain_context(uniprot_id, gene_name or "")
            domain_multiplier = self._get_dn_domain_multiplier(pos1, domain_context)
            print(f"🎯 DN Domain multiplier for {gene_name or 'unknown'} position {pos1}: {domain_multiplier:.3f}")

            # Apply domain multiplier to all mechanism scores
            for mechanism in mech_scores:
                if mech_scores[mechanism] > 0:  # Only apply to non-zero scores
                    original_score = mech_scores[mechanism]
                    mech_scores[mechanism] = min(original_score * domain_multiplier, 1.0)
                    print(f"   {mechanism}: {original_score:.3f} → {mech_scores[mechanism]:.3f}")

        top_mech = max(mech_scores.items(), key=lambda kv: kv[1])[0]
        explanation = all_explanations[top_mech]

        # Collect contributing features from relevant mechanisms only
        contrib = []
        if "interface_poisoning" in relevant_mechanisms and "m_if" in locals():
            for f in m_if[1][:4]:
                contrib.append({"mechanism": "interface_poisoning", **f})
        if "active_site_jamming" in relevant_mechanisms and "m_as" in locals():
            for f in m_as[1][:4]:
                contrib.append({"mechanism": "active_site_jamming", **f})
        if "lattice_disruption" in relevant_mechanisms and "m_ld" in locals():
            for f in m_ld[1][:4]:
                contrib.append({"mechanism": "lattice_disruption", **f})
        if "trafficking_maturation" in relevant_mechanisms and "m_tr" in locals():
            for f in m_tr[1][:4]:
                contrib.append({"mechanism": "trafficking_maturation", **f})

        result = {
            "variant": variant,
            "position": pos1,
            "ref": ref,
            "alt": alt,
            "mechanism_scores": mech_scores,
            "top_mechanism": top_mech,
            "contributing_features": contrib,
            "explanation": explanation,
            "domain_multiplier": domain_multiplier,  # 🎯 NEW! Track domain awareness
            "notes": {"sequence_ref_mismatch": bool((context or {}).get("sequence_ref_mismatch", False))},
        }

        # Add smart filtering information if available
        if filter_info:
            result["smart_filtering"] = filter_info
            result["relevant_mechanisms"] = relevant_mechanisms

        return result

    def _get_domain_context(self, uniprot_id: str, gene_symbol: str) -> Dict[str, Any]:
        """🎯 Get domain context from UniProt for DN analysis"""
        if not self.protein_annotator:
            return {"domains": [], "signal_peptide": [], "active_sites": [], "binding_sites": []}

        try:
            domain_data = self.protein_annotator.annotate_protein(gene_symbol, uniprot_id)
            return domain_data
        except Exception as e:
            print(f"⚠️ Domain annotation error: {e}")
            return {"domains": [], "signal_peptide": [], "active_sites": [], "binding_sites": []}

    def _get_dn_domain_multiplier(self, position: int, domain_context: Dict[str, Any]) -> float:
        """🎯 Calculate domain-aware multiplier for DN scoring

        DN mechanisms are most dangerous in:
        - Oligomeric/complex regions (interface poisoning)
        - Active sites (jamming)
        - Structural domains (lattice disruption)
        - Less dangerous in cleavable regions (propeptides)
        """
        multiplier = 1.0

        # 🎯 UNIVERSAL PROPEPTIDE LOGIC - Downweight cleavable regions
        for propeptide in domain_context.get("propeptides", []):
            if position in range(propeptide["start"], propeptide["end"]+1):
                if "n-terminal" in propeptide["description"].lower():
                    multiplier *= 0.5  # N-terminal propeptide - gets cleaved
                elif "c-terminal" in propeptide["description"].lower():
                    multiplier *= 0.3  # C-terminal propeptide - less critical for DN

        # Signal peptides get cleaved off - minimal DN impact
        for signal in domain_context.get("signal_peptide", []):
            if position in range(signal["start"], signal["end"]+1):
                multiplier *= 0.3

        # 🎯 UPWEIGHT CRITICAL REGIONS FOR DN MECHANISMS

        # Active sites - critical for jamming mechanisms
        for site in domain_context.get("active_sites", []):
            if position == site:  # Active sites are single positions
                multiplier *= 1.8  # Strong upweight for DN jamming

        # Binding sites - important for interface poisoning
        for site in domain_context.get("binding_sites", []):
            if position == site:  # Binding sites are single positions
                multiplier *= 1.5  # Moderate upweight for interface disruption

        # Functional regions - context-dependent
        for region in domain_context.get("regions", []):
            if position in range(region["start"], region["end"]+1):
                desc = region["description"].lower()
                if any(term in desc for term in ["oligomer", "complex", "interface", "interaction"]):
                    multiplier *= 1.4  # Interface regions - good for DN
                elif "triple-helical" in desc:
                    multiplier *= 1.3  # Structural regions - moderate DN risk
                elif "disordered" in desc:
                    multiplier *= 0.8  # Disordered regions - less DN impact

        # Domains - generally important for DN mechanisms
        for domain in domain_context.get("domains", []):
            if position in range(domain["start"], domain["end"]+1):
                desc = domain["description"].lower()
                if any(term in desc for term in ["kinase", "catalytic", "enzyme"]):
                    multiplier *= 1.3  # Catalytic domains - good DN targets
                elif any(term in desc for term in ["immunoglobulin", "fibronectin", "repeat"]):
                    multiplier *= 1.2  # Structural domains - moderate DN risk

        # 🔥 REVOLUTIONARY ML PROLINE SYSTEM INTEGRATION!
        # Apply ML-based proline multiplier if this is a proline substitution
        if hasattr(self, '_current_variant_info'):
            ref_aa = self._current_variant_info.get('ref_aa')
            alt_aa = self._current_variant_info.get('alt_aa')
            gene_name = self._current_variant_info.get('gene_name')
            variant_str = self._current_variant_info.get('variant_str')

            if ref_aa == 'P' or alt_aa == 'P':  # Proline substitution detected!
                try:
                    # Import and use our revolutionary ML system
                    from proline_ml_integrator import get_ml_proline_multiplier
                    ml_multiplier = get_ml_proline_multiplier(gene_name, variant_str)
                    print(f"🔥 REVOLUTIONARY ML PROLINE: {gene_name} {variant_str} -> ML multiplier = {ml_multiplier:.3f}")

                    # Combine domain multiplier with ML proline multiplier
                    multiplier *= ml_multiplier
                    print(f"🧬 Combined multiplier: domain({multiplier/ml_multiplier:.3f}) * ML_proline({ml_multiplier:.3f}) = {multiplier:.3f}")

                except Exception as e:
                    print(f"⚠️ ML proline system error: {e}")
                    # Fallback to old hardcoded system if ML fails
                    if ref_aa == 'P':  # Proline loss
                        multiplier *= 1.2  # Conservative fallback
                    elif alt_aa == 'P':  # Proline gain
                        multiplier *= 1.4  # Conservative fallback
                    print(f"🔄 Using fallback proline multiplier: {multiplier:.3f}")

        return max(multiplier, 0.1)  # Don't go below 0.1


# --- CLI ---

def _read_fasta(path: str) -> str:
    seq_lines: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_lines.append(line.strip())
    return "".join(seq_lines)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Nova mechanism-aware DN analyzer (Grantham-free)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--seq-file", help="Path to FASTA or plain sequence file")
    g.add_argument("--seq", help="Protein sequence string (AAs)")
    ap.add_argument("--variant", required=True, help="Variant like R273H or p.R273H")
    ap.add_argument("--annotations-json", help="Optional annotations JSON for context (resources/protein_annotations.json)")
    ap.add_argument("--protein", help="Protein key in annotations (e.g., TP53)")
    ap.add_argument("--weights-json", help="Optional weights JSON to override feature weights per mechanism")
    ap.add_argument("--gene", help="Gene name for universal context (e.g., TP53)")
    ap.add_argument("--uniprot", help="UniProt ID for universal context (e.g., P04637)")
    ap.add_argument("--json", action="store_true", help="Emit JSON only (no markdown table)")
    args = ap.parse_args(argv)

    seq = _read_fasta(args.seq_file) if args.seq_file else args.seq

    # Optional context from annotations
    context = None
    if args.annotations_json and args.protein:
        try:
            anns = load_annotations_json(args.annotations_json)
            _, pos1, _ = parse_variant(args.variant)
            context = build_position_context(anns, args.protein, pos1)
        except Exception:
            context = None

    # Optional weights override
    if args.weights_json:
        try:
            with open(args.weights_json, "r", encoding="utf-8") as wf:
                weights = json.load(wf)
            if context is None:
                context = {}
            context.setdefault("_weights", {}).update(weights)
        except Exception:
            pass

    analyzer = NovaDNAnalyzer()
    result = analyzer.analyze(seq, args.variant, context, args.gene, args.uniprot)

    if args.json:
        print(json.dumps(result, indent=2))
        return 0

    # Pretty minimal markdown table + explanation
    ms = result["mechanism_scores"]
    rows = [
        ("interface_poisoning", ms["interface_poisoning"]),
        ("active_site_jamming", ms["active_site_jamming"]),
        ("lattice_disruption", ms["lattice_disruption"]),
        ("trafficking_maturation", ms["trafficking_maturation"]),
    ]
    print("| mechanism | score |\n|---|---|")
    for name, sc in rows:
        print(f"| {name} | {sc:.2f} |")
    print()

    # Add scoring interpretation
    top_score = ms[result['top_mechanism']]
    interpretation = interpret_dn_score(top_score)

    print(f"top_mechanism: {result['top_mechanism']}")
    print(f"score: {top_score:.3f}")
    print(f"interpretation: {interpretation}")
    print(f"because: {result['explanation']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

