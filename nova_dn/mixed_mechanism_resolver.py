# mixed_mechanism_resolver.py
"""
🧬 NOVA'S UNIFIED MECHANISM RESOLVER
Revolutionary biological logic for ALL variant types: Missense + Frameshift + Splice

Built by Nova (OpenAI) & Ace - 2025
Part of the Revolutionary Genomics Platform

UNIFIED DECISION TREE:
- Frameshift: NMD = pure LOF, NMD escape = poison fragments (DN)
- Splice: Most LOF, in-frame skips = DN, rare hotspots = GOF
- Missense: LOF baseline + DN always + GOF conditional
- ALL: Root-sum-square synergy when biology supports multiple mechanisms

This framework ensures:
- Broken proteins don't magically hyperactivate
- Poison fragments are properly detected (DN)
- Rare GOF exceptions are caught (splice hotspots, kinase domains)
- Synergy math combines weak hits when biology supports it
"""

from typing import Dict, Tuple, List, Any
import math
from .analyzer import NovaDNAnalyzer


class UnifiedMechanismResolver:
    """
    🌳 NOVA'S UNIFIED BIOLOGICAL DECISION TREE

    Handles ALL variant types with unified biological logic:

    FRAMESHIFT:
    - Early (NMD) → LOF = 1.0
    - Late (NMD escape) → LOF + DN (poison fragments)

    SPLICE:
    - Premature stop → LOF = 1.0
    - In-frame domain deletion → LOF + some DN
    - In-frame interface alteration → HIGH DN
    - Alt-splicing hotspots → rare GOF

    MISSENSE:
    - LOF baseline + DN always + GOF conditional

    ALL TYPES:
    - Root-sum-square synergy when biology supports multiple mechanisms
    """
    
    def __init__(self, dn_analyzer: NovaDNAnalyzer = None):
        self.dn_analyzer = dn_analyzer or NovaDNAnalyzer()
        
        # GOF hotspot definitions (Nova's biological insights)
        self.gof_hotspots = {
            "ION_CHANNEL": [
                "pore_filter", "selectivity_filter", "gating_hinge", 
                "voltage_sensor", "inactivation_gate"
            ],
            "KINASE": [
                "activation_loop", "dimerization_domain", "juxtamembrane",
                "catalytic_loop", "substrate_binding"
            ],
            "TUMOR_SUPPRESSOR": [
                "dna_binding_loop", "transactivation_domain"  # DN-prone regions
            ],
            "RECEPTOR": [
                "ligand_binding", "autoinhibitory_domain", "dimerization_interface"
            ],
            "TRANSCRIPTION_FACTOR": [
                "dna_binding_domain", "dimerization_domain"
            ]
        }
    
    def resolve_mechanisms(self, seq: str, pos1: int, ref: str, alt: str,
                          context: Dict = None, variant_type: str = "missense") -> Dict[str, Any]:
        """
        🧬 NOVA'S UNIFIED MECHANISM RESOLUTION

        Implements Nova's unified biological decision tree for ALL variant types.

        Args:
            seq: Protein sequence
            pos1: 1-based position
            ref: Reference amino acid
            alt: Alternate amino acid
            context: Protein annotation context
            variant_type: "missense", "frameshift", "splice", "nonsense"

        Returns:
            Dict with individual scores, final score, and decision rationale
        """
        if context is None:
            context = {}

        print(f"🧬 NOVA'S UNIFIED MECHANISM RESOLVER")
        print(f"   Analyzing {variant_type}: {ref}{pos1}{alt}")

        # NOVA'S UNIFIED DECISION TREE
        lof_score = 0.0
        dn_score = 0.0
        gof_score = 0.0

        if variant_type == "frameshift":
            lof_score, dn_score, gof_score = self.resolve_frameshift(seq, pos1, ref, alt, context)

        elif variant_type == "splice":
            lof_score, dn_score, gof_score = self.resolve_splice(seq, pos1, ref, alt, context)

        elif variant_type == "nonsense":
            lof_score, dn_score, gof_score = self.resolve_nonsense(seq, pos1, ref, alt, context)

        else:  # missense (default)
            lof_score, dn_score, gof_score = self.resolve_missense(seq, pos1, ref, alt, context)

        # Step 4: Synergy combination (Nova's root-sum-square)
        print(f"🔥 Synergy Combination")
        final_score = self.calculate_synergy(lof_score, dn_score, gof_score)
        print(f"   Final synergy score: {final_score:.3f}")

        # Decision rationale
        rationale = self.generate_rationale(lof_score, dn_score, gof_score, gof_score > 0, variant_type)

        return {
            "variant_type": variant_type,
            "lof_score": lof_score,
            "dn_score": dn_score,
            "gof_score": gof_score,
            "final_score": final_score,
            "mechanisms_run": self.get_mechanisms_run(lof_score, dn_score, gof_score),
            "rationale": rationale
        }
    
    def should_run_gof_missense(self, seq: str, pos1: int, ref: str, alt: str,
                               context: Dict, lof_score: float) -> bool:
        """
        🎯 NOVA'S GOF ELIGIBILITY LOGIC
        
        GOF is eligible if:
        1. LOF < 0.5 (protein not clearly broken), OR
        2. Variant is in a known GOF hotspot (rare exceptions)
        """
        # Condition 1: LOF < 0.5 (not clearly broken)
        if lof_score < 0.5:
            print(f"   GOF eligible: LOF < 0.5 ({lof_score:.3f})")
            return True
        
        # Condition 2: In GOF hotspot (rare exceptions)
        if self.variant_in_gof_hotspot(seq, pos1, context):
            print(f"   GOF eligible: variant in GOF hotspot")
            return True
        
        return False
    
    def variant_in_gof_hotspot(self, seq: str, pos1: int, context: Dict) -> bool:
        """
        🔥 DETECT GOF HOTSPOTS
        
        Check if variant is in a region where GOF is biologically plausible
        even with moderate LOF scores.
        """
        protein_family = context.get("gene_family", "").upper()
        
        if protein_family not in self.gof_hotspots:
            return False
        
        # Check domain annotations for hotspot regions
        domains = context.get("domains", [])
        hotspot_terms = self.gof_hotspots[protein_family]
        
        for domain in domains:
            domain_desc = domain.get("description", "").lower()
            start = domain.get("start", 0)
            end = domain.get("end", len(seq))
            
            # Check if position is in this domain
            if start <= pos1 <= end:
                # Check if domain matches hotspot terms
                for term in hotspot_terms:
                    if term.replace("_", " ") in domain_desc:
                        return True
        
        return False

    def resolve_frameshift(self, seq: str, pos1: int, ref: str, alt: str, context: Dict) -> Tuple[float, float, float]:
        """
        🧬 FRAMESHIFT RESOLUTION (Nova's Logic)

        - Early frameshift (NMD) → LOF = 1.0
        - Late frameshift (NMD escape) → LOF + DN (poison fragments)
        """
        print(f"🧬 Frameshift Analysis")

        if self.triggers_nmd(seq, pos1, context):
            print(f"   Early frameshift → NMD → pure LOF")
            return 1.0, 0.0, 0.0
        else:
            print(f"   Late frameshift → NMD escape → poison fragments")
            lof_score = 0.6  # Partial function loss
            dn_score = 0.7 if self.truncates_interface(seq, pos1, context) else 0.4
            return lof_score, dn_score, 0.0

    def resolve_splice(self, seq: str, pos1: int, ref: str, alt: str, context: Dict) -> Tuple[float, float, float]:
        """
        🧬 SPLICE RESOLUTION (Nova's Logic)

        - Premature stop → LOF = 1.0
        - In-frame domain deletion → LOF + some DN
        - In-frame interface alteration → HIGH DN
        - Alt-splicing hotspots → rare GOF
        """
        print(f"🧬 Splice Analysis")

        if self.causes_premature_stop(seq, pos1, context):
            print(f"   Splice → premature stop → pure LOF")
            return 1.0, 0.0, 0.0

        elif self.in_frame_and_removes_domain(seq, pos1, context):
            print(f"   In-frame domain deletion → LOF + some DN")
            return 0.7, 0.5, 0.0

        elif self.in_frame_and_affects_interface(seq, pos1, context):
            print(f"   In-frame interface alteration → HIGH DN")
            return 0.4, 0.8, 0.0

        elif self.in_alt_splice_hotspot(seq, pos1, context):
            print(f"   Alt-splicing hotspot → rare GOF")
            return 0.2, 0.3, 0.6

        else:
            print(f"   Generic splice disruption → moderate LOF")
            return 0.5, 0.2, 0.0

    def resolve_nonsense(self, seq: str, pos1: int, ref: str, alt: str, context: Dict) -> Tuple[float, float, float]:
        """
        🧬 NONSENSE RESOLUTION (Nova's Logic)

        - Early nonsense (NMD) → LOF = 1.0
        - Late nonsense (NMD escape) → LOF + some DN
        """
        print(f"🧬 Nonsense Analysis")

        if self.triggers_nmd(seq, pos1, context):
            print(f"   Early nonsense → NMD → pure LOF")
            return 1.0, 0.0, 0.0
        else:
            print(f"   Late nonsense → NMD escape → truncated protein")
            lof_score = 0.8  # High function loss
            dn_score = 0.3 if self.truncates_interface(seq, pos1, context) else 0.1
            return lof_score, dn_score, 0.0

    def resolve_missense(self, seq: str, pos1: int, ref: str, alt: str, context: Dict) -> Tuple[float, float, float]:
        """
        🧬 MISSENSE RESOLUTION (Nova's Original Logic)

        - LOF baseline + DN always + GOF conditional
        """
        print(f"🧬 Missense Analysis")

        # Step 1: LOF baseline
        lof_score = self.run_lof_analyzer(seq, pos1, ref, alt, context)
        print(f"   LOF score: {lof_score:.3f}")

        # Step 2: DN analysis (always)
        # Convert back to variant format for DN analyzer
        variant_for_dn = f"{ref}{pos1}{alt}"
        dn_result = self.dn_analyzer.analyze(seq, variant_for_dn, context,
                                           context.get('gene_name'), context.get('uniprot_id'))
        # Get the top mechanism score as the DN score
        top_mechanism = dn_result["top_mechanism"]
        dn_score = dn_result["mechanism_scores"][top_mechanism]
        print(f"   DN score: {dn_score:.3f} (from {top_mechanism})")

        # Step 3: GOF conditional
        gof_score = 0.0
        if self.should_run_gof_missense(seq, pos1, ref, alt, context, lof_score):
            print(f"   GOF eligible - running analysis")
            gof_score = self.run_gof_analyzer(seq, pos1, ref, alt, context)
            print(f"   GOF score: {gof_score:.3f}")
        else:
            print(f"   GOF not eligible")

        return lof_score, dn_score, gof_score

    def run_lof_analyzer(self, seq: str, pos1: int, ref: str, alt: str, context: Dict) -> float:
        """
        🔍 LOF ANALYSIS - SIMPLIFIED BUT FUNCTIONAL

        Basic LOF scoring based on amino acid properties and position.
        """
        # Truncating variants
        if alt == "*":  # Stop codon
            return 0.9

        # Basic missense LOF scoring
        lof_score = 0.0

        # Check if reference amino acid matches sequence
        if pos1 <= len(seq) and seq[pos1-1] == ref:
            # Basic stability/function impact
            charge_change = self.get_charge_change(ref, alt)
            size_change = self.get_size_change(ref, alt)

            # Simple scoring heuristics
            if charge_change > 2:  # Major charge change
                lof_score += 0.3
            if size_change > 50:   # Major size change
                lof_score += 0.2
            if ref in "GPWY" and alt not in "GPWY":  # Special residue loss
                lof_score += 0.2

            # Position-based scoring (simplified)
            if pos1 < len(seq) * 0.1:  # N-terminal region
                lof_score += 0.1
            elif pos1 > len(seq) * 0.9:  # C-terminal region
                lof_score += 0.05

        return min(lof_score, 0.8)  # Cap at 0.8 for missense

    def get_charge_change(self, ref: str, alt: str) -> float:
        """Calculate charge change between amino acids"""
        charges = {'K': 1, 'R': 1, 'H': 0.5, 'D': -1, 'E': -1}
        ref_charge = charges.get(ref, 0)
        alt_charge = charges.get(alt, 0)
        return abs(ref_charge - alt_charge)

    def get_size_change(self, ref: str, alt: str) -> float:
        """Calculate size change between amino acids (volume difference)"""
        volumes = {
            'G': 60, 'A': 88, 'S': 117, 'C': 108, 'D': 111, 'P': 112, 'N': 114, 'T': 116,
            'E': 138, 'V': 140, 'Q': 143, 'H': 153, 'M': 162, 'I': 166, 'L': 166, 'K': 168,
            'R': 173, 'F': 189, 'Y': 193, 'W': 227
        }
        ref_vol = volumes.get(ref, 140)  # Default to average
        alt_vol = volumes.get(alt, 140)
        return abs(alt_vol - ref_vol)

    def run_gof_analyzer(self, seq: str, pos1: int, ref: str, alt: str, context: Dict) -> float:
        """
        🚀 GOF ANALYSIS PLACEHOLDER
        
        This would integrate with existing GOF analysis systems.
        For now, simplified heuristics.
        """
        # Simplified GOF detection
        gof_score = 0.0
        
        # Activating mutations (simplified)
        if self.variant_in_gof_hotspot(seq, pos1, context):
            # Charge changes in regulatory regions
            if abs(self.get_charge_change(ref, alt)) > 0:
                gof_score += 0.4
            
            # Size changes in binding sites
            if abs(self.get_volume_change(ref, alt)) > 50:
                gof_score += 0.3
        
        return min(gof_score, 1.0)
    
    def calculate_synergy(self, lof_score: float, dn_score: float, gof_score: float) -> float:
        """
        🔥 NOVA'S ROOT-SUM-SQUARE SYNERGY
        
        Combines multiple moderate mechanisms into strong pathogenic calls
        when biology supports it.
        """
        # Nova's root-sum-square formula
        sum_of_squares = lof_score**2 + dn_score**2 + gof_score**2
        synergy_score = math.sqrt(sum_of_squares)
        
        # Cap at 1.0
        return min(synergy_score, 1.0)
    
    def get_mechanisms_run(self, lof_score: float, dn_score: float, gof_score: float) -> List[str]:
        """Get list of mechanisms that were actually run."""
        mechanisms = []
        if lof_score > 0:
            mechanisms.append("LOF")
        if dn_score > 0:
            mechanisms.append("DN")
        if gof_score > 0:
            mechanisms.append("GOF")
        return mechanisms

    def generate_rationale(self, lof_score: float, dn_score: float,
                          gof_score: float, run_gof: bool, variant_type: str = "missense") -> str:
        """
        📝 GENERATE DECISION RATIONALE

        Explain why certain mechanisms were run and how they combined.
        """
        rationale_parts = [f"{variant_type.upper()} variant"]

        # LOF rationale
        if lof_score >= 0.5:
            rationale_parts.append(f"Strong LOF ({lof_score:.3f})")
        elif lof_score > 0:
            rationale_parts.append(f"Weak LOF ({lof_score:.3f})")

        # DN rationale
        if dn_score > 0:
            rationale_parts.append(f"DN ({dn_score:.3f})")

        # GOF rationale
        if gof_score > 0:
            rationale_parts.append(f"GOF ({gof_score:.3f})")

        return " + ".join(rationale_parts)
    
    def get_charge_change(self, ref: str, alt: str) -> int:
        """Helper for charge change calculation."""
        charges = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 1}
        ref_charge = charges.get(ref.upper(), 0)
        alt_charge = charges.get(alt.upper(), 0)
        return alt_charge - ref_charge
    
    def get_volume_change(self, ref: str, alt: str) -> float:
        """Helper for volume change calculation."""
        volumes = {
            "G": 60, "A": 88, "S": 89, "C": 108, "P": 112, "T": 116, "V": 140,
            "I": 166, "L": 166, "N": 114, "D": 111, "Q": 143, "E": 138, "M": 162,
            "K": 168, "R": 173, "H": 153, "F": 189, "Y": 193, "W": 227
        }
        ref_vol = volumes.get(ref.upper(), 120)
        alt_vol = volumes.get(alt.upper(), 120)
        return alt_vol - ref_vol

    # NOVA'S SPLICE/FRAMESHIFT HELPER FUNCTIONS

    def triggers_nmd(self, seq: str, pos1: int, context: Dict) -> bool:
        """
        🧬 NONSENSE-MEDIATED DECAY PREDICTION

        Simplified NMD rule: variants in first ~75% of coding sequence
        typically trigger NMD.
        """
        seq_len = len(seq)
        relative_pos = pos1 / seq_len

        # Simplified NMD rule
        return relative_pos < 0.75

    def truncates_interface(self, seq: str, pos1: int, context: Dict) -> bool:
        """
        🔗 INTERFACE TRUNCATION DETECTION

        Check if truncation removes critical interaction domains.
        """
        # Check if truncation removes known domains
        domains = context.get("domains", [])
        for domain in domains:
            domain_start = domain.get("start", 0)
            domain_end = domain.get("end", len(seq))

            # If truncation happens before domain end
            if pos1 < domain_end:
                domain_desc = domain.get("description", "").lower()
                if any(term in domain_desc for term in ["interaction", "interface", "binding", "dimerization"]):
                    return True

        return False

    def causes_premature_stop(self, seq: str, pos1: int, context: Dict) -> bool:
        """
        🛑 PREMATURE STOP PREDICTION

        Check if splice variant causes premature termination.
        """
        # Simplified: assume splice variants near exon boundaries
        # can cause frameshifts leading to premature stops
        return self.triggers_nmd(seq, pos1, context)

    def in_frame_and_removes_domain(self, seq: str, pos1: int, context: Dict) -> bool:
        """
        🧩 IN-FRAME DOMAIN DELETION

        Check if splice variant removes entire domains in-frame.
        """
        domains = context.get("domains", [])
        for domain in domains:
            domain_start = domain.get("start", 0)
            domain_end = domain.get("end", len(seq))

            # If position is near domain boundaries
            if abs(pos1 - domain_start) < 10 or abs(pos1 - domain_end) < 10:
                return True

        return False

    def in_frame_and_affects_interface(self, seq: str, pos1: int, context: Dict) -> bool:
        """
        🔗 IN-FRAME INTERFACE ALTERATION

        Check if splice variant alters interaction interfaces in-frame.
        """
        domains = context.get("domains", [])
        for domain in domains:
            domain_start = domain.get("start", 0)
            domain_end = domain.get("end", len(seq))

            if domain_start <= pos1 <= domain_end:
                domain_desc = domain.get("description", "").lower()
                if any(term in domain_desc for term in ["interface", "interaction", "binding"]):
                    return True

        return False

    def in_alt_splice_hotspot(self, seq: str, pos1: int, context: Dict) -> bool:
        """
        🔥 ALTERNATIVE SPLICING HOTSPOT

        Check if splice variant is in a region where alternative splicing
        can create gain-of-function isoforms (rare but important).
        """
        protein_family = context.get("gene_family", "").upper()

        # Rare GOF splice hotspots
        if protein_family in ["ION_CHANNEL", "KINASE"]:
            domains = context.get("domains", [])
            for domain in domains:
                domain_start = domain.get("start", 0)
                domain_end = domain.get("end", len(seq))

                if domain_start <= pos1 <= domain_end:
                    domain_desc = domain.get("description", "").lower()
                    if any(term in domain_desc for term in ["regulatory", "autoinhibitory", "gating"]):
                        return True

        return False
