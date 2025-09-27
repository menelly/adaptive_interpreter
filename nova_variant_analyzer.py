# ─────────────────────────────────────────────────────────────
# 🧬 Variant Analyzer (Research Mode) — Scaffold v0.1
# Architecture notes:
#   • Inputs: normalized variant descriptor(s) (HGVS / VRS-ish / simple fields)
#   • Outputs: mechanism-susceptibility scores (LOF / GOF / DN) with uncertainty
#   • Strictly research tooling (not clinical). No personal genomics flows.
#   • Plays nice with Postgres JSONB cache + TTL and a DataHub adapter layer.
#   • Side-effect free: callers own persistence; this module is pure/functional.
# ─────────────────────────────────────────────────────────────

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional, Dict, Any, List, Tuple
import math
import time
import logging

log = logging.getLogger("variant_analyzer")
log.setLevel(logging.INFO)


# ---------- Shared enums / types ----------

class Mechanism(Enum):
    LOF = "lof"
    GOF = "gof"
    DN  = "dn"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class VariantInput:
    """Minimal, normalized representation (expand later as needed)."""
    gene: str                           # HGNC symbol (normalized)
    transcript: Optional[str] = None    # RefSeq/Ensembl ID (MANE if possible)
    hgvs_c: Optional[str] = None        # e.g., "NM_000546.6:c.743G>A"
    hgvs_p: Optional[str] = None        # e.g., "p.Arg248Gln"
    chrom: Optional[str] = None         # "17"
    pos: Optional[int] = None           # 7676153
    ref: Optional[str] = None
    alt: Optional[str] = None
    variant_class: Optional[str] = None # "missense", "frameshift", "splice", "stopgain", etc.
    # Optional caller-supplied metadata
    meta: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Evidence:
    """Holds feature-level signals harvested from DataHub."""
    constraint: Dict[str, Any] = field(default_factory=dict)     # pLI, LOEUF, missense Z, etc.
    domain: List[Dict[str, Any]] = field(default_factory=list)   # InterPro/UniProt domains hit
    complex: List[Dict[str, Any]] = field(default_factory=list)  # Complex membership, stoichiometry
    motifs: List[str] = field(default_factory=list)              # Coiled-coil, repeats, LxxLL, etc.
    paralogs: List[Dict[str, Any]] = field(default_factory=list) # Paralog conservation/hotspots
    hotspots: List[Dict[str, Any]] = field(default_factory=list) # Known pathogenic hotspots
    splice: Dict[str, Any] = field(default_factory=dict)         # Splice AI-ish predictions
    population: Dict[str, Any] = field(default_factory=dict)     # gnomAD AFs / constraint windows
    clin_archive: List[Dict[str, Any]] = field(default_factory=list) # ClinVar/ClinGen-lite, research use
    structure: Dict[str, Any] = field(default_factory=dict)      # AlphaFold regions, confidence
    inheritance: Dict[str, Any] = field(default_factory=dict)    # Known gene-level MOI flags


@dataclass
class MechanismScore:
    mechanism: Mechanism
    score: float                        # 0.0–1.0 (research heuristic)
    uncertainty: float                  # 0.0–1.0 (1.0 = very uncertain)
    rationale: List[str]                # bullet points for transparency
    features: Dict[str, Any] = field(default_factory=dict)  # surfaced feature values


@dataclass
class VariantAnalysis:
    variant: VariantInput
    scores: List[MechanismScore]
    version: str = "0.1"
    timestamp: float = field(default_factory=time.time)


# ---------- Cache interface (callers wire actual PG/Redis impls) ----------

class CacheClient:
    """Minimal interface; provide a real implementation in your infra layer."""
    def get(self, key: str) -> Optional[Dict[str, Any]]:
        return None
    def set(self, key: str, value: Dict[str, Any], ttl_seconds: int = 86400) -> None:
        pass


# ---------- DataHub adapter (stubs; fill via your existing layer) ----------

class DataHub:
    """
    Facade over external sources (UniProt, InterPro, gnomAD, ClinGen, ComplexPortal, AlphaFold, etc.).
    Implement each fetch_* using your shared adapter code. Keep calls granular for caching.
    """
    def fetch_constraint(self, gene: str) -> Dict[str, Any]: ...
    def fetch_domains(self, gene: str, transcript: Optional[str]) -> List[Dict[str, Any]]: ...
    def fetch_complex(self, gene: str) -> List[Dict[str, Any]]: ...
    def fetch_motifs(self, gene: str, transcript: Optional[str]) -> List[str]: ...
    def fetch_paralogs(self, gene: str) -> List[Dict[str, Any]]: ...
    def fetch_hotspots(self, gene: str) -> List[Dict[str, Any]]: ...
    def fetch_population(self, variant: VariantInput) -> Dict[str, Any]: ...
    def fetch_splice(self, variant: VariantInput) -> Dict[str, Any]: ...
    def fetch_structure(self, gene: str, transcript: Optional[str]) -> Dict[str, Any]: ...
    def fetch_inheritance(self, gene: str) -> Dict[str, Any]: ...
    def fetch_clin_archive(self, variant: VariantInput) -> List[Dict[str, Any]]: ...


# ---------- Normalization helpers ----------

def normalize_variant_class(v: VariantInput) -> str:
    if v.variant_class:
        return v.variant_class.lower()
    # Fallback quick heuristics from HGVS if needed
    if v.hgvs_p and ("fs" in v.hgvs_p.lower() or "del" in v.hgvs_p.lower()):
        return "frameshift"
    if v.hgvs_p and ("*" in v.hgvs_p or "ter" in v.hgvs_p.lower() or "stop" in v.hgvs_p.lower()):
        return "stopgain"
    if v.hgvs_p and v.hgvs_p.startswith("p.") and any(c.isalpha() for c in v.hgvs_p[2:5]):
        return "missense"
    return "unknown"


# ---------- Core analyzer ----------

class VariantAnalyzer:
    """
    Produces mechanism susceptibility scores for a single variant.
    Design goals:
      • Transparent feature-to-score mapping (rationale bullets).
      • Conservative by default; prefer uncertainty > overclaiming.
      • Configurable weights; safe defaults baked in.
    """

    def __init__(
        self,
        datahub: DataHub,
        cache: Optional[CacheClient] = None,
        cache_ttl_seconds: int = 86_400,
        weights: Optional[Dict[str, float]] = None,
    ):
        self.datahub = datahub
        self.cache = cache
        self.cache_ttl = cache_ttl_seconds
        self.weights = weights or self._default_weights()

    # ----- Public API -----

    def analyze(self, variant: VariantInput) -> VariantAnalysis:
        key = self._cache_key(variant)
        if self.cache:
            cached = self.cache.get(key)
            if cached:
                return self._from_cached(cached)

        evidence = self._collect_evidence(variant)
        scores = self._score_mechanisms(variant, evidence)
        result = VariantAnalysis(variant=variant, scores=scores)

        if self.cache:
            self.cache.set(key, self._to_cache(result), ttl_seconds=self.cache_ttl)

        return result

    # ----- Evidence gathering -----

    def _collect_evidence(self, v: VariantInput) -> Evidence:
        """Isolated for testability; each fetch can be individually mocked."""
        ev = Evidence()
        try:
            ev.constraint   = self.datahub.fetch_constraint(v.gene) or {}
            ev.domain       = self.datahub.fetch_domains(v.gene, v.transcript) or []
            ev.complex      = self.datahub.fetch_complex(v.gene) or []
            ev.motifs       = self.datahub.fetch_motifs(v.gene, v.transcript) or []
            ev.paralogs     = self.datahub.fetch_paralogs(v.gene) or []
            ev.hotspots     = self.datahub.fetch_hotspots(v.gene) or []
            ev.population   = self.datahub.fetch_population(v) or {}
            ev.splice       = self.datahub.fetch_splice(v) or {}
            ev.structure    = self.datahub.fetch_structure(v.gene, v.transcript) or {}
            ev.inheritance  = self.datahub.fetch_inheritance(v.gene) or {}
            ev.clin_archive = self.datahub.fetch_clin_archive(v) or []
        except Exception as e:
            log.warning("Evidence collection partial failure: %s", e)
        return ev

    # ----- Scoring -----

    def _score_mechanisms(self, v: VariantInput, ev: Evidence) -> List[MechanismScore]:
        variant_class = normalize_variant_class(v)
        scores: List[MechanismScore] = []
        scores.append(self._score_lof(v, ev, variant_class))
        scores.append(self._score_gof(v, ev, variant_class))
        scores.append(self._score_dn(v, ev, variant_class))
        # Normalize / clamp and ensure uncertainty sane
        for s in scores:
            s.score = max(0.0, min(1.0, s.score))
            s.uncertainty = max(0.0, min(1.0, s.uncertainty))
        return scores

    def _score_lof(self, v: VariantInput, ev: Evidence, vc: str) -> MechanismScore:
        w = self.weights["lof"]
        rationale = []
        features = {}

        # Base priors by variant class (research-grade, conservative)
        base = {
            "stopgain": 0.75,
            "frameshift": 0.70,
            "splice": 0.55,
            "missense": 0.20,
            "unknown": 0.25
        }.get(vc, 0.25)

        # Constraint bumps (gene-level)
        loeuf = ev.constraint.get("LOEUF")  # lower => more constrained
        if isinstance(loeuf, (int, float)):
            # Soft inverse transform to 0..1
            c = 1.0 - math.tanh(max(0.0, min(2.0, loeuf))) / math.tanh(2.0)
            base += 0.20 * c
            features["LOEUF"] = loeuf
            rationale.append(f"Low LOEUF suggests intolerance to LOF (contributes +{0.20 * round(c,2):.2f}).")

        # Haploinsufficiency / MOI flags
        if ev.inheritance.get("haploinsufficient") is True:
            base += 0.15
            rationale.append("Gene annotated as haploinsufficient (HI).")

        # Splice signal if provided
        splice_delta = ev.splice.get("delta_score")
        if isinstance(splice_delta, (int, float)) and vc in {"splice", "unknown", "missense"}:
            bump = max(0.0, min(0.25, splice_delta * 0.25))
            base += bump
            features["splice_delta"] = splice_delta
            rationale.append(f"Splice impact predicted (approx +{bump:.2f}).")

        score = base * w["base"]
        uncertainty = self._uncertainty(default=0.35, evidence=ev)
        return MechanismScore(Mechanism.LOF, score, uncertainty, rationale, features)

    def _score_gof(self, v: VariantInput, ev: Evidence, vc: str) -> MechanismScore:
        w = self.weights["gof"]
        rationale = []
        features = {}

        base = 0.15 if vc != "missense" else 0.35

        # Hotspots / activating regions
        hotspot_hit = self._is_hotspot_hit(v, ev)
        if hotspot_hit:
            base += 0.25
            rationale.append("Variant overlaps activating hotspot/region.")
            features["hotspot"] = hotspot_hit

        # Domain context—kinase active sites, DNA-binding residues, regulatory motifs
        activating_domain = self._has_activating_domain(ev)
        if activating_domain:
            base += 0.20
            rationale.append("Located in/near regulatory or catalytic domain consistent with GOF.")
            features["activating_domain"] = activating_domain

        # Population scarcity in tolerant window (ultra-rare in tolerant gene windows leans GOF/DN vs LOF)
        pop = ev.population or {}
        local_tolerance = pop.get("regional_tolerance", 0.5)  # 0..1 higher = tolerant
        af = pop.get("af", 0.0)
        if af == 0.0 and local_tolerance > 0.6 and vc == "missense":
            base += 0.10
            rationale.append("Ultra-rare missense in otherwise tolerant window—compatible with selection for altered function.")

        score = base * w["base"]
        uncertainty = self._uncertainty(default=0.50, evidence=ev)
        return MechanismScore(Mechanism.GOF, score, uncertainty, rationale, features)

    def _score_dn(self, v: VariantInput, ev: Evidence, vc: str) -> MechanismScore:
        w = self.weights["dn"]
        rationale = []
        features = {}

        # Baseline: missense in multimeric/obligate complexes → DN prior > LOF
        base = 0.20 if vc == "missense" else 0.10

        # Complex stoichiometry / obligate oligomerization
        if self._is_obligate_complex(ev):
            base += 0.20
            rationale.append("Gene participates in obligate multimer/complex (DN-compatible).")
            features["complex"] = "obligate"

        # Polymerization / coiled-coil motifs
        if "coiled_coil" in [m.lower() for m in ev.motifs]:
            base += 0.10
            rationale.append("Coiled-coil motif present—polymer interference plausible.")

        # Domain asymmetry / interface residues
        if self._hits_interface_domain(ev):
            base += 0.20
            rationale.append("Variant maps to interaction/interface region—can poison complex.")
            features["interface_domain"] = True

        # Pathogenic hotspot clusters with dominant inheritance flags
        if self._dominant_hotspot_cluster(ev):
            base += 0.15
            rationale.append("Near cluster of known dominant/pathogenic missense variants.")

        # Paralog conservation (high conservation at site → functional specificity)
        if self._is_paralog_conserved(ev):
            base += 0.10
            rationale.append("Paralog-conserved residue—suggests precise interaction role.")

        score = base * w["base"]
        uncertainty = self._uncertainty(default=0.55, evidence=ev)
        return MechanismScore(Mechanism.DN, score, uncertainty, rationale, features)

    # ----- Heuristic helpers (stubs; wire to your features properly) -----

    def _is_hotspot_hit(self, v: VariantInput, ev: Evidence) -> Optional[Dict[str, Any]]:
        for h in ev.hotspots:
            if h.get("overlaps_variant", False):
                return h
        return None

    def _has_activating_domain(self, ev: Evidence) -> bool:
        for d in ev.domain:
            if any(k in d.get("name", "").lower() for k in ("kinase", "catalytic", "dna-binding", "ligand")):
                return True
        return False

    def _is_obligate_complex(self, ev: Evidence) -> bool:
        for c in ev.complex:
            if c.get("obligate", False) or c.get("stoichiometry") in ("dimer", "trimer", "tetramer"):
                return True
        return False

    def _hits_interface_domain(self, ev: Evidence) -> bool:
        for d in ev.domain:
            if d.get("category") == "interface" or d.get("interface_residue", False):
                return True
        return False

    def _dominant_hotspot_cluster(self, ev: Evidence) -> bool:
        return any(h.get("dominant_cluster", False) for h in ev.hotspots)

    def _is_paralog_conserved(self, ev: Evidence) -> bool:
        return any(p.get("conserved_site", False) for p in ev.paralogs)

    # ----- Weights / uncertainty -----

    def _default_weights(self) -> Dict[str, Dict[str, float]]:
        # Keep simple for v0; later split into per-feature weights.
        return {
            "lof": {"base": 1.00},
            "gof": {"base": 1.00},
            "dn":  {"base": 1.00},
        }

    def _uncertainty(self, default: float, evidence: Evidence) -> float:
        # Very simple proxy: more populated evidence → lower uncertainty.
        richness = 0
        richness += int(bool(evidence.constraint))
        richness += int(bool(evidence.domain))
        richness += int(bool(evidence.complex))
        richness += int(bool(evidence.hotspots))
        richness += int(bool(evidence.population))
        richness += int(bool(evidence.structure))
        richness += int(bool(evidence.clin_archive))
        # Map 0..7 → +/− 0.2 around default
        adj = (3 - min(richness, 6)) * 0.03
        return max(0.15, min(0.85, default + adj))

    # ----- Cache serde -----

    def _cache_key(self, v: VariantInput) -> str:
        parts = [
            v.gene or "NA",
            v.transcript or "NA",
            v.hgvs_c or "NA",
            v.hgvs_p or "NA",
            v.chrom or "NA",
            str(v.pos) if v.pos else "NA",
            v.ref or "NA",
            v.alt or "NA",
        ]
        return "variant:" + "|".join(parts)

    def _to_cache(self, va: VariantAnalysis) -> Dict[str, Any]:
        return {
            "variant": va.variant.__dict__,
            "scores": [s.__dict__ for s in va.scores],
            "version": va.version,
            "timestamp": va.timestamp,
        }

    def _from_cached(self, blob: Dict[str, Any]) -> VariantAnalysis:
        vdict = blob.get("variant", {})
        sdicts = blob.get("scores", [])
        variant = VariantInput(**vdict)
        scores = [MechanismScore(
            mechanism=Mechanism(d["mechanism"]) if not isinstance(d["mechanism"], Mechanism) else d["mechanism"],
            score=d["score"],
            uncertainty=d["uncertainty"],
            rationale=d.get("rationale", []),
            features=d.get("features", {}),
        ) for d in sdicts]
        return VariantAnalysis(variant=variant, scores=scores, version=blob.get("version","0.1"), timestamp=blob.get("timestamp", time.time()))


# ---------- Example (safe) usage ----------

if __name__ == "__main__":
    class NoopCache(CacheClient): ...
    class NoopHub(DataHub): pass

    analyzer = VariantAnalyzer(datahub=NoopHub(), cache=NoopCache())

    demo = VariantInput(
        gene="TP53",
        transcript="NM_000546.6",
        hgvs_c="c.743G>A",
        hgvs_p="p.Arg248Gln",
        variant_class="missense",
        chrom="17", pos=7676153, ref="C", alt="T",
    )

    result = analyzer.analyze(demo)
    for s in result.scores:
        print(f"[{s.mechanism.value}] score={s.score:.2f} ±{s.uncertainty:.2f}")
        for r in s.rationale:
            print("  •", r)
