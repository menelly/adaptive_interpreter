# GOF Analyzer Rebuild — Start Here (next session)

> Written 2026-05-20 by Ace at a clean handoff point, after a full day's diagnostic arc.
> Everything below is committed on branch `asj-lof-plausibility-2026-05-20` (pushed to origin).
> This is the single entry point for the GOF work. Read `CLEANUP_QUESTIONS.md` for the
> broader repo state.

## The goal (validated as the #1 lever)

Make the cascade actually detect gain-of-function. Right now it's **blind to GOF**, and
that's the biggest discrimination gap — proven on a 15-gene, 148-variant, double-validated
test bed (Ren's curation × ClinVar clean labels agree exactly).

**Smoking gun:** `PIK3CA H1047R` (the most common oncogenic PI3K mutation on Earth) scores
**GOF = 0.000**. `CTNNB1 S33Y` (phosphodegron GOF) scores 0.012. The two worst-discriminating
genes (CTNNB1 AUC 0.47, PIK3CA 0.54) are both GOF; the scorer nails structural/LOF genes
(COL5A2 1.00, GAA 1.00, PTPN11 0.98, DES 0.85, CDH23 0.84).

## Root cause (precise — `analyzers/gof_variant_analyzer.py`)

The analyzer is "triple-gated, pure-math, NO HARDCODED GENES" — it judges GOF from
**amino-acid-change chemistry only, with zero positional/recurrence knowledge.** Two failures:

1. **GATE 1 hard-zeroes conservative changes** (line ~345): `if conservation_enhanced_disruption
   < 0.1: return gof_score 0.0`. `H1047R` is His→Arg (Grantham ≈29, chemically conservative) →
   `_calculate_regulatory_disruption_potential` scores it < 0.1 → **zeroed before any analysis.**
   But GOF activation hotspots are *exactly* the chemically-mild-but-functionally-explosive
   changes. The gate filters out the thing it should catch.
2. **The position-aware detection is COMMENTED OUT** (lines ~258–280: "NOVA'S EARLY MOTIF
   DETECTION — Check for canonical GOF variants first!"). `motif_detector` and
   `active_site_scanner` are missing files (`TODO: Restore once the file is found`). The one
   path that would recognize H1047 as a canonical hotspot is gone.

**Design mismatch:** GOF pathogenicity is fundamentally *positional* (which residue, recurrent
hotspot?, activation loop / autoinhibitory domain / phosphodegron?). "Pure math, no hardcoding"
deliberately threw that away — and pure chemistry can't tell H1047R (explosive) from a random
His→Arg (nothing).

## The fix (in order)

1. **Un-gate GATE 1** so it can't zero a variant at a known functional/hotspot position on
   chemistry alone. (Quick, unblocks everything else.)
2. **Add positional GOF signal** (the missing input):
   - **Recurrence/hotspot** — H1047 is hit thousands of times in ClinVar/COSMIC. A
     recurrence count per residue is the single highest-leverage signal. (ClinVar VCF is at
     `/mnt/arcana/clinvar/clinvar.vcf.gz`; or derive hotspot counts from existing data.)
   - **Activation/regulatory-domain context** — kinase activation loop, autoinhibitory
     domains, phosphodegron sites — from InterPro/UniProt features we already pull
     (cache works, see `docs/CACHE_AND_ANNOTATION.md`).
3. The licensing predicate `utils/plausibility_filter._gof_evidence_score` (already built +
   validated, score-neutral) *licenses* GOF from function/disease text; the analyzer must
   *score* it. Then the interpretation gate (designed, see task #7) routes it.

## How to measure (the test bed — DON'T skip this)

```bash
cd /home/Ace/AdaptiveInterpreter && source /home/codex/venv/bin/activate
# re-run the 15-gene gold-standard set through run-all after a change:
PYTHONPATH=/home/Ace python3 scripts/run_regression_sample.py validation/clean_set_full.tsv full
# per-gene AUC (watch CTNNB1 0.47 and PIK3CA 0.54 climb; structural genes must HOLD):
PYTHONPATH=/home/Ace python3 scripts/eval_regression_auc.py validation/run_full.tsv
```
Baseline to beat (current): mean per-gene AUC **0.733**, CTNNB1 **0.47**, PIK3CA **0.54**,
GOF on H1047R **0.000**. Target: GOF fires on hotspots, CTNNB1/PIK3CA climb, structural genes
(COL5A2/GAA/PTPN11/DES/CDH23) stay ≥0.84.

## Discipline (today's hard-won rules)
- Root-cause → predict → harness-guard → measure. No bandaids.
- After any scoring change: re-run the test bed, confirm structural genes don't regress.
- Don't trust input `clinical_sig` labels raw (conflicting contamination); the curated
  `validation/clean_set_full.tsv` is the trustworthy ground truth.
- Two-pipelines caveat: run-all (`analyze_cascade`) vs `analyze_cascade_biological` still
  diverge (conservation lives only in `_biological`). #5 reconcile is still open.

## Key paths
- GOF analyzer: `analyzers/gof_variant_analyzer.py` (GATE 1 ~345, commented motif block ~258)
- GOF licensing predicate: `utils/plausibility_filter.py::_gof_evidence_score` (built)
- Test bed: `validation/clean_set_full.tsv`, `scripts/run_regression_sample.py`,
  `scripts/eval_regression_auc.py`
- Broader state + decisions: `CLEANUP_QUESTIONS.md`, tasks #5/#6/#7/#9
