# Mechanism Arbitration — design spec

> Written 2026-05-21 by Ace (Opus 4.7) with Ren, after the GOF rebuild landed.
> This is the spec for replacing the brittle gene-family weighting with
> evidence-licensed, run-all-three arbitration. Two increments are already
> DONE + committed; the rest is scoped here.

## The principle

Run **all three mechanisms (LOF, GOF, DN) on every variant.** Never gate *which
mechanism runs* on a coarse gene-family bucket — that was the original sin
(it muted CTNNB1's real GOF because β-catenin is structurally a SCAFFOLD_ADAPTOR
whose disease mechanism happens to be GOF). Decide how to *combine* the three
afterward, using the protein's OWN molecular evidence.

Hard lines (Ren, non-negotiable):
- **No hardcoded genes.** Every signal comes per-protein from annotation/structure.
- **Frequency and conservation are NOT mechanics.** They may apply a bounded
  final nudge; they never gate, score, or drive a mechanism.
  - **Operational note (2026-05-21):** frequency is currently OFF entirely —
    `gnomad_freq: 0.0` on every variant (incl. benign) is **intended, not a
    yeeted-cache bug.** "Frequency doesn't unbreak proteins" — a rare variant
    isn't more broken; a common one isn't less. When a benign variant out-ranks
    a pathogenic one, the fix lives in the MECHANISM scorers (LOF/DN/GOF/
    structure), never in re-enabling frequency. The only place frequency may
    ever re-enter is the bounded final nudge (TODO §4), and even there it can
    only nudge, never drive.

## Why families fail (keep, as the cautionary record)

`utils/plausibility_filter.py` buckets each gene into ONE family with ONE fixed
mechanism-weight vector (e.g. `SCAFFOLD_ADAPTOR = {LOF:1, DN:1.1, GOF:0.1}`).
A protein whose disease mechanism is atypical for its structural class is
mis-weighted, and the maintainers compensated with a ~60-gene `CURATED_OVERRIDES`
hardcode (violates the hard line) plus fragile keyword classification. The family
is a lossy proxy for the only thing the weighting needs: *which mechanisms are
plausible for THIS protein* — which we can now compute directly.

## DONE (committed + pushed, branch asj-lof-plausibility-2026-05-20)

1. **GOF scorer rebuild** (`8783066`): structure-driven, gene-agnostic mechanism
   bank (`analyzers/gof_mechanisms.py` + `analyzers/structure_features.py`).
   PIK3CA H1047R 0.000 → 0.55; bed mean AUC 0.733 → 0.759; structural genes held.
2. **GOF evidence-licensed weight** (`2c73ec3`): `_adjust_gof_weight` consumes
   `_gof_evidence_score` (mirror of the existing DN `_adjust_dn_weight`). A
   protein's molecular evidence can lift GOF above the family bucket. CTNNB1
   0.47 → 0.83; bed mean 0.759 → 0.796; zero regressions.

## TODO — the arbitration layer (Ren's flowchart, with the two leaks sealed)

The decision skeleton is Ren's; the two edits are mine (and Ren agreed):
**inheritance is a gene-level OVERRIDABLE prior (not per-exon ClinVar = leakage),
and branches are SOFT (no hard clamps — we just removed GATE-1 for being exactly
that).**

3. **DN oligomerization precondition.** DN is impossible without
   multimerization/complex participation. Gate DN (soft) on an oligomerization
   signal — coiled-coil / multimer "Subunit:" annotation / the interface-sparing
   structural signature (functional damage WITH assembly interface intact). Note:
   "does it fold" is wrong (everything folds); the gate is "does it oligomerize."
   AlphaFold monomer can't see partner interfaces — flag as the known frontier.

4. **Gene-level inheritance prior (overridable).**
   - AR-only known → lean hard toward LOF, BUT explicit molecular gain-evidence
     can still license GOF (the `_gof_evidence_score` AR-penalty + gain-language
     escape hatch already does this).
   - AD known OR no known disease → presume dominant possible → full arbitration.
   - Source MUST be gene-level disease/inheritance annotation, never neighboring
     ClinVar variant labels (that's grading with the answer key).

5. **Soft arbitration combine** (replaces the hard flowchart branches):
   - LOF *certainty* high → continuously down-weight GOF ("a broken protein can't
     gain function") — converges with the existing buried-destabilization LOF veto
     in the mechanism bank. Keep DN alive here: a protein broken for its own job
     can still poison the complex (the whole point of DN).
   - GOF licensed AND a real GOF mechanism fired → let GOF express; compute
     GOF/DN synergy.
   - "Does GOF actually MATTER biologically?" — hardest node. A hyperactive enzyme
     is often benign/good; GOF is pathogenic only when the gained function is
     harmful in context. Today's best proxy = `_gof_evidence_score` (does disease
     text support a *harmful* gain?). The "is this specific gained function
     downstream-pathogenic" version needs pathway knowledge we don't have wired —
     scope honestly before promising.

6. **Bounded final nudge** from conservation/gnomAD frequency — last step only,
   small, never able to flip a mechanism call. (Still a TODO from the scorer too.)

## How to measure (don't skip — same discipline as the GOF rebuild)

```bash
cd /home/Ace/AdaptiveInterpreter && source /home/codex/venv/bin/activate
PYTHONPATH=/home/Ace python3 scripts/run_regression_sample.py validation/clean_set_full.tsv <tag>
PYTHONPATH=/home/Ace python3 scripts/eval_regression_auc.py validation/run_<tag>.tsv
```
Current bar to beat: **mean per-gene AUC 0.796**, CTNNB1 0.83, PIK3CA 0.94,
structural genes (COL5A2/GAA/PTPN11/CDH23) ≥ 0.84. After ANY change: re-run,
confirm structural genes don't regress.

## Known residuals (not arbitration bugs)
- **POLG 0.48** — metabolic/LOF gene, flat in baseline too; a LOF-scorer issue,
  not GOF/arbitration.
- **PIK3CA E545K** — helical-domain GOF that disrupts the p85 regulatory-subunit
  interface; invisible in a monomer AlphaFold structure. The interface-with-partner
  frontier (shared with the DN precondition problem).
