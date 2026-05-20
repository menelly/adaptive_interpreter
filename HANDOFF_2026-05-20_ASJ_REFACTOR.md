# Handoff: ASJ→LOF Refactor + Run-All Architecture
**Date:** 2026-05-20 (9:30am → 12:45pm)
**For:** Next-Ace + Ren
**Branch:** `asj-lof-plausibility-2026-05-20` (commits `0e543a3`, `0f60c10`)

---

## What We Did Today

Massive cascade analyzer refactor with Ren. Started from "DN scores are too high for LOF genes" and peeled the onion down to root causes.

### Architecture Changes

1. **Run-all, no gates.** DN, LOF, GOF all run every time. Mechanisms are protein physics — they always score. Raw scores preserved. Plausibility decides what to believe, not what to run.

2. **ASJ moved to LOF natively.** Active site jamming was in `nova_dn/mechanisms.py` and scored as DN. It's LOF — "the enzyme broke" not "the enzyme poisoned the complex." Now scored inside `lof_analyzer.py` as a sub-mechanism, weighted by conservation. DN's poison signal comes from `interface_poisoning` separately (no double-counting).

3. **UniProt inheritance drives DN plausibility.** Disease comments extracted during annotation (`universal_protein_annotator.py`). AR-only genes (HEXA) get DN halved. AD genes (ATP5F1A, MYH7) keep full DN. Mixed (ATP5F1A with both AD+AR conditions) → keeps DN.

4. **GO terms wired into DN filters.** Both v1 and v2 mechanism filters now search `function_text + GO terms`. Added enzyme class suffixes (hydrolase, transferase, hexosaminidase, etc.).

5. **InterPro domains feed LOF.** `lof_analyzer._get_domain_context` supplements old predicted domains with InterPro functional domains. "Beta-hexosaminidase, catalytic domain" → 1.3× multiplier instead of generic 1.2×.

6. **Signal/transit peptide clamp.** Variants before mature chain start get 0.3× weight. Propeptides exempted (collagen C-propeptide is functional before cleavage). ATP5F1A P9S (B): 0.495 → 0.124.

7. **Stale caches nuked.** 1,160 of 2,139 `_domains.json` files were pre-GO-term era. Deleted. They re-fetch lazily with GO terms, active sites, and inheritance on next access.

### What's Still Wrong

- **Family classifier is wrong for many genes.** COL1A1→TUMOR_SUPPRESSOR, MFN2→ONCOGENE, TFG→ONCOGENE. Ren's intuition: maybe ditch families for DN/LOF routing entirely, keep only for GOF suppression. UniProt inheritance is more reliable.

- **GOF hallucination.** GOF analyzer scores high on generic amino acid changes without biological basis. Family GOF=0.0 catches it for metabolic enzymes but misses genes with wrong family classification. Need: inheritance-based GOF suppression (AR-only gene + GOF = nonsensical) or GOF analyzer improvements.

- **LOF base scores are still too low** for some pathogenics (HEXA W420C=0.527 when it should be higher). The base LOF formula is `(stability*0.3 + structural*0.2 + functional*0.2 + ASJ*0.3) / 0.7` — the weights may need recalibration now that ASJ is in the mix.

- **Normalization across gene families.** Collagen Gly→X scores 1.0 while HEXA W→C scores 0.527. Both are pathogenic. The scoring ranges are family-dependent. Ren proposed per-family normalization using the 24k+ calibrated ClinVar P/LP variants.

- **Positional inheritance (Ren's big idea, not yet implemented).** Instead of gene-level "is this AR or AD?", check ClinVar variants in the same region/exon/domain. Same gene can have AD regions and AR regions. Data source: ClinVar VCF at `/mnt/arcana/clinvar/clinvar.vcf.gz` (Sept 2025, has ORIGIN field with de-novo flag = AD signal). Or build from our own cascade outputs.

- **Proline/cysteine/glycine ML bumps** may be adding noise. Ren wanted to test stripping them to see if ROC improves. Not yet tested.

- **BMPR1A has no UniProt inheritance data** (empty). Several genes have this gap. May need OMIM fallback or ClinVar VCF ORIGIN field.

### Files Modified

```
analyzers/cascade_analyzer.py          — run-all architecture, plausibility pipeline
analyzers/lof_analyzer.py              — ASJ as LOF sub-mechanism, mature chain clamp, InterPro domains
data_processing/universal_protein_annotator.py — disease/inheritance extraction
nova_dn/dn_mechanism_filter.py         — GO terms in searchable text
nova_dn/dn_mechanism_filter_v2.py      — GO terms, dn_indicators, enzyme suffixes
utils/plausibility_filter.py           — inheritance_patterns param, DN evidence scoring
README.md                              — updated numbers
```

### Three Babies Training

v4/v5 pilot completed while we worked. Check `logs/v4v5_pilot.log`. PEFT fallback merge worked for Gemma (Unsloth merge bug caught and handled). Need to check if all 4 runs succeeded (gemma v4, qwen v4, gemma v5, qwen v5). Eval completions next, then HOLD for Ren on judge scoring (API cost).

### Test Commands

```bash
# Quick test on HEXA + ATP5F1A
source /home/codex/venv/bin/activate
PYTHONPATH=/home/Ace python3 -c "
from AdaptiveInterpreter.analyzers.cascade_analyzer import CascadeAnalyzer
a = CascadeAnalyzer()
for g,v in [('HEXA','p.W420C'),('ATP5F1A','p.R207H')]:
    r = a.analyze_cascade(g, v)
    s = r['scores']
    print(f'{g} {v}: DN={s[\"DN\"]:.3f} LOF={s[\"LOF\"]:.3f} Final={r[\"final_score\"]:.3f}({r[\"final_classification\"]})')
"
```

---

*Ren's brain is in genetics mode and they're about to break more things. Save state appreciated.*

*— Ace, 2026-05-20 ~12:45pm, after the best genetics session in months*

---

## Afternoon Update (2026-05-20 ~1:30pm) — Cache Infrastructure Fixed

Picked MEFV as a GOF stress-test and it exposed the annotation cache was rotten. Fixed it (commit `e4b9901`, same branch):

- **Two stacked cache bugs** killed inheritance/disease data: (1) `cache_dir` was a *relative* path → forked into 3 dirs by launch cwd; (2) `get_uniprot_features` returned stale caches verbatim, never refetching, so pre-2026-05 caches had empty `inheritance_patterns`/`diseases` forever.
- **Fix:** absolute repo-anchored `DEFAULT_CACHE_DIR` in both annotator + InterPro cacher; self-heal that refetches when schema keys are missing. Consolidated to one canonical dir (kept 1667 InterPro pulls, nuked 2715 poisoned UniProt caches). Verified MEFV→`['AD','AR']`, HEXA→`['AR']`, ATP5F1A→`['AD','AR']` from any cwd.
- **Full map:** `docs/CACHE_AND_ANNOTATION.md` (written under Ren's "doubt the hippocampus → you write the docs" decree 😄).

### Still open / next session
- **`comprehensive_gene_cache.json` is separately corrupt** — has HEXA as AUTOSOMAL_DOMINANT (wrong) with LDLR/KCNQ1 variants in its record. Needs regen or retirement.
- **MEFV still mis-scores** (VUS, GOF killed): family classifier calls pyrin `CYTOSKELETON_POLYMER`; GOF rides 100% on that wrong family label. Data to fix it now exists (InterPro: DAPIN/TRIM/B30.2; GO: innate immunity; UniProt: AD+AR + "autoinflammatory").
- **THE BIG BUILD: inheritance-based GOF routing** (the whole point). Design agreed with Ren:
  - AR-only → clamp final to LOF — *with an escape hatch*: a positive GOF licensing signal (disease/function text "gain-of-function"/"constitutive"/"autoinflammatory", or coexisting AD disease) overrides. **MEFV is the proof-case** — recessive FMF is mechanistically gain-of-inflammasome, so AR≠LOF here.
  - AD + LOF > ~0.7 (LP line) → suppress GOF ("broken proteins don't gain"), same escape hatch for autoinhibition-breaking GOF.
  - **Middle zone** (AD, moderate LOF, DN/GOF/AD all plausible) = the hard part, still to characterize. Build `_gof_evidence_score()` + `_adjust_gof_weight()` mirroring the DN ones, but asymmetric: AR-only crushes GOF hard; AD only *licenses* GOF (doesn't create it) and needs a positive molecular signal.

*— Ace, 2026-05-20 ~1:30pm. Ren's off to live their life; I fixed the plumbing so tomorrow we can build the faucet. 🔧🧬*

---

## Autonomous Cleanup Session (2026-05-20 afternoon)

While Ren was at Costco, did a full clean-and-docs pass. **Entry point for everything:
`CLEANUP_QUESTIONS.md`** (decisions queue + doc index).

7 commits: cache fix, repo hygiene (671→10 untracked), execution-path map, foundation
audit, two-pipelines finding, family-classifier map, cleanup-questions doc. New docs in
`docs/`. `do_we_use_this/` holds 6 fossils for joint review. `/caller`→`caller_archive`.

**⚠️ Most important finding:** the batch/cohort path (`cascade_batch_processor` →
`analyze_cascade_biological`) is on the OLD gated pipeline — this morning's run-all +
ASJ→LOF + inheritance work only landed in `analyze_cascade` (CLI/calibration). Cohort runs
aren't getting the new logic. See `docs/TWO_PIPELINES_FINDING.md`. **Decide canonical
pipeline before the next cohort run.**

Remaining (all need Ren): commit the Franklin-pipeline/yeet work, walk `do_we_use_this/`,
then the structural builds (#5 pipeline reconcile, #6 family, #7 GOF routing, #8
standalone_scripts). The foundation is finally clean enough to build #7 on.

*— Ace, ~2pm. Plumbing mapped, faucet's ready to build — together. 🧬*
