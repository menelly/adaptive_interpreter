# Cleanup — decisions that need Ren

Ace did an autonomous cleanup pass 2026-05-20. Everything below is either a decision
only Ren can make, or a structural task we agreed to do *together*. Walk this when back.

## A. `do_we_use_this/` — delete-for-real or restore? (joint review)
See `do_we_use_this/MANIFEST.md`. 6 items, all grep-confirmed unused + smoke-test-clean:
- `metrics/plot_figures.py`, `scripts/genecards_cascade_batch.py`,
  `utils/{logging_and_errors,protein_modeling,variant_filters}.py`
- `data_processing/clinvar_inheritance_cache.py` + corrupt `comprehensive_gene_cache.json`
  (the reader was orphaned; corruption never touched live scoring)

## B. Commit your in-progress work? (your call — I won't commit your code behind your back)
Untracked, all coherent (not scratch):
- **Franklin cohort pipeline:** `run_cohort.py`, `cumbursum/franklin_adapter.py`,
  `cumbursum/vcf_to_discovery.py`
- `cumbursum/panels_loader.py` (already imported by `cumbursum.py` — live), `cumbursum/panels/*.gmt`, `cumbursum/STRING_BRIDGE_NOTES.md`
- `yeet_synonymous.py` (the missense prefilter)
- `web_interface/satire/index.html`
- **119 modified `inputs_missense_only/*.tsv`** = the result of running `yeet_synonymous.py`
  (intronic/synonymous/benign rows stripped). Backup preserved on disk (gitignored
  `inputs_missense_only_PREFILTER_BACKUP/`). Commit the filtered corpus? (Want to eyeball
  the filter did the right thing first?)

## C. Structural tasks — agreed to do TOGETHER (need your brain / consent)
- **#5 ⚠️ CRITICAL — two divergent scoring pipelines (bigger than just the DN filters).**
  See `docs/TWO_PIPELINES_FINDING.md`. The **batch/cohort path** (`cascade_batch_processor`
  → `analyze_cascade_biological`) is on the OLD gated logic with v1; the morning's run-all
  + ASJ→LOF + inheritance work only landed in `analyze_cascade` (CLI/calibration/tests).
  **Cohort runs are currently NOT getting the new logic.** Decide canonical pipeline + make
  both entry points use it (this also resolves v1-vs-v2 — v1 likely becomes vestigial).
  Do NOT let Ace fix solo — scoring judgment.
- **#6 Family classification.** Live family logic confirmed in `lattice_disruption/`
  (NOT dead — corrected) + smeared across ~10 files. Consolidate, or replace
  family-for-routing with inheritance+InterPro+GO (your instinct).
- **#7 GOF routing** (the original goal). AR-only→LOF with escape hatch (MEFV proof),
  AD+LOF>0.7→suppress GOF, middle zone TBD. Build `_gof_evidence_score()`.
- **#8 `standalone_scripts/`.** 84 of 122 .py are runnable — separate scripts from
  library so things are findable. Careful (path/cwd assumptions); do incrementally.

## Done autonomously (committed, no action needed)
- Cache bug fix (`e4b9901`) — unified path + self-heal; inheritance reachable
- Repo hygiene 671→10 untracked (`92f5309`) — `.pyc`/outputs sorted, `/caller`→`caller_archive`
- Execution-path map + 6 fossils → `do_we_use_this/` (`409edd7`)
- Foundation-cache audit (`aff6b86`) — corrupt cache was orphaned, not live poison
- **#5 critical finding** (`6da7982`) — batch vs single use different pipelines
- **#6 family-classifier map** (`d30be39`)

**Docs to read (all in `docs/`):** `TWO_PIPELINES_FINDING.md` ⚠️ (read first),
`FAMILY_CLASSIFICATION_MAP.md`, `EXECUTION_PATH_MAP.md`, `CACHE_AND_ANNOTATION.md`.
This file (`CLEANUP_QUESTIONS.md`) is the entry point.
