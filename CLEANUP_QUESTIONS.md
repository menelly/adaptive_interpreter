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
- **#5 Reconcile the two DN brains.** v1 (`dn_mechanism_filter`) used by `BiologicalRouter`
  (routing); v2 (`dn_mechanism_filter_v2`) used by `NovaDNAnalyzer` (scoring). Do they
  agree? Collapse to one?
- **#6 Family classification.** Live family logic confirmed in `lattice_disruption/`
  (NOT dead — corrected) + smeared across ~10 files. Consolidate, or replace
  family-for-routing with inheritance+InterPro+GO (your instinct).
- **#7 GOF routing** (the original goal). AR-only→LOF with escape hatch (MEFV proof),
  AD+LOF>0.7→suppress GOF, middle zone TBD. Build `_gof_evidence_score()`.
- **#8 `standalone_scripts/`.** 84 of 122 .py are runnable — separate scripts from
  library so things are findable. Careful (path/cwd assumptions); do incrementally.

## Done autonomously (committed, no action needed)
Cache bug fix (`e4b9901`), repo hygiene 671→10 untracked (`92f5309`), execution-path map
+ fossil archive (`409edd7`), foundation-cache audit (`aff6b86`). Docs:
`docs/CACHE_AND_ANNOTATION.md`, `docs/EXECUTION_PATH_MAP.md`. `/caller`→`caller_archive`.
