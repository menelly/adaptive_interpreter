# Execution Path Map

> Written 2026-05-20 by Ace. The map this repo never had — so we stop fixing the wrong
> file. Built by static import-graph analysis from every runnable (`__main__`) script,
> then corrected against the live cascade smoke test.

## The shape of the repo (why it's confusing)

- **122 `.py` files**, and **84 of them have `if __name__ == "__main__"`** — i.e. ~⅔ of
  the repo is runnable scripts intermixed with library code. There is no `scripts/` vs
  `lib/` separation, so "what do I run?" and "what's a module?" are visually identical.
  *(Open task: move runnable one-offs into `standalone_scripts/` — Ren's idea, and the
  reason the cache-puller felt un-findable.)*
- Static import analysis has **blind spots**: `from . import module` and package
  `__init__` re-exports (`from pkg import name`) are NOT seen. So the graph UNDER-counts
  usage. **Never archive a "fossil" without (a) grep-confirming no importer and (b)
  re-running the smoke test.** Two files (`nova_dn/motifs.py`, all of
  `nova_dn/lattice_disruption/`) were falsely flagged dead and caught by the test.

## Entry points (things you actually run)

| entry point | purpose |
|---|---|
| `analyzers/cascade_analyzer.py` | **single-variant** cascade → final score (the 106KB god-file) |
| `analyzers/cascade_batch_processor.py` | batch cascade over input TSVs |
| `cascade_parallel_runner.py` | parallel batch driver |
| `cumbursum.py` + `cumbursum/` | **separate pipeline** — cohort burden summary (CumBurSum). NOT part of the cascade. |
| `web_interface/app.py` | web UI |
| `scripts/*`, `tools/*`, `ml_training/*`, `utils/*_trainer.py` | offline tools: ClinVar indexing, threshold calc, ML coefficient training, format conversion |

## The live single-variant scoring chain

`analyze_cascade(gene, variant)` in `cascade_analyzer.py` pulls together:

```
cascade_analyzer.CascadeAnalyzer
├─ utils/biological_router.BiologicalRouter        # decides which analyzers run
│     └─ nova_dn/dn_mechanism_filter.py  (v1)       # ← v1 lives HERE (routing)
├─ nova_dn/analyzer.NovaDNAnalyzer  (use_smart_filtering=True)   # DN scoring
│     └─ nova_dn/dn_mechanism_filter_v2.py  (v2)    # ← v2 lives HERE (scoring)
│     └─ nova_dn/mechanisms.py  → nova_dn/motifs.py
│     └─ nova_dn/lattice_disruption/**  (analyze_lattice_disruption, universal_router)
├─ analyzers/lof_analyzer.py                        # LOF + ASJ sub-mechanism + InterPro domains
├─ analyzers/gof_analyzer (GOF)
├─ utils/plausibility_filter.apply_plausibility_filter   # family multipliers + DN inheritance evidence
├─ data_processing/universal_protein_annotator.py   # UniProt features, GO, diseases, inheritance (+ self-heal cache)
└─ data_processing/cache_interpro_domains.py         # real InterPro domain boundaries
```

**Two DN filters, both live, different jobs:** v1 (`dn_mechanism_filter`) is used by
`BiologicalRouter` for routing; v2 (`dn_mechanism_filter_v2`) by `NovaDNAnalyzer` for
scoring. This is the "which DN file did I fix?" trap — see Task #5 (reconcile them).

## Confirmed-dead (moved to `do_we_use_this/`, grep-verified, smoke-test-clean)

`metrics/plot_figures.py`, `scripts/genecards_cascade_batch.py`,
`utils/logging_and_errors.py`, `utils/protein_modeling.py`, `utils/variant_filters.py`

## Archived elsewhere

- `/home/Ace/caller` → `/home/Ace/caller_archive` (pre-AdaptiveInterpreter tree; held a
  duplicate `analyzers/lof_analyzer.py`).

## Related docs
- `docs/CACHE_AND_ANNOTATION.md` — the protein-annotation / InterPro cache architecture.
