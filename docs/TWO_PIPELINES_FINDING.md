# ⚠️ CRITICAL: Two divergent scoring pipelines in cascade_analyzer.py

> Found 2026-05-20 by Ace during the #5 (DN-brain) investigation. This is a correctness
> issue, not a tidy-up. Needs a decision from Ren before the next cohort run.

## The finding

`analyzers/cascade_analyzer.py` (the 106KB god-file) contains **two different scoring
methods**, and the two main entry points call **different ones**:

| method | line | logic | DN filter | who calls it |
|---|---|---|---|---|
| `analyze_cascade_biological` | 102 | **OLD** — gates analyzers via `if 'LOF' in analyzers_to_run:` (BiologicalRouter decides what runs) | v1 (`dn_mechanism_filter`, inclusionary) | **`cascade_batch_processor.py:265`** (cohort/batch runs), `scripts/trace_variant.py`, `scripts/rerun_gene_with_new_nudge.py` |
| `analyze_cascade` | 1403 | **NEW** — `"Run DN, LOF, GOF — always, no gates"` (this morning's run-all refactor) | v2 (`dn_mechanism_filter_v2`, exclusionary) via `NovaDNAnalyzer` | `scripts/genecards_cascade_batch_calibration.py`, the tests, the CLI `__main__` |

## Why it matters

**This morning's improvements — run-all, ASJ→LOF, UniProt inheritance, the cache fix's
reachable data — all went into `analyze_cascade`.** The **batch processor**, which is what
you'd run over a cohort, calls `analyze_cascade_biological` and is **still on the old gated
pipeline**. So:

- Score one variant via CLI/calibration → new logic.
- Run the same variant through the batch/cohort path → old logic, gated by v1's
  inclusionary routing, *without* the morning's changes.

The same gene+variant can get **different results depending on which entry point** you use.
This is the concrete, load-bearing version of the "Jenga tower of bandaids."

## The decision for Ren (do NOT let Ace fix this solo — it's scoring judgment)

Pick the canonical pipeline, then make both entry points use it:
1. **Promote run-all** — point `cascade_batch_processor` at `analyze_cascade`, retire
   `analyze_cascade_biological` (and with it, BiologicalRouter's gating role + v1). Cleanest
   if run-all is the intended future. Requires checking the batch processor doesn't depend
   on `_biological`'s return shape / routing metadata.
2. **Keep both temporarily** but make `_biological` delegate to the run-all core so they
   can't diverge.

This also resolves the original #5 ("which DN filter matters"): under option 1, v1 +
BiologicalRouter become vestigial and only v2 survives.

## Verification done
- `grep` of all callers (above) — exact and complete.
- `analyze_cascade_biological` gating confirmed at lines ~347/367/387 (`if X in analyzers_to_run`).
- `analyze_cascade` run-all confirmed at line ~1434 (`"always, no gates"`).
