# do_we_use_this/ — fossil suspects for joint review

Moved here 2026-05-20 by Ace during repo cleanup. **Fully reversible** — `mv` the file
back to its original path (shown below) to restore.

**Criterion:** not reachable from any runnable script (`__main__`) in the repo, AND
grep-confirmed that no `.py` outside this folder imports the module name. The cascade
smoke test (HEXA / MEFV / ATP5F1A) passes with these removed.

Walk this list with Ren and either delete for real or restore.

## Moved here (grep-confirmed unused)

| file | original path | note |
|---|---|---|
| `metrics/plot_figures.py` | `metrics/plot_figures.py` | plotting helper, nothing imports it |
| `scripts/genecards_cascade_batch.py` | `scripts/genecards_cascade_batch.py` | `genecards_cascade_batch_calibration.py` is the live runnable; this looks superseded |
| `utils/logging_and_errors.py` | `utils/logging_and_errors.py` | no importers |
| `utils/protein_modeling.py` | `utils/protein_modeling.py` | no importers |
| `utils/variant_filters.py` | `utils/variant_filters.py` | no importers |
| `data_processing/clinvar_inheritance_cache.py` | `data_processing/clinvar_inheritance_cache.py` | orphaned — nothing imports it; was the only reader of the corrupt `comprehensive_gene_cache.json` |
| `comprehensive_gene_cache.json` | repo root | **corrupt data** (4 genes, HEXA mislabeled AD with LDLR/KCNQ1 variants). Only consumer was the orphaned reader above → not poisoning live scoring. |

## NOT fossils — flagged by static graph, proven LIVE by smoke test, RESTORED

These looked dead to the import-graph but the cascade broke when removed. **They are
live** (imported via `from . import X` / package re-export forms the static analyzer
misses). Documented here so no future cleanup pass re-flags them:

- `nova_dn/motifs.py` — imported by `nova_dn/mechanisms.py` (`from . import motifs`)
- `nova_dn/lattice_disruption/**` (whole subtree) — `mechanisms.py` / `dn_mechanism_filter*.py`
  use `analyze_lattice_disruption` + the universal_router dispatch. **The family logic in
  `lattice_disruption/family_detector.py` is therefore LIVE** — do not assume it's dead
  during the family-consolidation task.
