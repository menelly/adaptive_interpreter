# Family Classification — where it happens, why it misfires (Task #6 map)

> Mapped 2026-05-20 by Ace (analysis only — the fix is a with-Ren scoring decision).

## The authoritative live classifier

`utils/plausibility_filter.py:139` — **`classify_gene_family(gene_symbol, uniprot_function,
go_terms, override_family)`**, called at `:511` inside `apply_plausibility_filter`, printed
at runtime as `🧬 Family:` (`cascade_analyzer.py:1609`).

**Method:** keyword / GO-term text matching against UniProt function text ("Nova's weighted
scoring"). It assigns the family label that drives the GOF/DN/LOF multipliers in
`FAMILY_PLAUSIBILITY`.

**The tell that it's unreliable:** it carries a hardcoded `CURATED_OVERRIDES` dict (lines
167+) — a hand-maintained patch table for genes the keyword classifier gets wrong:
ATP7A/B (caught METABOLIC via catalytic keyword → forced TRANSPORTER), MYH9 (caught
CYTOSKELETON_POLYMER via 'actin' → forced MOTOR_PROTEIN), HEXA/GAA (→ METABOLIC_ENZYME),
cadherins (→ SCAFFOLD_ADAPTOR), etc. The override table grows every time a gene misfires.
**MEFV/pyrin is not in the table → it gets the wrong keyword answer (`CYTOSKELETON_POLYMER`)**,
which hands GOF a 0.1 multiplier and kills the (correct) gain-of-inflammasome mechanism.

## Competing family notions (the "smeared across files" problem)

These are *independent* classifiers with *different* taxonomies — no single source of truth:

| function | file | taxonomy / method | role |
|---|---|---|---|
| `classify_gene_family` | `utils/plausibility_filter.py:139` | keyword/GO + override table (CYTOSKELETON_POLYMER, METABOLIC_ENZYME, ...) | **live** — drives plausibility multipliers |
| `detect_protein_family` | `nova_dn/lattice_disruption/family_detector.py:11` | sequence-based (COLLAGEN/COILED_COIL/ION_CHANNEL) | live — structural routing in `universal_router` |
| `get_gene_family` | `utils/ensemble_scores.py:116` | (separate) | check usage |
| `get_gene_family` | `utils/family_aware_gly_cys_trainer.py:82` | training-time | offline ML |

`mixed_mechanism_resolver.py` and the ML integrators *consume* `context["gene_family"]` but
don't classify — they trust whatever was set upstream.

## Structural options for the fix (Ren decides)

1. **Ren's instinct (recommended):** stop using the family *label* to route DN/LOF. Route
   those by **UniProt inheritance + InterPro domains + GO** (all now reliably available post
   cache-fix). Keep `classify_gene_family` only for **GOF suppression** (where family
   genuinely helps: METABOLIC_ENZYME/STRUCTURAL → GOF implausible). Shrinks the blast radius
   of the unreliable classifier.
2. **Replace the keyword classifier with InterPro-domain-based classification** — we now pull
   real domains (DAPIN/PYD, B30.2/SPRY, TRIM, catalytic, etc.), which are far more reliable
   than function-text keywords. Bigger change.
3. **Bandaid:** add MEFV (+ other inflammasome/TRIM genes) to `CURATED_OVERRIDES`. Fast, but
   it's another row on the Jenga tower — the thing we're trying to stop doing.

Connects to **Task #7 (GOF routing):** the inheritance-based GOF gate is option 1's natural
partner — AR-only suppresses GOF regardless of family, with a positive-licensing escape
hatch (MEFV: gain-of-inflammasome under recessive inheritance).
