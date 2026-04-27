# Calibration Documents — Status

🚧 **This entire folder is WORK IN PROGRESS.** 🚧

These documents are the working scratchpad for the calibration extension of
AdaptiveInterpreter. The numbers and findings here are preliminary —
collected during an active cascade re-run and a series of bug fixes.
**Do not cite this folder; cite the published paper once final numbers land.**

## File status

| File | Status | What it is | When to retire/finalize |
|---|---|---|---|
| `HANDOFF_2026-04-26_CALIBRATION_DAY.md` | 📜 archival | Snapshot of one debugging session (3 bugs found + fixed; 11 dangerous flips revisited). Historical record. | Keep as-is for provenance. |
| `CALIBRATION_METHODOLOGY.md` | 📜 archival | Ren's water-park insight: use known-pathogenic single variants as the calibration reference scale. The IDEA is final; the implementation is the rest of these docs. | Keep as-is. |
| `CALIBRATION_RESULTS.md` | ⚠️ **PRELIMINARY** | First-pass numbers from partial cascade (61/127 genes done at time of writing). Subsequent re-cascade with bug fixes will supersede. | Replace when full cascade rerun completes. |
| `CALIBRATION_PAPER_v2_SCAFFOLD.md` | 🚧 **ACTIVELY DRAFTING** | Graft-points map showing where today's findings land in the existing FINAL paper. Working scaffold. | Retire once content is grafted into `AdaptiveInterpreter_Paper_Genetics_v2_FINAL.md`. |

## Known still-open items

- Full cascade rerun (~127 genes, ~150K variants) is in progress; final per-family thresholds depend on its completion
- Confidence-score implementation just landed; needs validation on the new cascade output
- 36 genes that emptied during prefilter (no v2 cascade to backfill, only `hgvs_g` inputs) need fresh GeneCards pulls
- Real `gnomad_freq` pipeline not yet wired (currently 0.0 placeholders)
- DNA_REPAIR family classifier needs JSON keywords added (BRCA1/2, MSH2/6, MLH1, MUTYH, RAD50 currently fall to GENERAL fallback)
- Frequency-flag (not nudge) per the design discussion — pending implementation
- Founder-mutation exception list partially implemented in confidence score — needs broader scope for cumulative-burden context

## Pointers

- The published paper draft lives at `/home/Ace/AdaptiveInterpreter_Paper_Genetics_v2_FINAL.md` (in main folder, not in this repo). Once the WIP content here is grafted in, this folder gets cleaned out.
- The confidence-score design is at `../../CONFIDENCE_SCORE_DESIGN.md` (one level up).

---
*Last updated: 2026-04-27. Maintained by Ace + Ren.*
