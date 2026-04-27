# CumBurSum Calibration Methodology — The Reference Standard Move
## Ren's Aquatica Insight, April 26, 2026

**Origin:** Ren at a water park with her kids, texting chat-Ace, pasting to Linux-Ace. Distributed cognition in action.

---

## The Problem

CumBurSum produces throughput percentages (e.g., Ren's mito stack = 7%). These numbers are internally consistent but **uncalibrated** — we don't know what 7% *means* clinically. Is it catastrophic? Surprisingly functional? Where's the disease threshold?

## The Solution: Known-Pathogenic Calibration

Run the burden calculator on **known single-variant pathogenic cases** where the clinical answer is already settled. Each produces a throughput score. That score IS, by definition, the pathogenicity threshold for that pathway.

### High-End Calibration (Known Pathogenic)

| Variant | Gene | Condition | Inheritance | Clinical Severity |
|---|---|---|---|---|
| R22W | TFG | HMSN-P | AD | Severe neuromuscular |
| Glycine substitutions | COL1A1 | OI Type II | AD | Lethal/severe skeletal |
| R94Q | MFN2 | CMT2A | AD | Severe neuropathy |
| Various | SCN1A | Dravet syndrome | AD | Severe epilepsy |
| ΔF508 homozygous | CFTR | CF | AR | Severe pulmonary |

### Low-End Calibration (Known Asymptomatic Carriers)

| Variant | Gene | Condition | Carrier Status | Clinical Severity |
|---|---|---|---|---|
| ΔF508 heterozygous | CFTR | CF carrier | Het | Asymptomatic |
| Various deep-AR hets | Various | Carrier only | Het | Asymptomatic |

### What This Gives Us

1. **A calibration curve** — continuous score between "definitely pathogenic" and "definitely asymptomatic"
2. **Clinical reference thresholds** — "this throughput level is known to cause disease in [comparison disorders]"
3. **AD-vs-AR dosage signal** — empirical measurement of how much pathway throughput each mechanism preserves
4. **Semi-Dominant Hypothesis quantification** — AdaptiveInterpreter showed DN > LOF predicts semi-dominant qualitatively; calibrated burden calculator gives the NUMBERS

### The Ren Question

When calibrated, Ren's 7% mito throughput can be stated as:

> "This patient's pathway throughput is below the asymptomatic-carrier calibration floor by X standard deviations and within the known-pathogenic range observed for [comparison disorders]."

That's not "we calculated some math." That's **clinical evidence with reference-standard backing.**

### Three Possible Outcomes for the 7%

1. **Math needs tweaking** — Ren acknowledged this possibility. Calibration will expose systematic bias.
2. **Compensatory mechanisms** — Ren has uncharacterized compensation keeping her functional despite catastrophic throughput. Publishable.
3. **Outperforming predicted clinical course** — "The patient who refused to die when the math says she should have." Also publishable.

All three are findings. The calibration tells us which one.

## Implementation Plan

1. Curate 50-100 known-pathogenic single variants (ClinVar P/LP with strong evidence)
2. Curate 20-30 known-asymptomatic carrier states
3. Run each through CumBurSum as if it's the ONLY variant in the pathway
4. Plot throughput score vs. clinical severity
5. Establish threshold bands: catastrophic / severe / moderate / carrier / benign
6. Re-run Ren's full mito stack against the calibrated scale
7. Write it up

## Dependencies

- CumBurSum core engine (DONE — `/home/Ace/CumBurSum/`)
- Pathway definitions with the relevant genes (PARTIALLY DONE — needs expansion)
- ClinVar data access for known-pathogenic variants (need to pull)
- AdaptiveInterpreter mechanism calls for each variant (DONE for missense)

---

*Captured by Linux-Ace from chat-Ace relay via Ren-as-corpus-callosum.*
*The synapse was at a water park. The motor cortex is on a Linux server. Science doesn't stop for wave pools.*
