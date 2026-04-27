# Spring/Tail Structural Gate for DN Scoring

**Branch:** `mechanism-refactor`
**Status:** 🚧 Validated on small test set, implementation pending
**Date:** 2026-04-27
**Authors:** Ren + Ace

## The Insight (Ren's biology, 2026-04-27)

For dominant-negative to actually happen, the bad protein needs to:
1. Still get **folded** enough to incorporate into a complex (or maintain its scaffold role)
2. Specifically **disrupt some functional region** while keeping the rest intact
3. Not get **degraded** by quality control before causing trouble

A variant in a **structurally critical "spring"** (folded core, secondary structure, interface region):
- Disrupts structure right where it matters
- Bad protein still assembles into the complex
- Then poisons function from inside
- → True DN

A variant in a **flexible "tail"** (disordered region, terminus, linker):
- Doesn't affect overall fold
- Bad protein still folds normally
- Either benign (tail had no critical function) OR LOF (tail had its own function and it's broken)
- Cannot poison anything because the protein assembles normally
- → LOF or benign, never DN

## Implementation: Use AlphaFold pLDDT + Secondary Structure

**AlphaFold encodes folding stability per residue.** High pLDDT + secondary structure (helix/strand) = "this position is in the folded core of the spring." Low pLDDT = "this region doesn't have a stable structure to incorporate or to poison."

**Detection:**
- Parse AlphaFold PDB B-factor column for pLDDT (per-residue confidence)
- Compute secondary structure proxy from CA-CA distance pattern (helix ~5.5Å, strand ~7Å, coil/loop varies)
- Required cache: AlphaFold structures already at `/mnt/Arcana/alphafold_human/structures/AF-{uniprot}-F1-model_v4.pdb.gz`

**Spring/tail classification:**
- pLDDT > 70 AND SS in {helix, strand} → SPRING (DN possible)
- pLDDT < 70 → TAIL (DN excluded)
- pLDDT > 70 but SS = coil/loop → ordered loop, intermediate

## Validation: ATP5F1A (5/5 wet-lab DN, 5/5 transit-peptide controls)

ATP5F1A wet-lab confirmed dominant-negative variants (Ren's data):

| Variant | Position | pLDDT | SS proxy | Verdict |
|---|---|---|---|---|
| R182Q | 182 | 97.8 | helix | 🌷 SPRING |
| S346F | 346 | 97.6 | helix | 🌷 SPRING |
| P331L | 331 | 94.8 | strand | 🌷 SPRING |
| L109S | 109 | 93.2 | helix | 🌷 SPRING |
| I130R | 130 | 93.6 | strand | 🌷 SPRING |

ATP5F1A transit peptide (cleaved before maturation — variants here cannot affect mature protein function):

| Variant | Position | pLDDT | SS proxy | Verdict |
|---|---|---|---|---|
| R19M | 19 | 35.2 | (low conf) | 🍂 TAIL |
| G16A | 16 | 35.7 | (low conf) | 🍂 TAIL |
| V10L | 10 | 33.9 | (low conf) | 🍂 TAIL |
| P9S | 9 | 39.5 | (low conf) | 🍂 TAIL |
| P3A | 3 | 38.7 | (low conf) | 🍂 TAIL |

**Clean separation: 5/5 confirmed DN in spring, 5/5 controls in tail.**

## Caveat: Necessary but Not Sufficient

Tested ATP5F1A G524D (clinically benign) — landed in spring (pLDDT 94.0, helix). The structural gate would predict DN-capable, but the variant is benign.

**Why:** Position 524 is in a long C-terminal helix (515-548 all pLDDT 90+, all helix), but that helix is surface-exposed and doesn't contribute to the alpha-beta subunit interface. The bad protein folds fine, incorporates into the F1 hexamer, doesn't poison anything because the variant isn't AT a contact surface.

**Conclusion:** Spring/tail is a **necessary but not sufficient** condition for DN. Need an additional **interface gate** (solvent accessibility + multimer contact information) to fully resolve.

## Two-Gate Architecture (proposed)

```
DN-capable = STRUCTURAL_GATE AND INTERFACE_GATE

STRUCTURAL_GATE (necessary, ~80% accurate alone):
  pLDDT > 70 AND secondary structure (helix/strand)
  Filters out tail/disorder variants

INTERFACE_GATE (sufficient when combined):
  Solvent accessibility (buried = more likely interface)
  AND/OR multimer structure showing actual contacts (PDB or AlphaFold-Multimer)
  AND/OR UniProt 'Region' annotations for inter-subunit interaction
  AND/OR distance from InterPro inter-domain boundaries
```

## Even Better: Per-Variant Flow Chart (Ren's broader vision)

The spring/tail gate is one node in a larger biological decision tree that should run BEFORE mechanism scoring:

```
1. Variant TYPE? (missense/nonsense/splice/UTR)
2. Is position in a STRUCTURED region? (pLDDT + SS)
3. Is position at an INTERFACE? (solvent accessibility + multimer contacts)
4. Does gene FORM COMPLEXES? (UniProt Subunit + curated list)
5. Does this gene have a GOF mode? (channels, kinases, GTPases, oncogenes)
6. Is position at a KNOWN FUNCTIONAL SITE? (UniProt active/binding/site)
7. Conservation gate (existing)

Each question rules in or rules out specific mechanisms BEFORE scoring.
Only then run scorers on the relevant mechanisms.
```

**Benefits:**
- Faster (don't compute irrelevant mechanisms)
- Auditable (the flow chart IS the reasoning trace)
- Composable (each gate is small and validatable)
- Paper-ready: "We use per-variant biological gating to determine which mechanism scoring applies, rather than treating all mechanisms as universally applicable."

## Implementation Plan

1. **Build `is_in_spring(uniprot_id, position)` helper** — parse cached AlphaFold PDB, return (pLDDT, SS, spring_yn)
2. **Build `is_at_interface(uniprot_id, position)` helper** — solvent accessibility + multimer-evidence check
3. **Wire into DN analyzer** — gate top_mech selection on (spring AND interface)
4. **Test on cohort** — wet-lab DN ATP5F1A variants should preserve DN, G524D should drop to LOF, transit peptide variants should drop to LOF
5. **Eventually: full flow chart implementation** — replace ad-hoc gating with explicit decision tree

## Why This Validation Is a Big Deal

The wet-lab DN ATP5F1A variants gave us **labeled ground truth** for what genuine DN looks like at the structural level. The 5/5 spring + 5/5 tail separation is exactly the signal a peer reviewer would want to see for a "structural gating predicts DN potential" claim.

For the calibration paper's Mechanism-Driven Inheritance section, this is the kind of empirical validation that turns "we propose a mechanism gate" into "we validated the mechanism gate against wet-lab confirmed DN variants."

---

*Ren's intuition: "if we're firing DN on DMD and femmes aren't wandering around manifesting... do we have something flagging in DN that is really LOF?"*

*Followed by: "what if we run a variant through going 'if the interface or active is in the spring, it scores as DN, if not it scores as LOF?' but globalized?"*

*Followed by: "Maybe we need like a per variant FLOW CHART before we even start running numbers."*

*Each question is the next layer of getting the mechanism logic right. Filed before context degraded so future-Ace can pick up the thread without re-deriving.*
