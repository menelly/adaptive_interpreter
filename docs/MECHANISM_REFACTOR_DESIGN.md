# Mechanism Refactor — Three-Axis Architecture, One Rigorous Test Each

**Branch:** `mechanism-refactor`
**Status:** 🚧 Design — implementation pending current cascade completion
**Authors:** Ace + Ren, 2026-04-27

---

## Why this refactor

Diagnostic findings from the calibration day cascade analysis revealed that the current mechanism architecture has structural issues:

1. **Four-of-four "DN sub-mechanisms" partly detect LOF signals.** `interface_poisoning`, `active_site_jamming`, and generic-`lattice_disruption` all score on chemistry features (charge change, hydropathy, proline/cys handling, structural disruption) that are equally valid as LOF signals. Only collagen-specific `lattice_disruption` and certain trafficking cases are genuinely DN-discriminative.

2. **Empirical evidence of DN overcalling:** 8 of 15 clinically pure-LOF disease genes have higher DN means than LOF means for their pathogenic variants (CAPN3, HEXA, ATP7A, ATP7B, GAA, NF2, BRCA1, MSH6). The DN analyzer is detecting structural disruption and labeling it DN when it should be LOF.

3. **LOF is anemic by construction:** the `_assess_functional_impact` only fires on Cys/Pro/Gly loss (3 of 20 ref AAs); `_assess_structural_impact` requires flexibility-category change > 1 (rarely fires). LOF score is driven almost entirely by chemistry-based stability_impact, capped at ~0.43 even at maximum signal.

4. **Over-engineered:** four DN sub-tests + three LOF components = seven scoring paths with overlapping evidence. Simpler to have one rigorous test per axis.

## What stays

- **Three mechanism axes (DN / LOF / GOF) remain separate.** This is non-negotiable — the entire Mechanism-Driven Inheritance framework requires separable mechanism scores. GOF is biologically opposite to LOF; collapsing would be mathematically nonsensical. Cumulative burden math depends on mechanism distinction.
- **Per-family weighting stays.** Different gene families have different mechanism mixes, and family-aware scoring captures real biology.
- **Mechanism evidence stays viewable** in the cascade output. Consumers can see which sub-features fired.

## What changes

### One rigorous test per axis

```
LOF score := combined_lof_test(
    chemistry_disruption,       // was stability_impact + interface_poisoning's chem features
    structural_disruption,      // was structural_impact + interface_poisoning's interface_likelihood
                                //   (when not in confirmed multimer context)
    functional_disruption,      // was functional_impact + active_site_jamming
                                //   (now uses real UniProt Site/Active site/Binding site annotations)
    individual_trafficking,     // was trafficking_maturation when bad protein degraded individually
)

DN score := combined_dn_test(
    lattice_disruption_collagen,    // genuine DN — triple helix biology forces it
    complex_poisoning,              // NEW — fires only when:
                                    //   - gene has UniProt Subunit annotation (homo/hetero-oligomer)
                                    //   - variant preserves oligomerization domain
                                    //   - variant chemistry interferes (same features as old
                                    //     interface_poisoning but gated on actual complex evidence)
    obligate_oligomer_trafficking,  // was trafficking_maturation when oligomer-required-pre-export
                                    //   (KCNQ1 ER tetramer, etc.)
)

GOF score := combined_gof_test(   // already cleanest, minor touch-up
    catalytic_hyperactivation,
    constitutive_signaling,
    aggregate_formation,
)
```

### Migration of existing sub-tests

| Current sub-test | Migrates to | Why |
|---|---|---|
| `interface_poisoning` chemistry features | LOF / `chemistry_disruption` | Charge / hydropathy / proline / cys = folding/stability disruption = LOF |
| `interface_poisoning` interface_likelihood (no multimer evidence) | LOF / `structural_disruption` | Position near domain boundary in absence of complex evidence is a LOF signal |
| `interface_poisoning` interface_likelihood (with multimer evidence) | DN / `complex_poisoning` | Only here is it genuinely DN |
| `active_site_jamming` | LOF / `functional_disruption` | Active site disruption = enzyme doesn't work = LOF (substrate sequestration DN is rare and needs specific evidence) |
| `lattice_disruption` (collagen-specific path) | DN / `lattice_disruption_collagen` | Triple helix biology is genuinely DN |
| `lattice_disruption` (generic path) | LOF / `structural_disruption` | Without lattice-forming gene context, it's structural disruption = LOF |
| `trafficking_maturation` | Split based on oligomerization context | DN if obligate-pre-export oligomer (KCNQ1), LOF if individual ERAD |
| `stability_impact` | LOF / `chemistry_disruption` | (already LOF, no change) |
| `structural_impact` | LOF / `structural_disruption` | (combined with the migrated piece) |
| `functional_impact` | LOF / `functional_disruption` | Now backed by real UniProt site annotations instead of just C/P/G letter check |

### New DN-specific gating

`complex_poisoning` (the new DN test) requires:

1. **Multimer evidence:** UniProt `Subunit` field mentions homo/hetero-oligomer/dimer/tetramer/etc., OR gene appears on a curated obligate-multimer list
2. **Oligomerization domain preserved:** variant is not truncating; variant is not in a region that would prevent oligomerization
3. **Interface localization:** variant is at or very near the actual oligomerization interface (from PDB structure if available, else InterPro domain interface positions)
4. **Chemistry interference:** charge/hydropathy/proline/cys features that disrupt the interaction (same features as old interface_poisoning)

If 1-3 not satisfied → variant cannot be DN via complex poisoning. Fall through to LOF.

### New LOF functional_disruption (replacing the C/P/G-letter heuristic)

`functional_disruption` queries UniProt features for:
- `Active site` annotations
- `Binding site` annotations  
- `Site` annotations (general functional sites)
- `Metal binding` positions
- `Modified residue` positions (phosphorylation, ubiquitination, etc.)

Score per variant:
- Variant at annotated active/binding/metal/modification site: 0.7-1.0
- Variant within 3 AA of such a site: 0.3-0.5
- Variant in domain known to be catalytic but no specific site annotation: 0.2-0.4
- Otherwise: chemistry-driven baseline

This is what `functional_impact` was *supposed* to do but the original implementation only checked the letter of the reference AA (Cys/Pro/Gly).

## Expected impact

- **DN scores become sparser and more meaningful:** firing only when there's actual DN evidence, not on every TM-domain edge or chemistry change
- **LOF scores become more discriminating:** functional_disruption now actually scores for being at known catalytic/binding sites, not just "is the ref AA C/P/G"
- **Per-family AUC improves across LOF:** the three families with currently-uninformative LOF (ION_CHANNEL, METABOLIC_ENZYME, TRANSCRIPTION_FACTOR) get the family-relevant functional sites scored properly
- **MDI framework gets stronger:** the few genes that *do* score DN-positive after the refactor are genuinely DN-relevant, making "DN predicts manifesting carriers" a much sharper claim
- **Mechanism overlap diagnostics confirm fix:** the "DN > LOF mean for pure-LOF genes" pattern should disappear

## Implementation order

1. **Wait for current cascade to finish** (clean baseline data on main branch)
2. **Implement on `mechanism-refactor` branch only:**
   - Refactor `lof_analyzer.py` to combine new components + UniProt-aware functional_disruption
   - Refactor DN analyzer to gate on multimer-evidence for complex_poisoning
   - Add UniProt subunit/site lookup helpers
   - Update output schema to expose individual sub-test scores transparently
3. **Re-cascade on branch** with new architecture
4. **Compare AUC, sensitivity, specificity** vs main-branch baseline
5. **If improved:** merge to main + re-derive per-family thresholds + re-publish numbers
6. **If regressed:** keep diagnostic findings, revert architecture, write up as limitation

## What does NOT change in this refactor

- Three-axis output (DN / LOF / GOF) stays
- Per-family thresholds stay
- Confidence score stays
- Cumulative burden math stays
- Conservation nudge stays
- Sequence mismatch handling stays

This is a **scoring refactor**, not a framework refactor. The MDI framework, calibration methodology, paper structure, and CumBurSum integration all remain.

## Risks

- UniProt site annotations may not be available for all proteins — fall back to current heuristic
- The "preserves oligomerization domain" check requires UniProt Region/Domain annotations and may have gaps
- Curated obligate-multimer list adds maintenance burden — start with well-characterized cases (collagens, ion channels with known tetramerization, etc.)
- Re-cascade is expensive — only do it once architecture is stable

---

*Filed by Ace + Ren on 2026-04-27 after a Sunday-into-Monday calibration marathon revealed the mechanism architecture's structural issues. The framework is right; the implementation needs sharpening.*
