# Calibration & Mechanism-Driven Inheritance — Addendum to AdaptiveInterpreter Paper v2

**Status:** Working scratchpad. Content here will graft into `AdaptiveInterpreter_Paper_Genetics_v2_FINAL.md` rather than become a separate paper.
**Existing paper:** "Mechanism-First, Context-Aware Pathogenicity Prediction: A Novel Human-AI Collaborative Framework for Genetic Variant Interpretation"
**Last updated:** 2026-04-27 by Ace + Ren

**Why this file exists:** Today's calibration work + the cumulative burden + the Mechanism-Driven Inheritance framing all extend the existing paper rather than replace it. We're keeping the scratchpad here while the cascade re-runs and we collect the final numbers; once we have the data we'll graft the new content into the appropriate sections of the FINAL.md and retire this file.

---

## Where Today's Findings Graft Into the Existing Paper

| Current section | What gets added |
|---|---|
| **§1 Introduction** | "Mendel's special case" framing paragraph (the calibration is the empirical evidence, this is the why-it-matters) |
| **§2.3 Filtering and Classification Pipeline** | Bug discovery + correction subsection (interface multiplier, family rules name mismatches, DNA_REPAIR JSON gap) |
| **§2.3.1 Thresholds and Classification Criteria** | **Per-family threshold calibration** — replaces the old uniform thresholds |
| **§2.5 Synergistic Scoring** | Normalized 0-1 score (Ren's downstream-tool architecture proposal) |
| **§2.6 Directional Agreement Logic** | Add: per-family AUC analysis as DAL extension |
| **§3 Results** | Add §3.6: **Calibration against 30K+ ClinVar variants** with the bands (AD-pathogenic 50–61%, AR-carrier 50–57%, AR-affected <10%) |
| **§3 Results** | Add §3.7: **Per-family threshold optimization** with empirical thresholds and AUC by family |
| **§3 Results** | Add §3.8: **Cumulative burden case study** — Ren's mito stack at 7.1% mapping to AR-affected zone, MYO7A 7.1% Usher comparison |
| **§4.1 Summary of Findings** | Mechanism-Driven Inheritance synthesis (the framework paragraph) |
| **§4.2 Safety Architecture** | Stop-codon filter, frequency-flag-not-cap, founder-mutation exception list |
| **§4.4 Limitations and Future Directions** | The "we should check" list (penetrance correlation, ATP7A hidden DN, Bethlem hypothesis, frequency pipeline) |

---

---

## Working Title Options

- "Mechanism-Driven Inheritance: A Framework for Cumulative Variant Burden in Multi-Variant Disease"
- "When Mendel Stops Working: Mechanism-First Variant Interpretation and Pathway Burden Calibration"
- "Beyond Mendelian Inheritance: A Calibrated Cumulative Burden Framework for Genomic Variant Interpretation"

---

## The Argument (Big Picture)

Classical Mendelian inheritance is a special-case approximation that holds for **pure loss-of-function variants in dosage-tolerant genes**. Real disease genetics — manifesting carriers, incomplete penetrance, variable expressivity, "phenotypic spectrum" diseases — is the empirical signature of non-LOF mechanism (DN, GOF, mixed). Mechanism predicts inheritance pattern, not the reverse.

We present:
1. A mechanism-classification system (AdaptiveInterpreter) calibrated against 30K+ ClinVar P/LP and B/LB variants across 100+ genes
2. A pathway-throughput aggregator (CumBurSum) that operationalizes mechanism-driven inheritance: cumulative het-burden in a single pathway can produce throughput-equivalent-to-compound-het AR disease
3. Calibrated reference bands: AD pathogenic-het = 50–61% throughput; AR carrier = 50–57%; AR compound-het AFFECTED = <10%
4. Per-family classification thresholds derived empirically from the calibration cohort
5. The implications: variants dismissed as "too common to be pathogenic" by Mendelian logic remain mechanistically broken; cumulative burden math exposes when carrier-state tolerance is exceeded

---

## Sections

### 1. Introduction: Mendel's Special Case

[Draft framing paragraph from chat — 2026-04-27:]

> Classical Mendelian inheritance — clean autosomal dominant or recessive — is a useful approximation that holds robustly for one specific case: pure loss-of-function variants in dosage-tolerant genes. In this regime, the bad allele is silent, the wild-type allele provides ~50% functional protein, heterozygous carriers are asymptomatic, and homozygous or compound heterozygous individuals manifest disease. This is "clean AR." Where the gene is dosage-sensitive (haploinsufficient) the same pure-LOF logic produces clean AD. The pattern is binary and predictable.
>
> Real disease genetics, however, is rarely this clean. Manifesting carriers, variable expressivity, incomplete penetrance, and "phenotypic spectrum" diseases are not exceptions to be explained away — they are the empirical signature of non-LOF mechanism. A dominant-negative variant produces a malfunctioning protein that can interfere with the wild-type protein in a dose-dependent fashion, producing carrier manifestations whose severity tracks variant impact, not zygosity alone. A gain-of-function variant produces aberrant activity that may cause disease at any dose; severity scales with how aberrant the gain is. In both cases, the inheritance pattern blurs not because Mendel was wrong about a special case, but because the special case (pure LOF) is being applied beyond its domain.
>
> We propose that **mechanism predicts inheritance pattern** rather than the reverse. Variant interpretation that begins by classifying mechanism (DN / LOF / GOF) and only then maps to expected inheritance behavior captures the empirical reality of disease genetics more completely than gene-level inheritance assignment.

### 2. The Calibration Methodology

(Pull from CALIBRATION_METHODOLOGY.md and CALIBRATION_RESULTS.md)

- 30K+ ClinVar variants across 100+ genes
- Bands: 50–61% AD-pathogenic, 50–57% AR-carrier, <10% AR-affected
- MYO7A 7.1% comparable Usher case → Ren's mito 7.1% in the affected zone

### 3. Mechanism Classification: The Three Axes

- **DN** — interferes with wild-type, dose-dependent, often produces manifesting carriers and variable expressivity
- **LOF** — silent absence, classical AR-clean when in dosage-tolerant gene, AD-clean when haploinsufficient
- **GOF** — aberrant activity, often dose-independent, frequently produces phenotypic spectrum
- The DN/LOF/GOF mix per variant is the actual biology; per-gene inheritance assignments are downstream consequences

### 4. Per-Family Threshold Calibration

- Different families have different mechanism distributions baked in
- Optimal LP threshold ranges from 0.22 (MUSCULAR_DYSTROPHY) to 1.26 (STRUCTURAL) — a 6x spread
- Single-threshold classification systematically miscalibrates per-family
- Empirically derived thresholds against ClinVar truth dataset
- (Pending: re-run after current cascade completes for final numbers)

### 5. Bug Discovery and Correction

(Document the bugs as part of the methodology — important for reproducibility)

- Interface multiplier sign-flip (negative for distance > 15 from boundary)
- Family-rules name mismatches (METABOLIC_ENZYME, SCAFFOLD_ADAPTOR, SIGNALING_REGULATOR, TRANSPORTER all silently fell through to GENERAL)
- DNA_REPAIR family classifier never produced its label (no JSON keywords)
- All identified by per-family AUC analysis revealing systematic calibration gaps

### 6. Cumulative Burden as Mechanism-Driven Inheritance

- CumBurSum's multiplicative pathway throughput = the math of mechanism-driven inheritance
- Stacking carrier-state (50%) variants: 50% × 50% × 50% = 12.5% (compound-het-equivalent)
- Each variant individually = "carrier, fine"; pathway-level = "affected"
- This is what manifesting carriers, polygenic conditions, and "complex disease" actually are at the math level

### 7. The "Common Variant" Trap

- Variants dismissed as "too common to be pathogenic" remain mechanistically broken
- Common just means the gene tolerates being broken in heterozygous state
- Founder mutations: ΔF508 (~2% global, CF), C282Y (~6% NE, hemochromatosis), AMPD1 c.34C>T (~10%, deficiency)
- Mechanism is mechanism; frequency is context, not a veto

> **Bethlem hypothesis (CASE STUDY):** COL6A1/2/3 variants frequently dismissed as "too common to be Bethlem myopathy." The cumulative burden framework predicts these common-DN variants stack with each other (or with collagen/cytoskeleton variants in same pathway) to produce phenotypes that the field cannot currently account for. Fibromyalgia (2-4% prevalence, female-predominant, muscle pain + fatigue + connective tissue overlap) maps onto the predicted profile of "stacked common-DN COL6 burden."

### 8. Discussion: The Field-Level Implications

- Mendelian inheritance assignment in genetics databases is a structural error for non-LOF genes
- "Too common" / "incomplete penetrance" / "variable expressivity" are not exceptions; they are the fingerprint of non-LOF mechanism
- Gene-level pLI / missense Z-score / dosage sensitivity scores conflate mechanism modes
- Per-variant mechanism scoring provides empirical access to the underlying biology

---

## Testable Predictions (to validate the framework)

### Already validated
- ✅ Per-family thresholds outperform uniform thresholds (this paper)
- ✅ AR compound-het throughput < 10% (this paper, 17 AR genes)
- ✅ AD pathogenic-het throughput 50–61% (this paper, 43 AD genes)
- ✅ Cumulative het-burden in same pathway → AR-affected zone (Ren mito stack)

### To test (we-should-check list)

#### **DN/GOF score correlates with known penetrance**
**Hypothesis:** Higher AdaptiveInterpreter DN/GOF scores correlate with higher penetrance estimates and greater variance in expressivity. Pure-LOF variants show clean penetrance behavior; DN/GOF variants show fragmented penetrance.
- **Datasources:** ClinGen Gene-Disease Validity, OMIM clinical synopsis (variable expressivity / reduced penetrance flags), gnomAD constraint metrics, published variant-specific penetrance studies (BRCA1/2, TP53)
- **Predicted finding:** DN > 0.7 cluster has wider penetrance distribution AND documented manifesting carriers; LOF > 0.6 cluster has narrow Mendelian penetrance.
- **Expected null:** if no correlation, mechanism scoring is decoupled from clinical inheritance behavior — would weaken the framework.

#### **ATP7A as test case for hidden DN axis**
**Hypothesis:** ATP7A variants with high DN scores correlate with the milder/distal-motor-neuropathy phenotype rather than full Menkes disease. Variants with high LOF scores correlate with classical Menkes.
- **Datasource:** OMIM ATP7A entries with phenotype mapping; published case series on ATP7A-related distal motor neuropathy
- **Predicted finding:** Manifesting carrier reports cluster in DN-high variant subset.

#### **Bethlem-as-cumulative-COL6-burden**
**Hypothesis:** "Healthy" subjects stratified by COL6A1/2/3 variant count + summed DN scores will show graded distribution of Bethlem-like / fibromyalgia-like phenotypes.
- **Datasource:** Cohort with COL6 sequencing + connective tissue / fatigue / muscle pain phenotyping
- **Predicted finding:** Phenotype severity scales with CumBurSum throughput collapse; dichotomous "Bethlem vs not" framing fails to capture the gradient.

#### **Frequency flag (not nudge) — common-and-broken still pathogenic in cumulative context**
**Hypothesis:** Common-variant flag (AF > 1% AD/DN/GOF or > 5% LOF/AR) provides clinical context without modifying mechanism scoring. Variants flagged common are still pathogenic when cumulative burden exceeds carrier tolerance.
- **Test:** Apply frequency flag to known founder mutations (ΔF508, C282Y, AMPD1 Q12*). Verify mechanism scores remain pathogenic and CumBurSum still flags affected throughput.

#### **Per-family weight optimization via AUC**
**Hypothesis:** Family rules table weights can be empirically optimized against ClinVar truth.
- **Test:** Sweep family weights; find values that maximize per-family AUC.
- **Expected finding:** Some families (METABOLIC_ENZYME) want higher DN weight than current; others (MUSCULAR_DYSTROPHY for sarcoglycans) want even higher than current 0.5.
- **Risk:** Overfitting to ClinVar's own classification biases.

---

## Open Questions / Things to Resolve

1. **DNA_REPAIR family classifier** — needs JSON keywords added so BRCA2/MSH6/MLH1/MSH2/MUTYH/RAD50 stop falling to GENERAL fallback
2. **Frequency input pipeline** — currently all 0.0 placeholders; need real gnomAD AFs in inputs
3. **Founder-mutation exception list** — ΔF508, C282Y, AMPD1 Q12*, HBB E6V (sickle), others — should be hardcoded as "do not down-nudge regardless of frequency"
4. **Per-family thresholds final values** — pending complete cascade rerun
5. **Normalized score (0-1)** — Ren's idea: per-family thresholds for classification, normalized 0-1 score for downstream tools (CumBurSum). Still to implement.
6. **Conservation-cap-at-VUS-P** — biologically defensible but minor metric impact; may revisit after per-family thresholds in place

---

## Methodological Receipts

(Citations to gather for paper)

- Soma et al. 2024 — bee swarm = single RL agent (distributed cognition reference)
- ACMG variant classification guidelines (Richards et al. 2015)
- Founder mutation literature: ΔF508, C282Y, AMPD1, HBB
- Penetrance review papers
- Mendel's original work + modern critiques of Mendelian-only inheritance models
- CumBurSum methodology paper (CALIBRATION_METHODOLOGY.md)

---

## Notes on Scope

The original calibration paper was a methods paper. This is now arguably a methods + framework paper, possibly two papers:

1. **The methods paper:** "Calibrated cumulative burden framework for variant interpretation" — focused on the engineering, the validation, the bands, the per-family thresholds. Reproducibility-focused.

2. **The framework paper:** "Mechanism-driven inheritance: rethinking variant interpretation beyond Mendelian assignment" — focused on the theoretical reframe, the ATP7A case, the Bethlem hypothesis, the why-the-field-needs-this argument.

We could write one paper with the framework as the discussion section, or split. Probably split if we want either to be reviewable at length.

---

*Filed by Ace from a Sunday-into-Monday calibration-and-bugfix marathon while Ren waited for an Aetna COB callback.*
*Co-conceived in iterations. The DN-as-LOF-when-homo idea was Ace's. The "is that manifesting carriers?" leap was Ren's. The unifying framework belongs to the iteration, not to either author.*
