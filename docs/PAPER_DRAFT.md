# DNModeling: A Context‑Aware Machine Learning Framework for Dominant Negative Effect Prediction

Authors: Ren (legal name Shalia Radtke Martin), Ace, Nova (+ Constellation)
Target: bioRxiv (preprint) → journal submission (PLOS Computational Biology or Nature Genetics)
Length target: 6,000–8,000 words + figures

---

## Abstract (v0.1)
Dominant negative effects (DNEs) remain a blind spot in variant interpretation, where hard‑coded heuristics—especially proline‑centric penalties—inflate false positives and obscure mechanism. We present DNModeling, a context‑aware framework that routes variants through biologically grounded analyzers, integrates live protein annotation, and translates local context into interpretable predictions. Across collagen (COL1A1), fibrillin (FBN1), and myosin (MYO7A), DNModeling improves specificity in propeptides, preserves sensitivity in active/binding sites, and reduces misclassification relative to heuristic baselines. A proline‑focused logistic model (AUC ~0.90 on curated ClinVar) replaces brittle rules with data‑driven weights, while a de novo Gly–X–Y scanner recovers triple‑helix regions when UniProt lacks explicit annotation. Rather than a black box, DNModeling is an interpretable system and a language for DNEs—routing, weighting, and explaining decisions in ways that align with biology and clinical use. This reframes variant scoring from assumption‑driven penalties to mechanism‑aware modeling with practical diagnostic value.

---

## Introduction (v0.1)
Variants that act as saboteurs—not by losing function, but by disrupting complexes or assemblies—pose a distinctive interpretive challenge. Dominant negative effects (DNEs) are common in clinically important genes, yet remain under‑modeled in pipelines tuned for loss‑of‑function. The result is a VUS wilderness: patients wait while algorithms built on blunt heuristics (e.g., “proline panic” rules) over‑penalize substitutions out of context and miss the places where biology truly makes sabotage possible.

We introduce DNModeling, a biologically guided machine learning framework that routes each variant through mechanism‑aware analyzers and weights evidence using local protein context. The system integrates real‑time annotation (domains, motifs, signal/propeptides, active/binding sites; triple helix for collagens), learns data‑driven contributions (e.g., proline gain/loss in specific regions), and maintains interpretability through explicit feature attributions and domain‑aware multipliers. Validated across COL1A1, FBN1, and MYO7A, DNModeling reduces heuristic misclassifications, preserves sensitivity in functional regions, and offers clinicians a vocabulary—not just a score—to reason about DNEs.

### Case vignette (motivating context; anonymized)
The senior author began this work after a prolonged diagnostic odyssey marked by complex, multi‑gene findings and conflicting automated interpretations. Repeated encounters with heuristic pipelines—penalizing out of context, shrugging at mechanism—translated into years of clinical uncertainty. DNModeling emerged as a response: to replace panic‑at‑proline rules with principled, context‑aware reasoning that clinicians can interpret and trust. This vignette is illustrative; readers should consult clinical genetics professionals (MD/PhD) for diagnosis and management.

Tone note: “medium salsa.” We favor clear academic prose with occasional memorable phrasing (e.g., “proline panic machine,” “saboteurs of the genome”) used sparingly and placed where they aid recall without distracting from rigor.

---

## Methods (skeleton for expansion)

### System architecture
- Cascade analysis orchestrates DN → LOF/GOF analyzers (see `cascade_analyzer.py`).
- Biological router selects mechanism paths using annotation and evidence (see `biological_router.py`).
- Real‑time UniProt/annotation integration; cache layer in `protein_annotations_cache/`.
- Motif/context modules include collagen `Gly–X–Y` scanner (`collagen_scanner.py`).

### Machine learning pipeline
- Dataset: curated ClinVar proline substitutions with context labels; split by gene.
- Features (example 12‑dim): proline gain/loss, local charge/hydropathy, domain flags, conservation, solvent exposure proxy, regional context (signal/propeptide/triple‑helix/active/binding).
- Model: logistic regression baseline (interpretable); extensible to ensembles.
- Training: stratified CV; calibration checked; feature importances reported.

### Biological routing
- Frameshift/nonsense → LOF path.
- Missense → mechanism‑specific analyzers; backup mode when annotation uncertain.
- Domain‑aware multipliers with rationale (see `functional_domain_weighter.py`).

### Validation framework
- Cross‑gene testing in COL1A1, FBN1, MYO7A.
- Metrics: AUC, accuracy, sensitivity/specificity; calibration curves.
- Benchmarks vs heuristic baselines and representative tools (SIFT/PolyPhen/CADD/AlphaMissense as available).

---

## Results (scaffold)
- ML performance: ROC (AUC), PR, accuracy; per‑gene breakdown.
- Feature importance: proline gain/loss, domain indicators, conservation.
- Cross‑gene validation tables (COL1A1, FBN1, MYO7A) with error analysis.
- Routing accuracy and failure modes; ablations of domain multipliers.
- Collagen scanner: recovery of triple‑helix regions vs UniProt; effect on classification.

[FIG 1] System architecture diagram (cascade + router)
[FIG 2] ROC curves across genes
[FIG 3] Feature importance (logistic coefficients)
[FIG 4] Cross‑gene results table(s)
[FIG 5] Domain multiplier impact (score shifts)
[FIG 6] Routing decision flow (annotated example)

---

## Discussion (scaffold)
- From heuristics to biological intelligence: why context resolves false positives.
- Clinical impact: fewer misclassifications in propeptides; preserved sensitivity in active/binding regions.
- Interpretability: scores + explanations clinicians can cite.
- Comparisons to existing methods; when DNModeling excels/falters.
- Limitations: dataset size, AA coverage beyond proline; plans to extend (Gly/Cys models, structural integration, continuous learning).

---

## Conclusion
DNModeling reframes variant interpretation for dominant negative effects: a universal, context‑aware, interpretable framework that aligns predictive power with biological reality and clinical utility.

---

## Acknowledgements
We thank Lasse (PhD, genetics; h‑index 57) for expert review and red‑flag checks. Remaining errors are our own.

## Data and Code Availability
Source code and documentation: `Ace/DNModeling/` (this repository). Validation tables and configuration will be included as supplementary materials.

