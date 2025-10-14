# Supplementary Table 1 — Directional Agreement Logic (DAL)

This supplement formalizes the comparison framework used throughout the manuscript.

Definitions
- ClinVar normalization (CV): map to {P, LP, VUS-P, VUS, LB, B}
- Model normalization (AI): map to {P, LP, VUS-P, VUS, LB, B}
- Buckets
  - Positive bucket (POS): {P, LP, VUS-P}
  - Negative bucket (NEG): {VUS, LB, B}
- Special equivalence: VUS ↔ VUS-P are considered directionally equivalent (neutral) for agreement and lenient binary metrics.

Agreement categories (DAL)
- AGREE: CV and AI fall in the same bucket (POS or NEG), counting VUS↔VUS-P as AGREE
- BETTER_DATA: CV is VUS and AI is any definitive class in {P, LP, LB, B}
- DISAGREE: CV and AI fall in opposite buckets (POS vs NEG)

Binary metrics (reported in two modes)
- Strict clinical mode
  - Positive = {P, LP}
  - Negative = {B, LB}
  - Rows with VUS or VUS-P on either side are excluded from TP/FP/TN/FN tallies
  - Outputs: sensitivity, specificity, PPV, NPV
- Lenient directional mode
  - Positive = {P, LP, VUS-P}
  - Negative = {VUS, LB, B}
  - VUS↔VUS-P pairs are excluded (neutral) to avoid inflating FP
  - Outputs: sensitivity, specificity, PPV, NPV

Pre-processing notes
- Protein-synonymous (p.=, p.XnX, p.XxxnXxx) entries are excluded from protein-effect metrics and reported separately; splice-driven pathogenic synonymous variants are handled outside protein-effect analyses.

Rationale
- DAL separates true model error (DISAGREE) from database-driven ambiguity (VUS) and captures cases where the model provides strictly more actionable information (BETTER_DATA). Reporting both strict and lenient binaries mirrors clinical practice (actionable vs directional).

