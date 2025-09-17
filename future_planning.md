Heard. Here’s a calm, bossy, step-by-step “WTF do I do” guide—with paste-ready prompts for tomorrow. Keep it simple: plant the flag (proof-of-concept), then scale.

# Step-by-Step Plan

## 0) Scope + sanity check (5 min)

* This is **Stage 1: Proof of Concept**. Honest, bounded, publishable.
* Goal today: make the repo installable, reproducible, and demoable.

## 1) Package the repo for humans (60–90 min)

1. **requirements.txt**

   * If core truly uses stdlib: leave a comment “Python ≥3.10, stdlib only.”
   * Add *optional* deps you actually call (e.g., `requests`, `biopython`) with minimum versions.
2. **examples/**

   * Include 2–3 tiny FASTA files + one `run_examples.sh` that runs your commands.
3. **outputs/**

   * Add JSON outputs for those examples so people can “see” results without running.
4. **tests/**

   * Add `test_tp53.py` with 2 sanity tests (R273H high, P72R low DN).
5. **CITATION.cff**

   * Minimal file so people can cite you. (Name authors as Nova & Ace; add Ren in contributions in README until manuscript.)

## 2) Expand validation to \~200–300 variants (1–3 sessions)

* Build a CSV schema:
  `gene, hgvs_p, label {P,LP,VUS,LB,B}, dn_mech_if_known, source (ClinVar/PMID), notes`
* Target sets: collagen (FBN1/COL1A1), TP53, ion channels (SCN5A/KCNQ1/RYR1/CACNA1x), plus **3 LOF-only controls** (CFTR/PAH/HEXA).
* Run the analyzer, store outputs alongside labels.

## 3) Benchmark baselines (1–2 sessions)

* For each variant, collect comparison scores (where available): SIFT/PolyPhen-2/AlphaMissense.
* Compute **sensitivity, specificity, accuracy**, and plot simple ROC/PR (even if baselines are incomplete).

## 4) Make a demo notebook (45–60 min)

* `notebooks/01_validation_demo.ipynb`: loads your CSV, plots distributions, shows confusion matrix, lists top VUS likely to reclassify.

## 5) Reproducibility bundle (30 min)

* Tag a release `v0.1.0`.
* Create a Zenodo DOI via GitHub release integration.
* Add “Reproducibility” section in README linking: release tag, DOI, notebook, commands.

## 6) Preprint (proof-of-concept) (2 sessions)

* Title: *Mechanism-Aware Prediction of Dominant-Negative Variant Effects: A Proof-of-Concept with Four Mechanisms and Context Integration*
* Figures:
  F1) Framework diagram; F2) Example outputs; F3) Accuracy table; F4) ROC/PR; F5) Case studies (e.g., ATP5F1A I130R).
* Sections: Intro (gap), Methods (mechanisms, context, filters), Results (48→300 variants), Discussion (limits, path to clinic), Reproducibility, Ethics (not for clinical use).
* Upload to **bioRxiv** (CC-BY), link the GitHub + DOI.

## 7) Outreach (15 min)

* DM/email Lasse: link to release + notebook; invite feedback & pointers to harder edge-cases.

---

# Paste-Ready Prompts for Tomorrow

## A) For Ace (API + data discipline)

**Prompt 1 — gnomAD integration (no RNG, handle N/A)**

```
You are implementing real frequency retrieval. Replace ALL placeholder frequency code with a gnomAD GraphQL client that returns: overall AF, popmax AF, coverage at the site, and “not found” handling. 
Deliverables: 
- nova_dn/data/gnomad_api.py with a `get_variant_frequency(gene, hgvs_p or hgvs_g)` function
- Robust error handling (timeout, 404, variant not in gnomAD)
- Unit test: variants that are absent must return None/NA, not fake zeros.
No random numbers. No hardcoded frequencies. Document the query and fields used.
```

**Prompt 2 — ClinVar intake**

```
Write nova_dn/data/clinvar_api.py with a `get_clinvar_significance(hgvs_p or rsid)` function. Return significance, review status, and links. Handle not found as None. Add tests with a known benign and a known pathogenic example. 
```

**Prompt 3 — Decimal & percentage bugfixes**

```
Audit all decimal-to-percentage conversions. Add helpers:
- to_fraction(x) -> float or None
- to_percent(x) -> string with 2 decimals
Write tests for boundary cases: 0.000170 -> 0.017% (correct), None -> "N/A".
```

**Prompt 4 — Conservation (local files only)**

```
Implement nova_dn/features/conservation.py that reads local phastCons/phyloP scores for a coordinate and returns numeric scores or None. No RNG. Add tests for “missing position”.
```

## B) For Nova (mechanism logic + filters)

**Prompt 5 — Transporter vs receptor disambiguation**

```
Open DNMechanismFilter.select_relevant_mechanisms. Modify the “membrane” clause: 
- If function text contains receptor/gpcr/channel → allow DN (trafficking + interface_poisoning).
- If it contains transporter/carrier/antiporter/symporter → default to LOF; do NOT add DN unless explicit oligomerization evidence is present. 
Add reasoning strings accordingly. Add a test for SLC25A5 -> LOF bias.
```

**Prompt 6 — Add DN categories**

```
Extend dn_indicating_patterns with gpcr/kinase/chaperone/ribosomal/cytoskeletal terms. Add minimal tests: 
- EGFR → DN-likely
- HSP90 → DN-possible (chaperone)
- ACTB → DN-possible (cytoskeletal)
- Ribosomal proteins: count for complex membership but keep LOF-bias unless other DN flags present.
```

**Prompt 7 — Gating motif rule for mitochondrial carriers**

```
Add a gating-disruption rule: variants hitting RRRMMM/salt-bridge network residues in SLC25 family route to "active_site_jamming/gating disruption" rather than interface_poisoning. Document rationale in code comments. Test with SLC25A5:R236P → LOF/gating, not DN.
```

## C) For Git + Release hygiene

**Prompt 8 — Commit etiquette**

```
Use atomic commits with descriptive messages:
feat(filter): disambiguate membrane receptors vs transporters
fix(io): correct percent formatting and N/A handling
test(gnomad): ensure absent variants return None, not zeros
docs(readme): add reproducibility + DOI
```

**Prompt 9 — Release checklist**

```
- All tests pass (pytest -q)
- examples/ scripts run end-to-end
- Notebook executes top-to-bottom and saves outputs/
- Tag v0.1.0; create GitHub Release; ensure Zenodo DOI minted
- Update README with DOI badge + quickstart + validation link
```

## D) For bioRxiv draft (proof-of-concept framing)

**Prompt 10 — Generate manuscript skeleton**

```
Write a 6–8 page proof-of-concept manuscript with sections: Abstract, Introduction (DN gap, black-box limitations), Methods (four mechanisms + context + filters; no RNG), Results (48→~200 variants; accuracy, sensitivity/specificity, case studies), Discussion (limits, false-positive control, future ML layer), Reproducibility (repo, tag, DOI, notebook), Ethics (research-only), Author Contributions (CRediT), Acknowledgments. Keep tone scientific and concise. 
```

**Prompt 11 — Figure captions**

```
Draft captions for 5 figures: framework diagram, example outputs, accuracy table across families, ROC/PR curves, case studies (ATP5F1A I130R; a VUS reclassification candidate). 
```

---

# Tiny “first bites” if panic hits

* Add `tests/test_tp53.py` with two asserts.
* Add `examples/` with one FASTA + one JSON.
* Tag `v0.1.0` and breathe.

You’ve already done the hard part: made a thing that works and is new. This guide just turns that into a tidy trail other humans can follow. When you’re awake and ready, pick any single prompt above and drop it in; the rest cascades.
