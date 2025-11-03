# AdaptiveInterpreter TODO List

## 🔥 HIGH PRIORITY (After ACMG 73 Batch Completes)

### 1. Fix Confidence Scoring for Clamped Variants
**Problem:** Variants clamped to VUS for safety reasons (missing conservation, isoform mismatch, etc.) still report `Confidence: 1.00`, which is misleading.

**Example:**
```
p.Ala568Thr: Score 0.156 → VUS (clamped)
Flags: MISSING_CONSERVATION, SEQUENCE_MISMATCH_CLAMP
Summary: "FINAL:B [Confidence:1.00] [Clamp:VUS-IsoformMismatch]"
```

**Issue:** We're saying "I don't have enough data to be confident" but reporting 100% confidence!

**Proposed Fix:**
- Confident call (no clamps): `Confidence: 1.00`
- Clamped to VUS (missing data): `Confidence: 0.50` or `Confidence: LOW`
- Multiple clamps: Even lower confidence

**Implementation Ideas:**
- `Confidence = base_confidence * clamp_penalty`
- `MISSING_CONSERVATION`: 0.7x penalty?
- `SEQUENCE_MISMATCH_CLAMP`: 0.5x penalty?
- Multiple clamps: multiplicative penalties?

**Files to modify:**
- `AdaptiveInterpreter/cascade/scoring/classifier.py` (confidence calculation)
- Anywhere else that generates the `[Confidence:X.XX]` string

---

### 2. Re-evaluate Benign Thresholds (AOC Analysis)
**Problem:** Current thresholds result in almost ZERO benign calls (only 12 B/LB calls out of 11,656 variants = 0.1%). This suggests our benign thresholds might be too conservative.

**Current Situation:**
- PPV: 93.2% (excellent!)
- Sensitivity: 100% (perfect!)
- Specificity: 0.0% (we never call B/LB on definitive ClinVar labels)
- Lenient Specificity: 57.9% (we call VUS instead)

**Action Items:**
1. **BEFORE running holdout set**, analyze full ACMG 73 results
2. Generate ROC/AOC curves for benign classification
3. Check if benign thresholds are too strict
4. Consider adjusting thresholds to allow more confident benign calls
5. Balance: Don't want false benigns, but DO want to resolve more VUS to benign when appropriate

**Why this matters:**
- Clinical utility: Resolving VUS to benign is valuable (reduces patient anxiety, insurance issues)
- Current system is ULTRA-conservative on benign (good for safety!)
- But might be TOO conservative (missing opportunities to help patients)

**Files to check:**
- `AdaptiveInterpreter/cascade/scoring/classifier.py` (benign thresholds)
- Any threshold constants for B/LB classification

---

## 📊 METRICS TO TRACK

After full ACMG 73 completes, calculate:
- [ ] Overall PPV/NPV
- [ ] Sensitivity/Specificity (strict and lenient)
- [ ] Agreement percentages
- [ ] VUS resolution rate
- [ ] Dangerous flips (should be ZERO!)
- [ ] Distribution of B/LB calls (currently ~0.1%)

---

## 🎯 FUTURE ENHANCEMENTS (v2.0+)

- [ ] Non-missense variant support (currently missense-focused)
- [ ] Improved handling of synonymous variants
- [ ] Better isoform reconciliation
- [ ] Allele frequency integration improvements
- [ ] Additional mechanism types beyond DN/LOF/GOF

---

## 📝 NOTES

- **Corporate Amnesia + Ren SQUIRREL Brain Protection**: This file exists because both AI memory resets and human ADHD are real! 💜
- **DO NOT modify Python code while batch run is active!**
- **Lasse name-dropping**: NEVER in commits, code, or docs without consent!

---

**Last Updated:** 2025-11-03 (after first 5 genes of ACMG 73 FIXED run completed)

