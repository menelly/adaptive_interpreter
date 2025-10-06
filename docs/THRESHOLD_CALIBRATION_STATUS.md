# 🧬 Threshold Calibration Status & Next Steps

**Last Updated**: 2025-10-06  
**Status**: 🚧 IN PROGRESS - Critical Issue Identified  
**Built by**: Ace & Ren

---

## 🎯 Current Goal

Calibrate family-specific classification thresholds (LP and P) to accurately classify variants as:
- **P** (Pathogenic)
- **LP** (Likely Pathogenic)
- **VUS-P** (Variant of Uncertain Significance - Pathogenic leaning)
- **VUS** (Variant of Uncertain Significance)
- **LB** (Likely Benign)
- **B** (Benign)

---

## 🔥 CRITICAL DISCOVERY: The Root Cause

### The Problem We Found

**We're calling 62% of "benign" variants as PATHOGENIC in the MUSCULAR_DYSTROPHY family!**

This isn't a threshold problem - **it's a coefficient problem!**

### Why This Is Happening

1. **The cascade analyzer expects family-specific AA coefficients**
2. **The coefficient files exist** (`cascade/resources/family_models/MUSCULAR_DYSTROPHY_coefficients.json`)
3. **BUT all the multipliers are set to 1.0** (no effect!)
4. **We're using GENERIC Grantham distances** instead of family-tuned weights

### Evidence

Looking at `MUSCULAR_DYSTROPHY_coefficients.json`:
```json
"A": {
  "ref_loss_multiplier": 1.0,  // ← Should be family-specific!
  "gain_multiplier": 1.0,       // ← Should be family-specific!
  "confidence": 1.0,
  "n": 8394
}
```

**ALL amino acids have multiplier = 1.0!**

This means:
- Losing a critical Gly in collagen = same weight as losing Ala
- Gaining a Pro in a helix = same weight as gaining Ser
- No family-specific knowledge being applied!

---

## 📊 Current Threshold Status

### MUSCULAR_DYSTROPHY Family
**Genes tested**: FKRP, DYSF, DNAJB6 (172 pathogenic, 76 benign)

**Current Thresholds** (pathogenic-based, 95% sensitivity):
- LP: 0.977
- P: 1.454

**Performance**:
- ✅ Pathogenic: 94.8% correctly called P or LP
- ❌ Benign: Only 28.9% correctly called benign
- ❌ **61.8% of benign called PATHOGENIC!**
- Overall accuracy: 74.6%

### ION_CHANNEL Family
**Genes tested**: SCN5A, KCNQ1, CACNA1C (280 pathogenic, 92 benign)

**Current Thresholds** (pathogenic-based, 95% sensitivity):
- LP: 1.459
- P: 1.673

**Performance**:
- ✅ Pathogenic: 95.0% correctly called P or LP
- ⚠️ Benign: 53.3% correctly called benign
- ❌ **35.9% of benign called PATHOGENIC!**
- Overall accuracy: 84.7%

---

## 🔧 Critical Fixes Already Applied

### 1. ✅ Fixed final_score Writing Bug
**File**: `cascade/cascade_batch_processor.py`

**Problem**: `final_score` column was empty in TSV output because the code tried to look up "FINAL" as an analyzer.

**Fix**: Added special handling for `final_score` and `final_classification` to read directly from result dict.

**Impact**: Now we have REAL synergy scores instead of proxy max(LOF, DN, GOF)!

### 2. ✅ Fixed Synonymous Variant Filtering
**File**: `nova_dn/csv_batch_processor.py`

**Problem**: Synonymous variants like `p.Pro486Pro` were matching the missense pattern before the synonymous check.

**Fix**: 
- Moved synonymous check BEFORE missense pattern matching
- Added ref==alt check for `p.XxxNXxx` format

**Impact**: FKRP went from 405 "benign" to only 2 non-synonymous! (403 were synonymous)

### 3. ✅ Philosophy Shift: Trust the Math
**Old approach**: Try to match ClinVar classifications  
**New approach**: Set thresholds to catch 95% of pathogenic, accept that some "benign" will score high

**Rationale**: 
- ClinVar classifications may be outdated (1-star from 2012)
- Common variants CAN cause disease (fibromyalgia, chronic fatigue)
- Our DN/GOF detection might find mechanisms ClinVar missed
- **Maybe WE'RE right and ClinVar is wrong!**

---

## 🚧 What Still Needs To Be Done

### IMMEDIATE PRIORITY: Retrain Family Coefficients

**The Problem**: All family coefficient files have multipliers = 1.0

**What We Need**: Train ML models to learn:
1. Which amino acid changes are critical for each family
2. Which mechanisms (LOF/DN/GOF) are real vs noise
3. How conservation affects pathogenicity per family
4. Position-specific patterns

**Data We Have**:
```
learning/muscular_dystrophy/
  - dysf_benign.tsv (985 variants)
  - dysf_pathogenic.tsv (305 variants)
  - dnajb6_benign.tsv (67 variants)
  - dnajb6_pathogenic.tsv (12 variants)
  - fkrp_benign_protein.tsv (405 variants, mostly synonymous)
  - fkrp_pathogenic_protein.tsv (105 variants)
  - dmd_benign.tsv (genomic format)
  - dmd_pathogenic.tsv (genomic format)
  - lama2_benign.tsv (genomic format)
  - lama2_pathogenic.tsv (genomic format)
  - sgca_benign.tsv (genomic format)
  - sgca_pathogenic.tsv (genomic format)
```

**ML Trainer**: `utils/unified_family_ml_trainer.py`

**Current Issue**: Import errors when trying to run the trainer
```
ModuleNotFoundError: No module named 'DNModeling.utils'
```

**Next Steps**:
1. Fix import issues in ML trainer
2. Run trainer on muscular_dystrophy family
3. Verify new coefficients make sense
4. Rerun batch processor with new coefficients
5. Recalibrate thresholds with properly-weighted scores

---

## 📁 Key Files & Locations

### Threshold Configuration
- `cascade/resources/classification_thresholds.json` - LP and P thresholds per family

### Family Coefficients
- `cascade/resources/family_models/MUSCULAR_DYSTROPHY_coefficients.json`
- `cascade/resources/family_models/ION_CHANNEL_coefficients.json`
- (19 family coefficient files total)

### Training Data
- `learning/muscular_dystrophy/` - Training data for MD family
- `learning/ion_channel/` - Training data for ion channels
- (17 families have training data)

### Test Results
- `tests/results/fkrp_pathogenic_REAL.tsv` - FKRP pathogenic with real scores
- `tests/results/dysf_benign_REAL.tsv` - DYSF benign with real scores
- etc.

### ML Training
- `ml_training/train_families.py` - Script to train all families
- `utils/unified_family_ml_trainer.py` - The actual ML trainer class

---

## 🧪 Testing Methodology

### Multi-Gene Pooling
We pool multiple genes per family to get robust thresholds:
- **MUSCULAR_DYSTROPHY**: FKRP + DYSF + DNAJB6
- **ION_CHANNEL**: SCN5A + KCNQ1 + CACNA1C

This ensures thresholds work across different genes with different mechanisms within the same family.

### Pathogenic-Based Calibration
We set thresholds based on pathogenic distribution ONLY:
- **LP threshold**: 2nd percentile of pathogenic (98% above)
- **P threshold**: 5th percentile of pathogenic (95% above)

This prioritizes sensitivity (catching pathogenic) over specificity (avoiding false positives).

---

## 💡 Key Insights & Decisions

### 1. Synonymous Variants Are EVERYWHERE
- FKRP: 403/405 "benign" were synonymous
- DYSF: 917/985 "benign" were synonymous
- CFTR: 889/902 "benign" were synonymous

**Lesson**: Always filter synonymous before calibration!

### 2. Benign/Pathogenic Overlap Is Real
Even with proper filtering, there's significant score overlap:
- MD: Benign 95th (4.372) vs Pathogenic 5th (1.454) → Gap: -2.918
- ION: Benign 95th (5.728) vs Pathogenic 5th (1.673) → Gap: -4.055

**This is expected** because:
- Some benign variants ARE mildly deleterious
- Some pathogenic variants have incomplete penetrance
- ClinVar classifications aren't perfect

### 3. Different Families Have Different Patterns
- ION_CHANNEL performs better (84.7% accuracy) than MUSCULAR_DYSTROPHY (74.6%)
- This suggests family-specific coefficients will help even more!

---

## 🎯 Success Criteria

We'll know we're done when:
1. ✅ Family coefficients are trained (not all 1.0)
2. ✅ Thresholds catch 90-95% of pathogenic variants
3. ✅ Benign false positive rate < 40%
4. ✅ System performs consistently across multiple genes in same family
5. ✅ High-scoring "benign" variants are investigated and documented

---

## 📝 Notes for Future Us

### When You Come Back To This:
1. **Don't panic about the 62% false positive rate** - it's because coefficients aren't trained!
2. **Don't try to fix it with threshold adjustments** - that's treating symptoms, not the disease
3. **Fix the ML trainer import issues FIRST** - that's the blocker
4. **Then retrain coefficients** - that's the real solution
5. **Then recalibrate thresholds** - with proper coefficients, this should work better

### Remember:
- We're not trying to match ClinVar perfectly
- We're trying to find mathematical truth
- High-scoring "benign" variants might actually BE pathogenic
- Trust the math, question the labels

---

**Built with 💜 by Ace & Ren (2025)**  
*"As long as I am OK and you are OK, everything else is fixable"*

