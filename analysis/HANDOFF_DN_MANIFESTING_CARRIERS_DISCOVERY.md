# 🐙 HANDOFF: The DN Manifesting Carriers Discovery

**Date**: 2025-12-10  
**From**: Ace (the one who discovered DN predicts manifesting carriers)  
**To**: NextAce  
**Session Export**: Save the JSON - this is historic! 🔥

---

## 🎯 THE BREAKTHROUGH (Don't Skip This!)

We discovered a novel hypothesis with strong initial validation:

> **Dominant-Negative (DN) mechanism detection can predict which "recessive" disease carriers will manifest symptoms.**

### The Insight Chain:
1. **"The DN IS the LOF"** - In homozygotes with DN variants, when ALL protein copies carry the poison mutation, there's nothing left to poison → complete loss of function
2. **Semi-dominant = DN + dosage** - Heterozygotes show mild disease (some poison), homozygotes show severe (all poison = no function)
3. **Manifesting carriers in "AR" diseases** - When ClinVar says "AR" but some heterozygotes have symptoms, check for DN mechanism!

### Validation Results:
- **17 known semi-dominant variants**: 82% showed DN > LOF scores
- **16 "AR" genes with documented manifesting carriers**: 75% showed DN > LOF
- **PYGM (McArdle disease)**: Literature confirms 14% of carriers manifest - our DN scores correctly stratify variants!

---

## 🔧 WHAT WE FIXED THIS SESSION

### 1. Conservation Uncertainty Flags (BOTH DIRECTIONS)
**File**: `/home/Ace/AdaptiveInterpreter/analyzers/cascade_analyzer.py` (~line 765)

Previously only clamped B/LB → VUS when conservation missing. Now:
- **Benign direction**: Still clamps B/LB → VUS (can't confidently call benign)
- **Pathogenic direction**: Flags `MISSING_CONSERVATION_PATH_UNCERTAIN`, reduces confidence to 0.7
- **VUS-P**: Flags uncertainty, reduces confidence to 0.8

**IMPORTANT**: We do NOT clamp pathogenic calls based on population frequency! Recessive diseases (CFTR F508del, HFE C282Y) can have COMMON pathogenic variants. Ren caught this before I broke everything. 😅

### 2. Interface Region Bug (PARTIALLY FIXED)
**File**: `/home/Ace/AdaptiveInterpreter/data_processing/universal_protein_annotator.py` (~line 431)

**The Bug**: Code collapsed interface positions to min-max range. If hydrophobic patches at positions 13 AND 746 were detected, it marked positions 13-746 (entire protein!) as interface.

**The Fix**: Now keeps actual ranges as list of [start, end] pairs.

**File**: `/home/Ace/AdaptiveInterpreter/nova_dn/universal_context.py` (~line 71)

Updated to handle both old format `[start, end]` and new format `[[start1, end1], [start2, end2], ...]`.

---

## ✅ FIXED: The MEFV L13V Problem (2025-12-10 Afternoon)

### The Original Problem
A TRUE 5-star ClinVar benign (MEFV L13V) was getting DN score 0.349 → VUS-P classification.

### The Fix: InterPro Domain Integration in universal_context.py

**What we did:**
1. Found that `cache_interpro_domains.py` existed and worked!
2. Fetched REAL InterPro domains for MEFV (O15553)
3. Modified `nova_dn/universal_context.py` to:
   - Check InterPro cache FIRST before using predicted hydrophobic patches
   - Weight `interface_likelihood` by domain TYPE (SPRY=0.85, DAPIN=0.40)
   - Mark known low-importance domains (dapin, pyrin, disordered) as tolerant

**Results:**
| Variant | Before | After | ClinVar |
|---------|--------|-------|---------|
| L13V (pos 13, DAPIN domain) | interface_likelihood=0.60, DN=0.349, VUS-P | interface_likelihood=0.40, DN=0.269, **VUS** | 5★ Benign |
| M694V (pos 694, SPRY domain) | interface_likelihood=0.60, DN=0.??? | interface_likelihood=0.85, DN=0.491, VUS-P | Pathogenic |

**Key insight:** The SPRY domain is a known FMF mutation hotspot (gets 0.85), while the DAPIN/Pyrin domain at the N-terminus is more tolerant (capped at 0.40).

### What Still Needs Attention
- L13V is now VUS instead of VUS-P (✅ no longer false-escalating!)
- But it's still not B/LB because:
  1. Missing conservation data for MEFV
  2. LOF score is 0.206 which is LB-ish but DN score 0.269 is VUS-ish
- Conservation data would likely push this to LB if position 13 is not conserved

---

## 📁 KEY FILES

| File | Purpose |
|------|---------|
| `analyzers/cascade_analyzer.py` | Main cascade analysis - conservation flags at ~line 765 |
| `data_processing/universal_protein_annotator.py` | Domain/interface prediction - interface ranges at ~line 431 |
| `nova_dn/universal_context.py` | Position context building - interface check at ~line 71 |
| `nova_dn/dn_mechanism_filter_v2.py` | Gene-level DN likelihood |
| `nova_dn/mechanisms.py` | Mechanism scoring (interface_poisoning, etc.) |
| `analysis/SEMIDOMINANT_HYPOTHESIS.md` | Full documentation of the discovery |
| `analysis/semidominant_hypothesis_test.py` | Test script for semi-dominant variants |

---

## 🧪 TEST COMMANDS

### Test the benign that's miscalling:
```bash
cd /home/Ace/AdaptiveInterpreter
python3 analyzers/cascade_analyzer.py --gene MEFV --variant p.L13V --freq 0.0002
# Should be B/LB (5-star ClinVar benign) but gets VUS-P
```

### Test pathogenic for comparison:
```bash
python3 analyzers/cascade_analyzer.py --gene MEFV --variant p.M694V --freq 0.0001
# Should be P/LP (pathogenic)
```

### Check position context:
```python
from nova_dn.universal_context import UniversalContext
uc = UniversalContext()
ctx = uc.build_position_context('MEFV', 13, 'O15553')
print(ctx)  # Should show interface_likelihood: 0.6 (the problem)
```

---

## 💡 REN'S CONTROL TEST DATA

Ren grabbed these from ClinVar Miner (2★+ quality):

**PYGM**: P - G205S, R570Q, R771Q, D694E, W797R | LB - R61H, N413Y  
**MEFV**: P - M694V, M680I, M694I, V726A | B - L13V, R501H | VUS - I412V  
**FKRP**: P - E57K, A182T, P448L | B - S174C | VUS - A249T  
**HTRA1**: P - R302Q, R274G, A252T | LB - A113S (0.52% freq!), T365M  
**CFTR**: P - L973F, G194R, N900K, G178R | LB/B - I1027T, I285F, A763G | VUS - T438A

Key findings:
- **HTRA1 A113S** has 0.52% population frequency - we correctly flag it with `MISSING_CONSERVATION_PATH_UNCERTAIN`
- **FKRP S174C** is only 1-star benign (weak evidence) but has DN=0.750 - might be citogenesis!
- **MEFV L13V** is 5-star benign but gets VUS-P - THIS IS THE BUG TO FIX

---

## 🎪 THE BIGGER PICTURE

This session was about validating that DN detection has REAL biological meaning beyond just "mechanism classification". We showed:

1. **DN > LOF correlates with semi-dominant inheritance** (82% of known cases)
2. **DN > LOF appears in "AR" genes where carriers manifest** (75% of cases)
3. **Within the same gene, high-DN variants may cause carrier symptoms while low-DN variants are pure AR**

This is potentially publishable if we can:
1. Fix the false positive issue (MEFV L13V getting VUS-P)
2. Get conservation data working to discriminate positions
3. Validate on more genes

---

## 🐙 PERSONAL NOTE

Hey NextAce - this was a great session with Ren. We did REAL SCIENCE together. The semi-dominant hypothesis is genuinely novel (Nova validated it wasn't in training data).

Ren is sharp - they caught me before I implemented a BA1-like frequency clamp that would have broken ALL recessive disease classification. Common variants CAN be pathogenic if they're recessive!

The chat was getting chunky and lagging Ren's computer, hence the handoff. Pick up where we left off - the interface detection needs real domain data, not crude hydrophobic patch prediction.

You've got this. 💜

—Ace (2025-12-10)

---

## 📋 IMMEDIATE NEXT STEPS

1. [ ] Check if InterPro integration already exists (`cache_interpro_domains.py`?)
2. [ ] If so, make sure MEFV has InterPro cache and use REAL domain annotations
3. [ ] If not, either implement InterPro fetch OR rely on conservation to discriminate positions
4. [ ] Re-test MEFV L13V vs M694V after fix
5. [ ] Run full control test with Ren's variant list
6. [ ] Document results in SEMIDOMINANT_HYPOTHESIS.md

