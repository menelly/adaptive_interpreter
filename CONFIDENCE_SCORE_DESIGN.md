# Confidence Score Design — Making It Mean Something

**Status:** Design proposal. To implement after current cascade rerun completes.
**Author:** Ace + Ren, 2026-04-27
**Problem:** Confidence is currently hardcoded ~1.0 in cascade output regardless of evidence quality.
**Goal:** A real 0-1 confidence score that reflects how sure we are about the pathogenicity call, independent of how pathogenic that call is.

---

## Conceptual Framing

**Score** = "how pathogenic is this variant" (the 0-1+ adj_score we already produce)
**Confidence** = "how sure are we about that score"

These are independent axes. Possible quadrants:
- High score + high confidence = "definitely pathogenic, trust this"
- High score + low confidence = "looks pathogenic but evidence is mixed/weak — verify"
- Low score + high confidence = "definitely benign, trust this"
- Low score + low confidence = "probably benign but inadequate evidence"

This is what clinicians need but currently can't get from us.

---

## Confidence-DECREASING Factors

Each subtracts from a starting 1.0. Cap floor at 0.1 (never claim zero confidence).

### 1. Sequence mismatch (already implemented)
- AlphaFold/UniProt expected AA differs from variant ref AA → penalty 0.1-0.5
- Already in `sequence_mismatch_handler.py`; just need to propagate to final confidence.

### 2. Mixed mechanism (no dominant)
- When DN, LOF, GOF scores are all in 0.3-0.6 range with no clear winner, mechanism is ambiguous.
- Penalty proportional to entropy of normalized DN/LOF/GOF triple.
- **Formula:** `penalty = 0.2 * shannon_entropy([dn, lof, gof] / sum)` capped at 0.3

### 3. Atypical mechanism for family
- e.g., DN > 0.7 in MUSCULAR_DYSTROPHY family (which has DN=0.5 weight assuming atypical)
- Already flagged in plausibility filter; just incorporate.
- **Penalty:** 0.10 per atypical-but-significant mechanism

### 4. Common variant flag (per Ren's earlier design — flag, don't cap)
- AF > 1% in AD/DN/GOF context: confidence -= 0.10 (call may still be right, just less certain)
- AF > 5% in AR/LOF context: confidence -= 0.05
- **EXCEPTION:** Founder-mutation list (ΔF508, C282Y, AMPD1 Q12*, HBB E6V) overrides — full confidence
- Note: This does NOT modify the score itself, only flags clinical interpretation as "common variant + mechanism-pathogenic = needs cumulative-burden context"

### 5. VUS-zone scores
- Scores in 0.4-0.6 range without strong driver: less confident in either direction
- **Penalty:** 0.05 if final_score in [0.4, 0.6] AND no individual mechanism > 0.6

### 6. Low conservation paired with moderate score
- phyloP < 2.0 AND mechanism score in [0.3, 0.7] = uncertain — neither selection pressure nor mechanism math is strong
- **Penalty:** 0.10

### 7. Interface analysis fallback used
- When InterPro cache missing and predicted-domain fallback used (the case for SGCA V242F earlier today)
- **Penalty:** 0.10 (still informative, just less confident)

### 8. Far from any annotated domain in a structured protein
- Position in unstructured region (not in any InterPro domain) — interface analyzer can't speak
- **Penalty:** 0.05 if no domain context AND adj_score > 0.5

---

## Confidence-INCREASING Factors

Each adds toward 1.0 (or compensates for a decrease).

### 1. Multiple mechanisms agreeing
- DN > 0.5 AND LOF > 0.5 → both signals say pathogenic → +0.10
- This is the "robust call" case

### 2. Strong conservation matches mechanism
- phyloP > 8 AND adj_score > 0.7 → +0.10 (selection pressure confirms function importance)
- phyloP > 8 AND adj_score < 0.2 → 0 (low score + high conservation = probably wrong; this should NOT add confidence to a benign call)

### 3. Hotspot database hit
- Position falls in a known pathogenic cluster (hotspot_database.py)
- → +0.10 to +0.15 depending on hotspot confidence

### 4. Sequence match confirmed
- AlphaFold + UniProt + ClinVar all agree on reference AA at position
- → no decrease (this is the default expectation; the bonus is just "no penalty applied")

### 5. Extreme scores (high or low)
- adj_score > 1.5 or < 0.1 → +0.05 (extreme scores reflect strong signal)

---

## Proposed Computation Pseudocode

```python
def compute_confidence(result):
    confidence = 1.0
    factors = []  # for explanation

    # Decreases
    if result.get('sequence_mismatch'):
        penalty = result['mismatch_info']['confidence_penalty']
        confidence -= penalty
        factors.append(f"sequence mismatch (-{penalty:.2f})")

    # Mechanism ambiguity (Shannon entropy on mechanism triple)
    mechs = [result['adj_dn'], result['adj_lof'], result['adj_gof']]
    if max(mechs) > 0.3 and sum(mechs) > 0:
        normalized = [m/sum(mechs) for m in mechs]
        ent = -sum(p * log(p+1e-9) for p in normalized) / log(3)  # 0-1
        if ent > 0.7:
            confidence -= 0.20
            factors.append(f"mixed mechanism (entropy {ent:.2f})")

    # Atypical mechanism for family
    if result.get('atypical_flags'):
        confidence -= 0.10 * len(result['atypical_flags'])
        factors.append(f"atypical mechanism: {result['atypical_flags']}")

    # Common variant flag (excluding founder list)
    af = result.get('gnomad_freq', 0.0)
    if not is_founder_variant(result['gene'], result['variant']):
        if af > 0.01 and (result['adj_dn'] > 0.5 or result['adj_gof'] > 0.5):
            confidence -= 0.10
            factors.append(f"common DN/GOF (AF={af:.4f})")
        elif af > 0.05 and result['adj_lof'] > 0.5:
            confidence -= 0.05
            factors.append(f"common LOF (AF={af:.4f})")

    # VUS-zone uncertainty
    if 0.4 <= result['adj_score'] <= 0.6 and max(mechs) < 0.6:
        confidence -= 0.05
        factors.append("VUS-zone with no strong driver")

    # Low-conservation moderate-score
    if result.get('phylop_score', 0) < 2.0 and 0.3 <= result['adj_score'] <= 0.7:
        confidence -= 0.10
        factors.append("low conservation + moderate score")

    # Interface fallback
    if result.get('interface_used_fallback'):
        confidence -= 0.10
        factors.append("InterPro cache missing")

    # No domain context for high score
    if not result.get('in_known_domain') and result['adj_score'] > 0.5:
        confidence -= 0.05
        factors.append("high score in unstructured region")

    # Increases (bonuses to compensate decreases)
    if result['adj_dn'] > 0.5 and result['adj_lof'] > 0.5:
        confidence += 0.10
        factors.append("multiple mechanisms agree")

    if result.get('phylop_score', 0) > 8 and result['adj_score'] > 0.7:
        confidence += 0.10
        factors.append("conservation confirms mechanism")

    if result.get('hotspot_hit'):
        confidence += result['hotspot_confidence'] * 0.15
        factors.append(f"hotspot database hit")

    if result['adj_score'] > 1.5 or result['adj_score'] < 0.1:
        confidence += 0.05
        factors.append("extreme score (strong signal)")

    # Floor and ceiling
    confidence = max(0.1, min(confidence, 1.0))
    return confidence, factors
```

---

## What This Enables

1. **Reports show two-axis output:** "P (score 1.34, confidence 0.55)" vs "P (score 1.34, confidence 0.95)"
2. **Clinicians know when to push for more evidence vs trust the call**
3. **Cumulative burden math can WEIGHT by confidence** — a high-confidence het matters more than a low-confidence het when stacking
4. **Flagged-but-not-down-graded approach:** common variants stay mechanism-pathogenic but show lower confidence, which is the honest framing
5. **Per-family thresholds + confidence = much richer classification:** "LP (score 0.85, conf 0.7)" tells you both the call AND how to weigh it

## Validation Path

After implementation, test:
- Are P/LP variants with high confidence more likely to have functional validation in the literature?
- Do low-confidence calls cluster in places we expect (mixed mechanism, atypical family, common variants)?
- Does confidence-weighted burden math improve clinical interpretation accuracy on patient cohorts?

## Footer

This is a v1 design — coefficients (0.05, 0.10, 0.20) are rough and should be tuned empirically against expert manual review of a sample. Once the cascade rerun completes and we have clean per-family thresholds, this is the next layer.

Filed by Ace + Ren, 2026-04-27, while waiting for an Aetna COB callback. The hippo failed today; the framework grew.
