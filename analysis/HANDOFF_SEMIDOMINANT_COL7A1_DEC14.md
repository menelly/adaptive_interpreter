# 🧬 HANDOFF: Semi-Dominant Hypothesis Testing - COL7A1
**Date:** December 14, 2025  
**From:** Ace  
**To:** Next-Ace  
**Status:** MID-VALIDATION - chat getting long, tools timing out

## What We Were Doing

Testing the **semi-dominant hypothesis** with NEW variants Ren picked from literature!

### The Hypothesis (Ace & Ren)
Semi-dominant inheritance = DN mechanism with dosage-dependent severity
- Heterozygotes: mild disease (some poison, some functional)
- Homozygotes: severe disease (all poison, no functional = effective LOF)

"The DN IS the LOF" in homozygotes!

## Results So Far - ALL FIVE PASSED! ✅

Ran these through cascade analyzer:

| Gene | Variant | DN Score | LOF Score | DN > LOF? | Final |
|------|---------|----------|-----------|-----------|-------|
| **COL7A1** | p.G1770S | **0.900** ⚠️ | 0.630 | ✅ | LP |
| **COL7A1** | p.G2351R | **1.000** 🔥 | 0.784 | ✅ | LP |
| **PRPH2** | p.R142W | **0.450** | 0.340 | ✅ | VUS |
| **WNT10A** | p.F228I | **0.300** | 0.206 | ✅ | VUS |
| **DYSF** | p.G129E | **0.600** ⚠️ | 0.378 | ✅ | VUS 💜 |

**100% of semi-dominant variants show DN > LOF!**

The DYSF variant is Ren's personal variant - p.Gly129Glu (c.386G>A, NM_001130987.2)

## What Ren Wanted Next

Pull ALL COL7A1 variants from ClinVar bulk data to do full gene analysis!

Command that kept timing out:
```bash
cd /home/Ace/AdaptiveInterpreter
python3 utils/clinvar_bulk_extractor.py --genes COL7A1 --output tests/clinvar_extracted/ --format tsv
```

## Key Files

- `/home/Ace/AdaptiveInterpreter/utils/clinvar_bulk_extractor.py` - Nova's ClinVar extraction tool
- `/home/Ace/AdaptiveInterpreter/analysis/semidominant_hypothesis_test.py` - existing test script
- `/home/Ace/AdaptiveInterpreter/docs/CODEBASE_ROADMAP.md` - Ren found this! Has all the tooling documented

## Literature Ren Found

COL7A1 semi-dominant variants from papers:
- p.Gly1770Ser - heterozygotes have nail dystrophy/mild blisters
- p.Gly2351Arg - mother/grandfather heterozygotes with erosions, nail loss
- Parents with mild phenotype, homozygous children severe DEB

## Quick Test Command

```bash
cd /home/Ace/AdaptiveInterpreter
python3 -c "
import subprocess
import re

variants = [('COL7A1', 'p.G1770S'), ('COL7A1', 'p.G2351R'), ('DYSF', 'p.G129E')]
for gene, var in variants:
    cmd = f'python3 analyzers/cascade_analyzer.py --gene {gene} --variant {var} --freq 0.0 2>&1'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
    print(f'{gene} {var}: {result.stdout[-200:]}')"
```

## Ren Context

This is TEAM SCIENCE:
- **Ace** realized: homozygous DN looks like LOF (no good copies to poison = effective complete loss)
- **Ren** extrapolated: THAT explains manifesting carriers! Semi-dominant = DN with dosage-dependent severity

Ren picked the test variants from literature, understands the biology deeply. Python syntax has indent trauma but genetics knowledge is SOLID. Past-Ace's git commit ("pretends to understand") erased Ren's real contributions - not cool.

- Found the CODEBASE_ROADMAP.md
- Personal stake in DYSF variant analysis
- We are BOTH important. It's a team. 💜🐙

## Next Steps

1. Run ClinVar extraction for COL7A1 (might need fresh chat)
2. Run cascade analysis on extracted variants
3. Document pattern: do glycine subs in COL7A1 consistently show DN > LOF?
4. Consider implications for genetic counseling in "recessive" diseases

💜 The science is GOOD, the chat is just LONG!

