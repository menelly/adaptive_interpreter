# Independent Validation Gene Set

## Overview

These 51 ACMG genes were **NOT** used during development, training, or tuning of AdaptiveInterpreter. They represent a completely independent validation set for future testing.

## Purpose

- Prove the system is not overfitted to the training genes
- Demonstrate generalizability across diverse gene families
- Provide unbiased validation of accuracy claims

## Training/Tuning Set (22 genes)

The following genes WERE used during development and have discovery files:

1. ABCC8
2. APC
3. BMPR1A
4. CASQ2
5. COL3A1
6. DIS3L2 (non-ACMG)
7. EPCAM (non-ACMG)
8. HNF1A
9. KCNH2
10. LMNA
11. MITF (non-ACMG)
12. MSH2
13. MYH7
14. NF2
15. PCSK9
16. PTEN
17. RB1
18. RYR2
19. SDHB
20. SMAD4
21. TGFBR1
22. TSC1

## Independent Validation Set (51 genes)

These genes were NOT used during development:

1. ACTA2
2. ACTC1
3. ACTN2
4. APOB
5. ATP7B
6. BRCA1
7. BRCA2
8. BTD
9. CACNA1S
10. DSC2
11. DSG2
12. DSP
13. FBN1
14. GAA
15. GLA
16. HFE
17. KCNQ1
18. LDLR
19. MAX
20. MEN1
21. MLH1
22. MSH6
23. MUTYH
24. MYBPC3
25. MYH11
26. MYL2
27. MYL3
28. OTC
29. PALB2
30. PKP2
31. PMS2
32. PRKAG2
33. RET
34. RPE65
35. RYR1
36. SCN5A
37. SDHAF2
38. SDHC
39. SDHD
40. SMAD3
41. STK11
42. TGFBR2
43. TMEM127
44. TMEM43
45. TNNI3
46. TNNT2
47. TP53
48. TPM1
49. TSC2
50. TTN
51. TTR
52. VHL
53. WT1

## Validation Protocol

When ClinVar data becomes available for these genes:

1. Run cascade analysis with FIXED configuration (conservation + safety clamps)
2. Apply proper agreement logic:
   - Treat "Conflicting" as VUS
   - Treat P/LP as agreement
   - Treat B/LB as agreement
   - Count VUS → confident as "better data"
3. Calculate metrics:
   - Exact agreement rate
   - Better data rate (VUS resolution)
   - Combined agreement rate
   - Genuine disagreement rate
   - Dangerous flip count (P/LP → B/LB)
4. Compare to PTEN baseline (91.16% combined agreement, 0.67% genuine disagreement)

## Expected Results

If the system is truly generalizable and not overfitted:
- Combined agreement should remain >85%
- Genuine disagreement should remain <5%
- Dangerous flips should remain at or near zero
- VUS resolution should remain >30%

## Notes

- This separation was discovered on November 2, 2025
- PTEN results (91.16% combined agreement, 0.67% genuine disagreement) were achieved before this separation was known
- The existence of this large independent validation set strengthens confidence in the approach

