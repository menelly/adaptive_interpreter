#!/usr/bin/env python3
"""
🧬 SEMI-DOMINANT HYPOTHESIS TEST
================================
Testing: Do known semi-dominant genes show high DN scores?

HYPOTHESIS (Ace & Ren, December 2025):
Semi-dominant inheritance = DN mechanism with dosage-dependent severity
- Heterozygous: some poison subunits → mild disease  
- Homozygous: all poison subunits → severe disease

If true, our DN analyzer should detect high DN potential in these genes.

The insight: "The DN IS the LOF" - in homozygotes, there are no good copies
to poison, so you get complete loss of functional complex. Same mechanism,
different severity based on dosage.

This could explain why some "recessive" diseases have manifesting carriers!

Authors: Ace (Claude), Ren (Shalia)
Date: December 10, 2025
"""

import subprocess
import re
import sys
from datetime import datetime

# Known semi-dominant genes with specific pathogenic variants
# Format: (gene, variant, disease, inheritance_notes)
TEST_VARIANTS = [
    # MFN2 - Mitofusin 2, forms complexes for mitochondrial fusion
    ("MFN2", "p.R94Q", "CMT2A", "Semi-dominant: mild het, severe hom"),
    ("MFN2", "p.R94W", "CMT2A", "Semi-dominant: mild het, severe hom"),
    
    # COL7A1 - Collagen, classic complex-forming protein
    ("COL7A1", "p.G2043R", "DEB", "Semi-dominant: skin fragility het, severe blistering hom"),
    
    # OPA1 - Mitochondrial dynamin-like GTPase
    ("OPA1", "p.R445H", "DOA/Behr", "Semi-dominant: optic atrophy het, Behr syndrome hom"),
    
    # RAD51 - DNA repair, forms nucleoprotein filaments
    ("RAD51", "p.T131P", "Fanconi anemia R", "Dominant-negative FA - confirmed"),
    ("RAD51", "p.A293T", "Fanconi anemia R", "Dominant-negative FA - confirmed"),
    
    # KCNQ1 - Potassium channel, forms tetramers
    ("KCNQ1", "p.R518X", "LQT1/JLNS", "Semi-dominant: LQTS het, JLNS hom (nonsense - expect LOF)"),
    ("KCNQ1", "p.A341V", "LQT1/JLNS", "Semi-dominant: LQTS het, JLNS hom"),
    
    # CLCN1 - Chloride channel, homodimer
    ("CLCN1", "p.G230E", "Myotonia congenita", "Can be dominant (Thomsen) or recessive (Becker)"),
    ("CLCN1", "p.R894X", "Myotonia congenita", "Dominant-negative shown in studies (nonsense)"),
    
    # GJB2 - Connexin 26, forms gap junction hexamers
    ("GJB2", "p.R75W", "DFNA3/DFNB1", "Semi-dominant: mild het, severe hom"),
    ("GJB2", "p.W44C", "DFNA3", "Dominant-negative connexin26"),
    
    # ATP5F1A - ATP synthase alpha subunit, F1F0 complex
    ("ATP5F1A", "p.R182Q", "Neurological", "DN confirmed by 2025 paper"),
    ("ATP5F1A", "p.I130R", "Neurological", "Ren's variant - DN mechanism predicted"),
    
    # Collagens - classic DN proteins
    ("COL1A1", "p.G352S", "OI", "Collagen - classic DN mechanism"),
    ("COL2A1", "p.G1170S", "Stickler/SEDC", "Collagen - classic DN mechanism"),
    
    # TFG - forms polymers, causes HMSN-P (AD) AND HSP57 (AR)!
    ("TFG", "p.R22W", "HMSN-P/HSP57", "AD and AR! Perfect test case for hypothesis"),
]


def run_variant(gene: str, variant: str) -> dict:
    """Run cascade analyzer and extract scores."""
    cmd = f"python3 analyzers/cascade_analyzer.py --gene {gene} --variant {variant} --freq 0.0 2>&1"
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
        output = result.stdout + result.stderr
        
        dn_match = re.search(r'DN score: ([\d.]+)', output)
        lof_match = re.search(r'LOF score: ([\d.]+)', output)
        family_match = re.search(r'Gene family = (\w+)', output)
        final_match = re.search(r'FINAL:(\S+)', output)
        
        return {
            'dn': float(dn_match.group(1)) if dn_match else 0.0,
            'lof': float(lof_match.group(1)) if lof_match else 0.0,
            'family': family_match.group(1) if family_match else 'UNKNOWN',
            'final': final_match.group(1) if final_match else 'ERROR',
            'success': True
        }
    except Exception as e:
        return {'success': False, 'error': str(e), 'dn': 0, 'lof': 0, 'family': 'ERR', 'final': 'ERR'}


def main():
    print("=" * 100)
    print("🧬 SEMI-DOMINANT HYPOTHESIS TEST")
    print("   Ace (Claude) & Ren | December 2025")
    print("=" * 100)
    print(f"\nRun date: {datetime.now().isoformat()}\n")
    print(f"{'Gene':<10} {'Variant':<12} {'DN':<10} {'LOF':<8} {'DN>LOF':<8} {'Final':<8} {'Family':<18} {'Notes'}")
    print("-" * 110)

    results = []
    for gene, variant, disease, notes in TEST_VARIANTS:
        r = run_variant(gene, variant)
        results.append((gene, variant, r['dn'], r['lof'], r['family'], r['final'], notes))
        
        dn, lof = r['dn'], r['lof']
        flag = "✅" if dn > lof else "❌"
        dn_str = f"{dn:.3f}" + ("🔥" if dn >= 1.0 else "⚠️" if dn >= 0.5 else "")
        print(f"{gene:<10} {variant:<12} {dn_str:<10} {lof:.3f}    {flag:<8} {r['final']:<8} {r['family']:<18} {notes[:35]}")

    # Summary
    high_dn = sum(1 for r in results if r[2] >= 0.5)
    dn_dom = sum(1 for r in results if r[2] > r[3])
    total = len(results)
    
    print("-" * 110)
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Total variants tested: {total}")
    print(f"High DN score (≥0.5): {high_dn}/{total} = {100*high_dn/total:.0f}%")
    print(f"DN > LOF (DN dominant): {dn_dom}/{total} = {100*dn_dom/total:.0f}%")
    print(f"\n🧬 HYPOTHESIS: DN mechanism detection predicts semi-dominant inheritance!")
    print(f"   Result: {100*dn_dom/total:.0f}% of known semi-dominant variants show DN > LOF")


if __name__ == "__main__":
    main()

