#!/usr/bin/env python3
import json
from DNModeling.cascade.cascade_analyzer import CascadeAnalyzer

CASES = [
    ("TP53",  "p.R248W"),
    ("TP53",  "p.D208V"),
    ("KIT",   "p.D816Y"),
    ("KIT",   "p.V559D"),
    ("SCN2A", "p.G1634V"),
]


def main():
    ca = CascadeAnalyzer()
    out = []
    for gene, var in CASES:
        print(f"\n=== {gene} {var} ===")
        try:
            res = ca.analyze_cascade_biological(gene, var, gnomad_freq=0.0, variant_type='missense', sequence=None)
            summary = {
                'gene': gene,
                'variant': var,
                'family': res.get('gene_family'),
                'family_aa_multiplier_applied': res.get('family_aa_multiplier_applied', 1.0),
                'DN': res['scores'].get('DN', 0.0),
                'LOF': res['scores'].get('LOF', 0.0),
                'GOF': res['scores'].get('GOF', 0.0),
                'final_score': res.get('final_score'),
                'final_classification': res.get('final_classification'),
                'summary': res.get('summary'),
            }
            out.append(summary)
            print(json.dumps(summary, indent=2))
        except Exception as e:
            print(f"ERROR {gene} {var}: {e}")

    print("\nCOMPACT JSON:")
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()

