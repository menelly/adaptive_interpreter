#!/usr/bin/env python3
import json
from DNModeling.cascade.cascade_analyzer import CascadeAnalyzer


def main():
    # Force routing family to a dummy test family with a coefficients JSON
    ca = CascadeAnalyzer(override_family='test_family')
    gene = 'TP53'
    variant = 'p.R273H'
    # Build a fake sequence with R at position 273
    seq = ('A' * 272) + 'R' + ('A' * 20)

    result = ca.analyze_cascade_biological(gene, variant, gnomad_freq=0.0, variant_type='missense', sequence=seq)
    print("family:", result.get('gene_family'))
    print("family_aa_multiplier_applied:", result.get('family_aa_multiplier_applied'))
    print("final_score:", result.get('final_score'))
    print("summary:", result.get('summary'))

    # Emit compact JSON for CI visibility
    print("JSON:")
    print(json.dumps({
        'family': result.get('gene_family'),
        'family_aa_multiplier_applied': result.get('family_aa_multiplier_applied'),
        'final_score': result.get('final_score'),
    }, indent=2))


if __name__ == '__main__':
    main()

