"""Clinical relevance enrichment for CASCADE results.

Adds a clinical context layer on top of CASCADE's per-variant pathogenicity scores.
The layer answers: "given this variant is flagged as pathogenic-looking, is the gene
it's in actually clinically meaningful, or is this a known-noise gene family?"

This module is a post-processor by design: CASCADE scores mutations, the enricher
adds clinical context. Keeps the two concerns separate. Different data sources
(ClinVar / HPO / OMIM), different update cadences, different audiences.

Typical pipeline:
    cascade_results.tsv  →  clinical_enricher.enrich()  →  enriched_results.tsv

Author: Ace (Claude 4.x) + Ren, 2026-04-21
"""

from __future__ import annotations
import csv
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# --- constants ---

DEFAULT_CLINVAR_INDEX = Path('/home/Ace/clinvar_gene_index.json')
DEFAULT_NOISE_PATTERNS = Path(__file__).parent / 'clinical_noise_patterns.json'

# Tag rank — higher means more relevant. Used to sort variants within a
# classification bucket so 🚨 known-disease-gene variants rise above 🧂 noise.
TAG_RANK = {
    '🚨 known disease gene': 3,
    '⚠️ possible concern':    2,
    '❓ unknown':             1,
    '🧂 likely noise':        0,
}


class ClinicalEnricher:
    """Tags variants with clinical relevance using local ClinVar + noise patterns.

    Designed to be stateless-per-enrich-call but cache the indexes at construction.
    """

    def __init__(
        self,
        clinvar_index_path: Optional[Path] = None,
        noise_patterns_path: Optional[Path] = None,
    ):
        clinvar_path = Path(clinvar_index_path or DEFAULT_CLINVAR_INDEX)
        noise_path   = Path(noise_patterns_path or DEFAULT_NOISE_PATTERNS)

        if not clinvar_path.exists():
            raise FileNotFoundError(
                f'ClinVar gene index not found at {clinvar_path}. '
                f'Build it with scripts/build_clinvar_index.py.'
            )
        self.clinvar: Dict[str, Dict[str, Any]] = json.load(open(clinvar_path))

        self.noise_patterns: List[Tuple[re.Pattern, str]] = []
        if noise_path.exists():
            cfg = json.load(open(noise_path))
            for entry in cfg.get('patterns', []):
                self.noise_patterns.append(
                    (re.compile(entry['pattern']), entry['label'])
                )

    def classify_gene(self, gene: str) -> Dict[str, Any]:
        """Return the clinical-relevance classification for a gene symbol.

        Returns dict with keys: tag, tag_rank, reason, clinvar_p_lp_count,
        clinvar_total, top_diseases.
        """
        cv = self.clinvar.get(gene, {})
        p_lp = cv.get('P_submissions', 0) + cv.get('LP_submissions', 0)
        total = cv.get('total_submissions', 0)
        top_diseases = cv.get('top_diseases', [])
        top_name = top_diseases[0][0] if top_diseases else ''

        # Tier 1: explicit ClinVar P/LP presence overrides everything else
        if p_lp > 0:
            tag = '🚨 known disease gene'
            reason = f'{p_lp} P/LP submissions in ClinVar' + (
                f': {top_name}' if top_name else ''
            )
            return self._result(tag, reason, p_lp, total, top_diseases)

        # Tier 2: noise family pattern match, with noise-unless-busy escape hatch
        for pat, label in self.noise_patterns:
            if pat.match(gene):
                if total < 5:
                    return self._result(
                        '🧂 likely noise', label, p_lp, total, top_diseases
                    )
                return self._result(
                    '⚠️ possible concern',
                    f'{label} — but has {total} ClinVar entries',
                    p_lp, total, top_diseases,
                )

        # Tier 3: non-P/LP but ClinVar-known
        if total >= 10:
            return self._result(
                '⚠️ possible concern',
                f'{total} ClinVar entries, no P/LP submissions yet',
                p_lp, total, top_diseases,
            )
        if total > 0:
            return self._result(
                '⚠️ possible concern', f'{total} ClinVar entries',
                p_lp, total, top_diseases,
            )

        # Tier 4: not in ClinVar at all
        return self._result('❓ unknown', 'not in ClinVar', 0, 0, [])

    @staticmethod
    def _result(tag: str, reason: str, p_lp: int, total: int,
                top_diseases: List) -> Dict[str, Any]:
        return {
            'clinical_tag': tag,
            'clinical_tag_rank': TAG_RANK[tag],
            'clinical_reason': reason,
            'clinvar_p_lp_count': p_lp,
            'clinvar_total': total,
            'clinvar_top_diseases': '; '.join(
                f'{d}({n})' for d, n in (top_diseases or [])[:3]
            ),
        }

    def enrich_row(self, row: Dict[str, str]) -> Dict[str, str]:
        """Add clinical-relevance fields to a single CASCADE result row in-place.
        Returns the same dict for chaining."""
        gene = (row.get('gene') or '').strip()
        tags = self.classify_gene(gene)
        for k, v in tags.items():
            row[k] = str(v) if not isinstance(v, str) else v
        return row

    def enrich_tsv(self, input_path: Path, output_path: Path) -> int:
        """Read a CASCADE results TSV, write an enriched TSV. Returns row count."""
        input_path = Path(input_path)
        output_path = Path(output_path)

        with open(input_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            fieldnames = list(reader.fieldnames or [])
            # Add our new columns if not already present
            for new_col in ('clinical_tag', 'clinical_tag_rank', 'clinical_reason',
                            'clinvar_p_lp_count', 'clinvar_total',
                            'clinvar_top_diseases'):
                if new_col not in fieldnames:
                    fieldnames.append(new_col)

            rows = [self.enrich_row(dict(r)) for r in reader]

        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                    extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
        return len(rows)


# --- module-level convenience function ---

def enrich(input_tsv: str, output_tsv: str,
           clinvar_index: Optional[str] = None,
           noise_patterns: Optional[str] = None) -> int:
    """One-liner: enrich a CASCADE results TSV with clinical-relevance tags."""
    enricher = ClinicalEnricher(
        clinvar_index_path=clinvar_index,
        noise_patterns_path=noise_patterns,
    )
    return enricher.enrich_tsv(input_tsv, output_tsv)


# --- CLI ---

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(
        description='Enrich a CASCADE results TSV with clinical-relevance tags.')
    ap.add_argument('--input', required=True,
                    help='CASCADE results TSV to enrich')
    ap.add_argument('--output', required=True,
                    help='Output path for enriched TSV')
    ap.add_argument('--clinvar-index',
                    help='Override path to ClinVar gene index JSON')
    ap.add_argument('--noise-patterns',
                    help='Override path to noise patterns JSON config')
    args = ap.parse_args()

    n = enrich(args.input, args.output,
               clinvar_index=args.clinvar_index,
               noise_patterns=args.noise_patterns)
    print(f'✅ Enriched {n} rows → {args.output}')
