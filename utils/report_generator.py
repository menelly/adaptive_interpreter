"""Markdown report generator for CASCADE (+clinical-enriched) results.

Reads a CASCADE results TSV (optionally already clinical-enriched) and produces
a human-readable markdown report with:
  - Counts by classification
  - Per-class tables sorted by severity then clinical-relevance rank
  - Associated disease names per variant (from ClinVar index)
  - MT-* filtering (upstream mitochondrial alignment is notoriously unreliable)
  - Optional full variant list with consequence-based overrides for non-missense
    (CASCADE was designed for missense; its calls on frameshift/splice/stop-gain
    are unreliable and get replaced by sensible defaults)

Pipeline shape:
    cascade_batch_processor → clinical_enricher → report_generator → Markdown

Author: Ace (Claude 4.x) + Ren, 2026-04-21
"""

from __future__ import annotations
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional

from .clinical_enricher import ClinicalEnricher


# --- constants ---

SEVERITY_ORDER = {
    'Pathogenic': 6, 'P': 6,
    'Likely Pathogenic': 5, 'LP': 5,
    'VUS-P': 4, 'VUS_P': 4,
    'VUS': 3,
    'VUS-B': 2,
    'LB': 1, 'Likely Benign': 1,
    'B': 0, 'Benign': 0,
}

CLASS_ICONS = {
    'Pathogenic': '🔴', 'P': '🔴',
    'Likely Pathogenic': '🟠', 'LP': '🟠',
    'VUS-P': '🟡', 'VUS': '⚪', 'VUS-B': '🟢',
    'LB': '🟢', 'B': '🟢', 'Benign': '🟢', 'Likely Benign': '🟢',
}

CLASS_DISPLAY_ORDER = [
    'Pathogenic', 'P', 'Likely Pathogenic', 'LP',
    'VUS-P', 'VUS', 'VUS-B', 'LB', 'Likely Benign', 'Benign', 'B',
    'Unknown',
]

AA_THREE = ('Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|'
            'Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val')
RE_MISSENSE_3LETTER = re.compile(
    rf'^p\.({AA_THREE})\d+({AA_THREE})$'
)
RE_MISSENSE_1LETTER = re.compile(r'^p\.[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]$')
RE_DISEASE_WITH_COUNT = re.compile(r'^(.+?)\((\d+)\)$')


# --- helpers ---

def _sev(classification: str) -> int:
    return SEVERITY_ORDER.get((classification or '').strip(), -1)


def _is_missense_variant(variant: str) -> bool:
    v = (variant or '').strip()
    if not v or 'fs' in v or v.endswith('*') or 'Ter' in v:
        return False
    return bool(RE_MISSENSE_3LETTER.match(v) or RE_MISSENSE_1LETTER.match(v))


def _is_mt_gene(row: Dict[str, Any]) -> bool:
    return (row.get('gene') or '').upper().startswith('MT-')


def _format_diseases(raw: str, top_n: int = 2) -> str:
    """'Disease1(204); Disease2(77); Disease3(26)' → 'Disease1 (×204), Disease2 (×77), +1 more'"""
    if not raw:
        return ''
    parts = [p.strip() for p in raw.split(';') if p.strip()]
    out = []
    for p in parts[:top_n]:
        m = RE_DISEASE_WITH_COUNT.match(p)
        if m:
            out.append(f'{m.group(1).strip()} (×{m.group(2)})')
        else:
            out.append(p)
    if len(parts) > top_n:
        out.append(f'+{len(parts) - top_n} more')
    return ', '.join(out)


def _escape(s: str) -> str:
    return (s or '').replace('|', '\\|').replace('\n', ' ')


def _consequence_override(row: Dict[str, Any], tag: str) -> Optional[tuple]:
    """Return (new_classification, reason) for non-missense variants.
    CASCADE was designed for missense — its calls on frameshift/splice/stop-gain
    are unreliable. Override with consequence-based defaults.
    Returns None if no override applies."""
    consequence = (row.get('molecular_consequence') or '').strip().lower()
    location = (row.get('review_status') or '').strip().lower()
    variant = (row.get('variant') or '').strip()

    is_frameshift = 'frameshift' in consequence or 'fs' in variant
    is_stop_gain = ('stop gain' in consequence or variant.rstrip(')').endswith('*')
                    or 'Ter' in variant)
    is_splice = 'splice' in location or 'splice' in consequence
    is_utr = 'utr' in location or 'utr' in consequence
    known_disease = '🚨' in tag

    if is_frameshift or is_stop_gain:
        return ('LP' if known_disease else 'VUS-P',
                'frameshift/stop — consequence-based default')
    if is_splice:
        return ('LP' if known_disease else 'VUS-P',
                'splice variant — consequence-based default')
    if is_utr:
        return ('VUS-P' if known_disease else 'VUS',
                'UTR variant — default conservative')
    return None


# --- the generator ---

class ReportGenerator:
    """Turn an enriched CASCADE results TSV into a clinician-friendly markdown report."""

    def __init__(self, enricher: Optional[ClinicalEnricher] = None):
        # If the input TSV isn't enriched, we'll enrich on the fly.
        self.enricher = enricher or ClinicalEnricher()

    def generate(
        self,
        input_tsv: str,
        output_md: str,
        person_name: str = 'Unnamed',
        full_list: bool = False,
        exclude_mt: bool = True,
    ) -> int:
        """Read enriched (or un-enriched) TSV → write markdown report.

        Args:
            input_tsv:    path to CASCADE results (enriched or raw)
            output_md:    where to write the markdown report
            person_name:  name used in report title (e.g. "Winnie")
            full_list:    if True, include non-missense variants with
                          consequence-based overrides; if False (default),
                          missense only
            exclude_mt:   if True (default), drop MT-* genes (Dante MT
                          alignment issues make these unreliable)

        Returns: number of variants in the output report.
        """
        rows = list(csv.DictReader(open(input_tsv), delimiter='\t'))

        # Auto-enrich if needed
        if rows and 'clinical_tag' not in rows[0]:
            rows = [self.enricher.enrich_row(dict(r)) for r in rows]

        # Filter
        if exclude_mt:
            rows = [r for r in rows if not _is_mt_gene(r)]
        mt_count = 0  # we could track this but it's lost post-filter; skip

        missense_count = sum(1 for r in rows if _is_missense_variant(r.get('variant', '')))
        if full_list:
            kept = rows
        else:
            kept = [r for r in rows if _is_missense_variant(r.get('variant', ''))]

        nonmissense_count = len(kept) - sum(1 for r in kept if _is_missense_variant(r.get('variant', '')))

        # Apply consequence overrides for non-missense in full-list mode
        for r in kept:
            r['_is_missense'] = _is_missense_variant(r.get('variant', ''))
            if full_list and not r['_is_missense']:
                override = _consequence_override(r, r.get('clinical_tag', ''))
                if override:
                    r['_orig_cascade'] = r.get('adj_classification', '')
                    r['adj_classification'] = override[0]
                    r['_override_reason'] = override[1]
            r['_sev'] = _sev(r.get('adj_classification'))
            r['_tag_rank'] = int(r.get('clinical_tag_rank', 0) or 0)

        kept.sort(key=lambda r: (r['_sev'], r['_tag_rank'],
                                 float(r.get('adj_score') or 0)),
                  reverse=True)

        # Group by classification
        buckets: Dict[str, List[Dict]] = defaultdict(list)
        for r in kept:
            buckets[(r.get('adj_classification') or 'Unknown').strip()].append(r)

        ordered = ([k for k in CLASS_DISPLAY_ORDER if k in buckets]
                   + [k for k in buckets if k not in CLASS_DISPLAY_ORDER])

        # Emit markdown
        lines = [f'# 🧬 {person_name} — CASCADE + Clinical Relevance', '']
        if full_list:
            missense_in_kept = sum(1 for r in kept if r['_is_missense'])
            lines.append(
                f'**{len(kept)} variants** '
                f'(MT-* {"excluded" if exclude_mt else "included"}). '
                f'Missense: {missense_in_kept}. '
                f'Non-missense (consequence-based defaults): {len(kept) - missense_in_kept}.'
            )
        else:
            lines.append(
                f'**{len(kept)} missense variants** '
                f'(MT-* {"excluded" if exclude_mt else "included"}).'
            )

        lines += [
            '',
            'Clinical relevance tags (from local ClinVar index):',
            '- **🚨 known disease gene** — ClinVar has ≥1 P/LP submission for this gene',
            '- **⚠️ possible concern** — gene has ClinVar entries but no established P/LP, '
            'or is in a mostly-noise family with elevated submission count',
            '- **🧂 likely noise** — olfactory / HLA / keratin-associated / Ig V-region / '
            'etc., no ClinVar P/LP',
            '- **❓ unknown** — not in ClinVar',
            '',
        ]

        if full_list:
            lines += [
                'Frameshift / stop-gain / splice / UTR variants use consequence-based defaults '
                'because CASCADE was designed for missense and its calls on those classes are '
                'unreliable. Defaults applied:',
                '- Frameshift/Stop-gain in 🚨 gene → **LP** (LOF actionable)',
                '- Frameshift/Stop-gain elsewhere → **VUS-P**',
                '- Splice donor/acceptor in 🚨 → **LP**',
                '- Splice donor/acceptor elsewhere → **VUS-P**',
                '- UTR in 🚨 → **VUS-P**; UTR elsewhere → **VUS**',
                '',
            ]

        # Counts table
        lines += ['## Counts by classification', '', '| Class | Count |', '|---|---|']
        for k in ordered:
            lines.append(f'| {k} | {len(buckets[k])} |')
        lines.append('')

        # Per-class tables
        for k in ordered:
            icon = CLASS_ICONS.get(k, '❓')
            lines.append(f'## {icon} {k}  ({len(buckets[k])})')
            lines.append('')
            if full_list:
                lines.append('| Clinical | Type | Gene | Variant | Score | '
                             'Associated Disease(s) | Notes |')
                lines.append('|---|---|---|---|---|---|---|')
            else:
                lines.append('| Clinical | Gene | Variant | Score | Dom | '
                             'Associated Disease(s) | Reasoning |')
                lines.append('|---|---|---|---|---|---|---|')
            for r in buckets[k]:
                diseases = _format_diseases(r.get('clinvar_top_diseases', ''))
                if full_list:
                    if r['_is_missense']:
                        notes = _escape(r.get('explanation', ''))[:100]
                    else:
                        orig = r.get('_orig_cascade', '')
                        ov = r.get('_override_reason', '')
                        suffix = f' (CASCADE→{orig})' if orig else ''
                        notes = f'{ov}{suffix}'[:100]
                    cells = [
                        r.get('clinical_tag', ''),
                        'missense' if r['_is_missense'] else 'non-missense',
                        _escape(r.get('gene', '')),
                        _escape(r.get('variant', '')),
                        (r.get('adj_score') or '')[:6],
                        _escape(diseases)[:80],
                        _escape(notes),
                    ]
                else:
                    reason = _escape(r.get('explanation', ''))[:120]
                    cells = [
                        r.get('clinical_tag', ''),
                        _escape(r.get('gene', '')),
                        _escape(r.get('variant', '')),
                        (r.get('adj_score') or '')[:6],
                        _escape((r.get('gene_family') or '')[:12]),
                        _escape(diseases)[:80],
                        reason,
                    ]
                lines.append('| ' + ' | '.join(cells) + ' |')
            lines.append('')

        Path(output_md).write_text('\n'.join(lines))
        return len(kept)


# --- convenience function ---

def generate_report(
    input_tsv: str,
    output_md: str,
    person_name: str = 'Unnamed',
    full_list: bool = False,
    exclude_mt: bool = True,
) -> int:
    """One-liner: enriched TSV → markdown report."""
    gen = ReportGenerator()
    return gen.generate(input_tsv, output_md, person_name=person_name,
                        full_list=full_list, exclude_mt=exclude_mt)


# --- CLI ---

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(
        description='Generate a markdown clinical report from a CASCADE results TSV.')
    ap.add_argument('--input', required=True,
                    help='Input CASCADE results TSV (will be auto-enriched if needed)')
    ap.add_argument('--output', required=True,
                    help='Output markdown path')
    ap.add_argument('--person', default='Unnamed',
                    help='Name for the report title')
    ap.add_argument('--full-list', action='store_true',
                    help='Include non-missense variants with consequence-based '
                         'overrides (default: missense only)')
    ap.add_argument('--include-mt', action='store_true',
                    help='Include mitochondrial MT-* genes (default: excluded, '
                         'because upstream MT alignment is unreliable)')
    args = ap.parse_args()

    n = generate_report(
        args.input, args.output,
        person_name=args.person,
        full_list=args.full_list,
        exclude_mt=not args.include_mt,
    )
    print(f'✅ Wrote report with {n} variants → {args.output}')
