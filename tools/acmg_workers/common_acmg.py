#!/usr/bin/env python3
"""
ACMG actionable genes batch runner (chunkable + threaded)

- Loads ACMG SF gene list from a file (one gene per line), or uses a small built-in
  starter list if none is provided.
- Uses the local ClinVar bulk VCF to extract variants for the chunk of genes.
- Builds discovery-format TSVs (gene, hgvs_p, hgvs_g, gnomad_freq) for Pathogenic/Likely Pathogenic variants.
- Runs the Cascade batch processor on each gene's TSV in parallel workers.

Run from repo root:
  python3 AdaptiveInterpreter/tools/acmg_workers/acmg_worker_1.py --genes-file /path/to/acmg_sf_genes.txt --outdir analysis/acmg_sf --clinvar-dir /mnt/Arcana/clinvar --max-workers 4

Notes:
- If --genes-file is omitted, a minimal fallback list is used. Provide the official ACMG SF list for full coverage.
- This script does not install anything; it uses local data and existing analyzers.
"""
from __future__ import annotations
import argparse
import json
import re
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Dict, Any, Tuple

# Ensure repo root is importable
REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

# Import extractor class directly
from AdaptiveInterpreter.utils.clinvar_bulk_extractor import ClinVarBulkExtractor

# ---------- Helpers ----------
ONE_LETTER = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G','His':'H',
    'Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V',
}

AA3_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")
AA1_RE = re.compile(r"p\.([A-Z])([0-9]+)([A-Z])")
G_RE   = re.compile(r":g\.")

def to_one_letter_p(hgvs_p: str) -> str:
    if AA1_RE.search(hgvs_p):
        return hgvs_p  # already one-letter
    m = AA3_RE.search(hgvs_p)
    if not m:
        return hgvs_p
    a1, pos, a2 = m.groups()
    return f"p.{ONE_LETTER.get(a1,a1[0])}{pos}{ONE_LETTER.get(a2,a2[0])}"

def pick_hgvs_p_and_g(hgvs_list: List[str]) -> Tuple[str|None, str|None]:
    hgvs_p, hgvs_g = None, None
    for h in hgvs_list:
        if 'p.' in h and hgvs_p is None:
            hgvs_p = h
        if ':g.' in h and hgvs_g is None:
            hgvs_g = h
        if hgvs_p and hgvs_g:
            break
    return hgvs_p, hgvs_g

def load_genes(genes_file: Path|None) -> List[str]:
    """Load genes from a simple text file (one per line) or from a Markdown file.
    - Text: one symbol per line (comments starting with # allowed)
    - Markdown: accepts fenced code blocks, tables, or bullet lists; extracts likely HGNC symbols
    """
    def _normalize_tokens(tokens: List[str]) -> List[str]:
        out, seen = [], set()
        for t in tokens:
            u = t.strip().upper()
            if not u:
                continue
            if u in seen:
                continue
            # Heuristic: typical HGNC symbols are 2–15 chars, start with a letter, allow digits and hyphens (e.g., NKX2-5)
            if not re.match(r"^[A-Z][A-Z0-9-]{1,14}$", u):
                continue
            seen.add(u)
            out.append(u)
        return out

    if genes_file and genes_file.exists():
        text = genes_file.read_text()
        # Fast path: treat as simple text list if it's not markdown-y
        if genes_file.suffix.lower() not in {'.md', '.markdown'}:
            genes = [ln.strip().split()[0] for ln in text.splitlines() if ln.strip() and not ln.strip().startswith('#')]
            return _normalize_tokens(genes)

        # Markdown parsing
        tokens: List[str] = []
        lines = text.splitlines()
        in_code = False
        for ln in lines:
            s = ln.strip()
            # fenced code block detection
            if s.startswith('```'):
                in_code = not in_code
                continue
            if in_code:
                # inside code block: assume one per line or whitespace separated
                tokens.extend(s.replace('|', ' ').split())
                continue
            # table rows starting with |
            if s.startswith('|') and not set(s) <= {'|', '-', ' '}:
                # take first cell as candidate
                parts = [p.strip() for p in s.strip('|').split('|')]
                if parts:
                    tokens.append(parts[0].split()[0])
                continue
            # bullet/numbered lists
            if s.startswith(('-', '*')) or re.match(r"^\d+\.\s+", s):
                tokens.extend(s.replace('-', ' ').replace('*', ' ').split())
                continue
            # bare gene symbols on their own line (like your file):
            if re.match(r"^[A-Za-z0-9_-]+$", s):
                tokens.append(s)
                continue
            # otherwise, ignore prose to avoid scooping acronyms
        normalized = _normalize_tokens(tokens)
        if normalized:
            return normalized
        # Fallback: treat as simple one-per-line text
        genes = [ln.strip().split()[0] for ln in text.splitlines() if ln.strip() and not ln.strip().startswith('#')]
        return _normalize_tokens(genes)

    # Minimal fallback starter (replace with official ACMG SF list for full run)
    return [
        'BRCA1','BRCA2','TP53','PTEN','STK11','MLH1','MSH2','MSH6','PMS2','APC',
        'LDLR','APOB','PCSK9','KCNQ1','KCNH2','SCN5A','RYR2','MYBPC3','MYH7','LMNA',
    ]

def chunk_list(items: List[str], num_chunks: int, index: int) -> List[str]:
    return [g for i, g in enumerate(items) if i % num_chunks == (index % num_chunks)]

# ---------- Core runner ----------

def run_acmg_chunk(chunk_index: int, num_chunks: int, genes_file: str|None, outdir: str,
                    clinvar_dir: str, max_workers: int = 4) -> None:
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    genes_all = load_genes(Path(genes_file) if genes_file else None)
    genes = chunk_list(genes_all, num_chunks, chunk_index)
    if not genes:
        print(f"No genes in chunk {chunk_index}/{num_chunks}.")
        return

    print(f"Chunk {chunk_index}/{num_chunks}: {len(genes)} genes → {genes}")

    # Extract variants for all genes in this chunk with a single VCF pass
    extractor = ClinVarBulkExtractor(clinvar_dir)
    gene_variants_map = extractor.extract_gene_variants(set(genes))  # {gene: [variant dicts]}

    # Build per-gene discovery TSVs (Pathogenic/Likely_pathogenic are mapped to 'pathogenic')
    discovery_paths: Dict[str, Path] = {}
    for gene, variants in gene_variants_map.items():
        rows: List[Dict[str, Any]] = []
        for var in variants:
            # Include ALL significance categories (pathogenic, likely_pathogenic, benign, LB, VUS, conflicting, etc.)
            hgvs_list = var.get('hgvs_list') or []
            hgvs_p, hgvs_g = pick_hgvs_p_and_g(hgvs_list)
            # Keep rows if we have either protein or genomic HGVS; skip only if we have neither
            if not hgvs_p and not hgvs_g:
                continue
            rows.append({
                'gene': gene,
                'hgvs_p': to_one_letter_p(hgvs_p) if hgvs_p else '',
                'hgvs_g': hgvs_g or '',
                'gnomad_freq': var.get('frequency') or 0.0,
                'clinical_sig': var.get('significance_raw',''),
                'molecular_consequence': var.get('molecular_consequence',''),
                'review_status': var.get('review_status',''),
            })
        if not rows:
            continue
        disc_path = outdir_p / f"{gene}.clinvarP.discovery.tsv"
        with disc_path.open('w') as f:
            header = ['gene','hgvs_p','hgvs_g','gnomad_freq','clinical_sig','molecular_consequence','review_status']
            f.write('\t'.join(header)+'\n')
            for r in rows:
                f.write('\t'.join([str(r.get(k,'')) for k in header])+'\n')
        discovery_paths[gene] = disc_path

    if not discovery_paths:
        print("No discovery TSVs built for this chunk (no P/LP variants found?).")
        return

    # Run cascade on each gene TSV in parallel
    def run_cascade(g: str, in_path: Path) -> Tuple[str, int]:
        out_path = outdir_p / f"{g}.clinvarP.cascade.tsv"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(REPO_ROOT/ 'AdaptiveInterpreter' / 'analyzers' / 'cascade_batch_processor.py'),
               '--input', str(in_path), '--output', str(out_path)]
        try:
            cp = subprocess.run(cmd, cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=3600)
            if cp.returncode != 0:
                print(f"[FAIL] {g}: {cp.stderr.strip()}")
            else:
                print(f"[OK]   {g}: → {out_path}")
            return g, cp.returncode
        except Exception as e:
            print(f"[EXC]  {g}: {e}")
            return g, 99

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(run_cascade, g, p) for g, p in discovery_paths.items()]
        for fut in as_completed(futures):
            fut.result()

# ---------- CLI ----------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description='Run ACMG actionable genes (chunked) through ClinVar→Cascade pipeline')
    p.add_argument('--genes-file', default=None, help='Text file of gene symbols (one per line). If omitted, uses a small fallback set.')
    p.add_argument('--outdir', default='analysis/acmg_sf', help='Output directory for discovery + cascade results')
    p.add_argument('--clinvar-dir', default='/mnt/Arcana/clinvar', help='Directory with ClinVar bulk files (expects clinvar.vcf.gz)')
    p.add_argument('--max-workers', type=int, default=4, help='Max parallel cascade runs per worker')
    p.add_argument('--num-chunks', type=int, default=4, help='Total number of chunks')
    p.add_argument('--chunk-index', type=int, default=0, help='This worker\'s chunk index (0-based)')
    return p

if __name__ == '__main__':
    ap = build_arg_parser()
    args = ap.parse_args()
    run_acmg_chunk(args.chunk_index, args.num_chunks, args.genes_file, args.outdir, args.clinvar_dir, args.max_workers)

