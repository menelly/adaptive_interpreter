#!/usr/bin/env python3
"""
🌊 CASCADE BATCH PROCESSOR - Multi-Analyzer CSV Processing
Process CSV files through the complete DN → LOF → GOF cascade system!
"""
import argparse
import csv
import json
import os
import sys
from typing import Dict, List
from pathlib import Path
import re

# This is necessary because the script is run from a different directory
# and needs to find the AdaptiveInterpreter package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from AdaptiveInterpreter.utils.gnomad_frequency_fetcher import GnomADFrequencyFetcher
from AdaptiveInterpreter.utils.conservation_fetcher import ConservationFetcher
from AdaptiveInterpreter.analyzers.cascade_analyzer import CascadeAnalyzer
from AdaptiveInterpreter.nova_dn.csv_batch_processor import CSVBatchProcessor
from AdaptiveInterpreter.nova_dn.mixed_mechanism_resolver import UnifiedMechanismResolver
from AdaptiveInterpreter.utils.genomic_to_protein import GenomicToProteinConverter
from AdaptiveInterpreter import config

_CHROM_OK_RE = re.compile(r"^(chr)?(\d+|[XYxy]|[Mm][Tt]?)$")


def _is_real_chrom(chr_val: str) -> bool:
    """True iff chr_val looks like a real chromosome identifier.

    Accepts 'chr7', '7', 'chrX', 'X', 'chrMT', 'MT' etc. Rejects the
    Franklin metadata-footer junk values ('scv', 'Taya', 'High, Medium,
    Low', blank, etc.).
    """
    return bool(_CHROM_OK_RE.match((chr_val or "").strip()))


class CascadeBatchProcessor:
    def __init__(self, alphafold_path: str = "/mnt/Arcana/alphafold_human/structures/",
                 override_family: str = None, conservative_mode: bool = False,
                 ignore_quality_floor: bool = False):
        self.ignore_quality_floor = ignore_quality_floor
        self.cascade_analyzer = CascadeAnalyzer(alphafold_path, override_family, conservative_mode)
        self.csv_processor = CSVBatchProcessor(alphafold_path)
        self.unified_resolver = UnifiedMechanismResolver()
        self.frequency_fetcher = GnomADFrequencyFetcher()
        self.frequency_fetcher.load_cache()
        try:
            # Use absolute path from config
            conservation_path = config.CONSERVATION_DATA_PATH / "hg38.phyloP100way.bw"
            self.conservation_fetcher = ConservationFetcher(str(conservation_path))
            print(f"✅ Conservation data loaded from: {conservation_path}")
        except Exception as e:
            print(f"⚠️ WARNING: Could not initialize ConservationFetcher: {e}")
            print(f"⚠️ Variants without conservation data will be classified as VUS for safety")
            self.conservation_fetcher = None
        # Dynamic gene-to-chromosome cache - fetches from Ensembl on demand
        self._gene_chromosome_cache = {}
        # Per-Ren: autosomal recessive genes can have high population AF even for pathogenic variants.
        # Use a higher benign short-circuit threshold for AR genes (default 5%).
        self.autosomal_recessive_genes = {
            'HFE', 'BTD', 'CFTR', 'FKRP'
        }
        self.ar_benign_threshold = 0.05  # 5% minimum for AR genes

    def get_gene_chromosome(self, gene: str) -> str:
        """Dynamically fetch chromosome for gene from Ensembl REST API, with caching."""
        if gene in self._gene_chromosome_cache:
            return self._gene_chromosome_cache[gene]

        try:
            import requests
            url = f'https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json'
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                data = r.json()
                chrom = f"chr{data.get('seq_region_name', '')}"
                if chrom and chrom != 'chr':
                    self._gene_chromosome_cache[gene] = chrom
                    print(f"📍 Fetched chromosome for {gene}: {chrom}")
                    return chrom
        except Exception as e:
            print(f"⚠️ Could not fetch chromosome for {gene}: {e}")

        return None

    def prepare_genes(self, input_path: str) -> set:
        """
        🧬 PREP STEP: Scan input for unique genes and pre-fetch all per-gene data.

        This runs ONCE before variant analysis so we don't re-fetch the same
        InterPro domains, UniProt annotations, and GO terms for every variant
        in the same gene. Ren's idea — stop doing redundant API calls!
        """
        print("\n🔧 PREP STEP: Scanning input for unique genes...")
        genes = set()
        try:
            # Auto-detect CSV/TSV and the column name housing the gene symbol.
            # Discovery/Ren formats use 'gene'; Franklin uses 'Gene' (and the
            # file may be CSV).
            with open(input_path, 'r') as f:
                first_line = f.readline()
            delim = '\t' if '\t' in first_line else ','
            with open(input_path, 'r') as f:
                reader = csv.DictReader(f, delimiter=delim)
                fieldnames = reader.fieldnames or []
                gene_col = 'gene' if 'gene' in fieldnames else ('Gene' if 'Gene' in fieldnames else None)
                if gene_col is None:
                    return genes
                for row in reader:
                    gene = (row.get(gene_col) or '').strip()
                    if not gene:
                        continue
                    # Skip mtDNA — unreliable in WGS due to reference drift.
                    if gene.startswith('MT-'):
                        continue
                    # Skip Franklin's metadata footer rows where 'Chr' is blank
                    # or not a real chromosome (e.g. "scv" on a "Filter logic"
                    # row, or "Taya" on an "Analysis Name" row), and mtDNA
                    # by chromosome. Accept both 'chr7' (Dante-source Franklin)
                    # and '7' (Nebula-source Franklin) styles.
                    chr_val = (row.get('Chr') or '').strip()
                    if 'Chr' in fieldnames:
                        if not _is_real_chrom(chr_val):
                            continue
                        if chr_val.lower() in ('chrm', 'chrmt', 'm', 'mt'):
                            continue
                    genes.add(gene)
        except Exception as e:
            print(f"⚠️ Could not scan input for genes: {e}")
            return genes

        print(f"   Found {len(genes)} unique genes: {', '.join(sorted(genes))}")

        # 1. Resolve UniProt IDs for all genes
        from AdaptiveInterpreter.data_processing.universal_protein_annotator import UniversalProteinAnnotator
        annotator = UniversalProteinAnnotator()
        gene_uniprot = {}
        for gene in sorted(genes):
            uid = self.cascade_analyzer.gene_to_uniprot.get(gene)
            if not uid:
                try:
                    uid = annotator._find_uniprot_id(gene)
                    if uid:
                        self.cascade_analyzer.gene_to_uniprot[gene] = uid
                except Exception:
                    pass
            if uid:
                gene_uniprot[gene] = uid
                print(f"   ✅ {gene} → {uid}")
            else:
                print(f"   ⚠️ {gene} → no UniProt ID found")

        # 2. Cache InterPro domains for all genes (fetch once, use forever)
        from AdaptiveInterpreter.data_processing.cache_interpro_domains import InterProDomainCacher
        cacher = InterProDomainCacher()
        print(f"\n🔗 Caching InterPro domains for {len(gene_uniprot)} genes...")
        for gene, uid in gene_uniprot.items():
            cacher.fetch_interpro_domains(uid)

        # 3. Pre-fetch protein annotations (function + GO terms) for gene family classification
        print(f"\n📋 Pre-caching protein annotations...")
        for gene, uid in gene_uniprot.items():
            cache_file = f"protein_annotations_cache/{uid}_domains.json"
            if os.path.exists(cache_file):
                print(f"   ⚡ {gene} annotations already cached")
            else:
                try:
                    annotator.annotate_protein(gene, uid)
                    print(f"   ✅ {gene} annotations fetched and cached")
                except Exception as e:
                    print(f"   ⚠️ {gene} annotation fetch failed: {e}")

        print(f"\n✅ PREP COMPLETE: {len(gene_uniprot)} genes ready for analysis")
        print("=" * 60)
        return genes

    def process_csv(self, input_path: str, output_path: str, freq_threshold: float = 0.05) -> Dict:
        results = []
        stats = self.get_initial_stats()
        print(f"🌊 Processing variants through CASCADE SYSTEM")
        print(f"📊 Frequency threshold: {freq_threshold}")
        print("=" * 60)

        # 🧬 PREP STEP: Pre-fetch all per-gene data before variant loop
        self.prepare_genes(input_path)

        try:
            # Auto-detect CSV vs TSV — Franklin exports either. Sniff from first line.
            with open(input_path, 'r') as f:
                first_line = f.readline()
            delim = '\t' if '\t' in first_line else ','
            with open(input_path, 'r') as f:
                reader = csv.DictReader(f, delimiter=delim)
                fieldnames = reader.fieldnames or []
                is_discovery_format = 'hgvs_g' in fieldnames and 'hgvs_p' in fieldnames
                is_ren_format = 'Chrpos' in fieldnames and 'AA Chg' in fieldnames
                is_franklin_format = 'AA Change' in fieldnames and 'Effect' in fieldnames

                for row in reader:
                    stats['total_variants'] += 1
                    try:
                        parsed_info = self.parse_variant_row(row, is_discovery_format, is_ren_format, freq_threshold, is_franklin_format)
                        if parsed_info.get('skipped'):
                            reason = parsed_info['skip_reason']
                            stats[reason] = stats.get(reason, 0) + 1
                            continue
                        # Universal mtDNA skip — mtDNA alignment in WGS is
                        # unreliable (aligned to a separate reference in most
                        # pipelines). Filter regardless of input format.
                        parsed_gene = (parsed_info.get('gene') or '')
                        parsed_chrom = (parsed_info.get('chrom') or '').lower()
                        if parsed_gene.startswith('MT-') or parsed_chrom in ('chrm', 'chrmt', 'm', 'mt'):
                            stats['skipped_mtdna'] = stats.get('skipped_mtdna', 0) + 1
                            continue

                        gene = parsed_info['gene']
                        variant = parsed_info['variant']
                        gnomad_freq = float(parsed_info.get('gnomad_freq') or 0.0)
                        hgvs = parsed_info['hgvs']
                        clinical_sig = parsed_info['clinical_sig']
                        conservation_score = parsed_info.get('conservation_score')

                        # Early benign short-circuit: use AR-aware dynamic threshold
                        effective_threshold = freq_threshold
                        if gene in self.autosomal_recessive_genes:
                            effective_threshold = max(freq_threshold, self.ar_benign_threshold)
                        if gnomad_freq >= effective_threshold:
                            print(f"\n🧬 Processing {gene} {variant} (freq: {gnomad_freq:.4f})")
                            print(f"🧹 High population frequency (≥ {effective_threshold:.4f}) — classifying as Benign and skipping heavy analysis.")
                            result = {
                                'gene': gene,
                                'variant': variant,
                                'status': 'SUCCESS',
                                'final_classification': 'LB',
                                'final_score': 0.0,
                                'explanation': f'High population frequency ({gnomad_freq:.4%}) ≥ threshold {effective_threshold:.2%} — benign by frequency',
                                'scores': {'DN': 0.0, 'LOF': 0.0, 'GOF': 0.0},
                                'classifications': {'DN': 'NOT_RUN', 'LOF': 'NOT_RUN', 'GOF': 'NOT_RUN'},
                                'analyzers_run': ['FREQUENCY_ONLY']
                            }
                            result['hgvs'] = hgvs
                            result['variant_type'] = parsed_info.get('variant_type')
                            result['gnomad_freq'] = gnomad_freq
                            result['clinical_sig'] = clinical_sig
                            result['review_status'] = parsed_info.get('review_status','')
                            result['molecular_consequence'] = parsed_info.get('molecular_consequence','')
                            result['rcv'] = parsed_info.get('rcv','')
                            result['variation_id'] = parsed_info.get('variation_id','')
                            result['zygosity'] = parsed_info.get('zygosity','')
                            result['confidence'] = parsed_info.get('confidence','')
                            results.append(result)
                            stats['processed'] += 1
                            continue

                        print(f"\n🧬 Processing {gene} {variant} (freq: {gnomad_freq:.4f})")

                        # Determine variant type from ClinVar if available; fallback to heuristic
                        variant_type = parsed_info.get('variant_type')
                        if not variant_type:
                            variant_type = 'splice_suspect' if (variant or '').strip().lower() in {'p.?', 'p.?*', 'p.?'} else 'missense'

                        result = self.cascade_analyzer.analyze_cascade_biological(
                            gene, variant, gnomad_freq, variant_type, sequence=None, conservation_score=conservation_score
                        )
                        # Preserve original inputs and ClinVar context for downstream consumers
                        result['hgvs'] = hgvs
                        result['variant_type'] = variant_type
                        result['gnomad_freq'] = gnomad_freq
                        result['clinical_sig'] = clinical_sig
                        result['review_status'] = parsed_info.get('review_status','')
                        result['molecular_consequence'] = parsed_info.get('molecular_consequence','')
                        result['rcv'] = parsed_info.get('rcv','')
                        result['variation_id'] = parsed_info.get('variation_id','')
                        result['zygosity'] = parsed_info.get('zygosity','')
                        result['confidence'] = parsed_info.get('confidence','')
                        results.append(result)
                        stats['processed'] += 1
                    except Exception as e:
                        stats['failed'] += 1
                        print(f"❌ Unhandled error in row {row}: {e}")
        except Exception as e:
            print(f"❌ Error reading CSV: {e}")
            return {'error': str(e)}
        
        self.finalize_processing(results, output_path, stats)
        return {'stats': stats, 'results_file': output_path}

    def parse_variant_row(self, row, is_discovery, is_ren, freq_threshold, is_franklin=False):
        if is_franklin:
            return self.parse_franklin_format(row)
        if is_discovery:
            return self.parse_discovery_format(row)
        if is_ren:
            return self.parse_ren_format(row)
        return {'skipped': True, 'skip_reason': 'skipped_unparseable'}

    def parse_franklin_format(self, row):
        """Parse Franklin/Genoox CSV/TSV export row (missense-only).

        Franklin exports both missense and non-missense in the same file.
        This parser scores ONLY missense — non-missense consequences
        (frameshift, stop-gain, splice, etc.) are handled downstream by
        cumbursum.franklin_adapter via consequence-override, which doesn't
        need AI scoring. The same source file can be fed to both paths.

        Key columns:
            Gene, Chr, Start Position, Ref, Alt, Transcript, AA Change,
            Effect, Zygosity, gnomAD (Exome), gnomAD (Genome),
            Aggregated Frequency, Genoox Classification, ClinVar Classification
        """
        gene = (row.get('Gene') or '').strip()
        effect = (row.get('Effect') or '').strip()
        aa_change = (row.get('AA Change') or '').strip()
        chr_val = (row.get('Chr') or '').strip()

        # Franklin appends a metadata footer (Analysis Name, User email, Date
        # of Export, ...) where Chr is not a real chromosome. Accept both
        # 'chr7' (Dante-source) and '7' (Nebula-source) styles; reject
        # everything that isn't clearly a chromosome identifier.
        if not _is_real_chrom(chr_val):
            return {'skipped': True, 'skip_reason': 'skipped_metadata_footer'}

        # Skip mtDNA — WGS pipelines (Dante et al) usually align mtDNA to a
        # different reference than nuclear DNA, so mtDNA calls are riddled
        # with alignment artifacts. Filter at both chr and gene-symbol level.
        if chr_val.lower() in ('chrm', 'chrmt', 'm', 'mt') or gene.startswith('MT-'):
            return {'skipped': True, 'skip_reason': 'skipped_mtdna'}

        # Only missense goes through AI scoring.
        # Franklin uses 'Missense' (Dante-source title case) OR 'NON_SYNONYMOUS'
        # (Nebula-source upper_snake) — both mean 'missense variant'.
        eff_norm = effect.lower().replace('_', ' ').strip()
        if eff_norm not in ('missense', 'non synonymous'):
            return {'skipped': True, 'skip_reason': 'skipped_non_missense'}

        if not gene or not aa_change:
            return {'skipped': True, 'skip_reason': 'skipped_unparseable'}

        # --- Prepass filter ------------------------------------------------
        # (a) Max population allele frequency > 10%. Floor is 10% not 5% to
        #     keep pathogenic-but-common variants like AMPD1 Q12* (8% AF,
        #     frequently tolerated in homs) and similar AR-tolerated alleles.
        # (b) Technical QUAL < 20 — below GATK's standard reliability floor.
        # Note: `Internal Sample Count` is NOT used — in Franklin it counts
        # the uploading user's own samples (family), so high values are an
        # inheritance signal, opposite of a benign signal.
        def _f(s):
            s = (s or '').strip()
            if not s or s.upper() in ('N/A', 'NA', 'NONE'):
                return 0.0
            try:
                return float(s)
            except ValueError:
                return 0.0

        freq_cols = ('gnomAD (Exome)', 'gnomAD (Genome)', 'ExAC (All)',
                     '1000 genomes', 'ESP 6500', 'UK10K TwinsUK',
                     'UK10K ALSPAC', 'Aggregated Frequency')
        max_pop_af = max((_f(row.get(c)) for c in freq_cols), default=0.0)
        if max_pop_af > 0.10:
            return {'skipped': True, 'skip_reason': 'skipped_common'}

        qual = _f(row.get('Quality'))
        if 0 < qual < 20 and not self.ignore_quality_floor:
            return {'skipped': True, 'skip_reason': 'skipped_low_quality'}
        # -------------------------------------------------------------------

        # Genomic coords
        pos = None
        try:
            pos_raw = (row.get('Start Position') or '').strip()
            pos = int(pos_raw) if pos_raw else None
        except ValueError:
            pos = None
        ref = (row.get('Ref') or '').strip()
        alt = (row.get('Alt') or '').strip()
        chrom = chr_val  # 'chr1', 'chrX', ...

        # Build a synthetic genomic HGVS (good enough for downstream regex parser)
        hgvs_g = None
        if pos and ref and alt and len(ref) == 1 and len(alt) == 1 and set(ref + alt) <= set("ACGT"):
            cn_raw = chrom.lower().replace('chr', '')
            cn_map = {'x': 23, 'y': 24, 'm': 25, 'mt': 25}
            cn = None
            if cn_raw.isdigit():
                cn = int(cn_raw)
            elif cn_raw in cn_map:
                cn = cn_map[cn_raw]
            if cn is not None:
                hgvs_g = f"NC_{cn:06d}.11:g.{pos}{ref}>{alt}"

        # gnomAD freq — prefer exome > genome > aggregated
        gnomad_freq = 0.0
        for col in ['gnomAD (Exome)', 'gnomAD (Genome)', 'Aggregated Frequency']:
            val = (row.get(col) or '').strip()
            if val and val.upper() not in ('N/A', 'NA', ''):
                try:
                    gnomad_freq = float(val)
                    break
                except ValueError:
                    continue

        # Conservation (if available)
        conservation_score = None
        if self.conservation_fetcher and pos:
            conservation_score = self.conservation_fetcher.get_conservation_score(chrom, pos)
            if conservation_score is not None:
                print(f"✅ Conservation score for {gene} {aa_change} ({chrom}:{pos}): {conservation_score:.4f}")

        return {
            'gene': gene,
            'variant': aa_change,
            'hgvs': hgvs_g or aa_change,
            'hgvs_g': hgvs_g,
            'chrom': chrom,
            'pos': pos,
            'ref': ref if len(ref) == 1 else None,
            'alt': alt if len(alt) == 1 else None,
            'gnomad_freq': gnomad_freq,
            'clinical_sig': (row.get('Genoox Classification') or '').strip()[:80],
            'clinvar_sig': (row.get('ClinVar Classification') or '').strip()[:80],
            'rsid': (row.get('dbSNP') or '').strip(),
            'variant_type': 'missense',
            'review_status': '',
            'molecular_consequence': effect,
            'conservation_score': conservation_score,
            # Franklin-native context so downstream consumers (cumbursum
            # adapter) can apply zygosity weighting and low-confidence
            # hom downgrades without re-reading the Franklin CSV.
            'zygosity': (row.get('Zygosity') or '').strip(),
            'confidence': (row.get('Confidence') or '').strip(),
        }

    def parse_ren_format(self, row):
        """Parse GeneCards/Ren TSV format with Variation Name and AA Chg columns"""
        variation_name = row.get('Variation Name', '')
        aa_chg = row.get('AA Chg', '')

        if not aa_chg or not aa_chg.strip():
            return {'skipped': True, 'skip_reason': 'skipped_no_protein_consequence'}

        # Extract gene from Variation Name like "NM_000088.4(COL1A1):c.2678G>C (p.Gly893Ala)"
        gene_match = re.search(r'\(([A-Z0-9]+)\):', variation_name)
        gene = gene_match.group(1) if gene_match else None

        if not gene:
            return {'skipped': True, 'skip_reason': 'skipped_unparseable'}

        # Convert AA Chg (like G893A) to variant format analyzer expects
        aa_chg = aa_chg.strip()

        # Build hgvs from variation name or construct from aa_chg
        hgvs = variation_name if variation_name else f"p.{aa_chg}"

        # Parse Chrpos like "50,189,528(-)" to get position
        chrpos_raw = row.get('Chrpos', '')
        pos = None
        if chrpos_raw:
            # Remove commas and extract just the number
            pos_match = re.search(r'([\d,]+)', chrpos_raw)
            if pos_match:
                pos = int(pos_match.group(1).replace(',', ''))

        # Get chromosome from gene dynamically via Ensembl API (cached)
        chrom = self.get_gene_chromosome(gene)

        # Fetch conservation score if we have position and chromosome
        conservation_score = None
        if self.conservation_fetcher and chrom and pos:
            conservation_score = self.conservation_fetcher.get_conservation_score(chrom, pos)
            if conservation_score is not None:
                print(f"✅ Conservation score for {gene} {aa_chg} ({chrom}:{pos}): {conservation_score:.4f}")

        return {
            'gene': gene,
            'variant': aa_chg,
            'hgvs': hgvs,
            'hgvs_g': None,
            'chrom': chrom,
            'pos': pos,
            'ref': None,
            'alt': None,
            'clinical_sig': row.get('Clinical significance and condition', '')[:80],
            'clinvar_sig': row.get('Clinical significance and condition', '')[:80],
            'rsid': row.get('rsID', ''),
            'variant_type': 'missense',
            'review_status': '',
            'molecular_consequence': row.get('Type', ''),
            'conservation_score': conservation_score
        }

    def parse_discovery_format(self, row):
        gene, variant_p, hgvs_g = row.get('gene'), row.get('hgvs_p'), row.get('hgvs_g')
        # Accept rows with either protein or genomic HGVS (or both)
        if not gene or (not variant_p and not hgvs_g):
            return {'skipped': True, 'skip_reason': 'skipped_unparseable'}

        # Determine if provided protein HGVS is parseable/useful; treat p.? and similar as uninformative
        protein_uninformative = False
        if variant_p:
            vp = str(variant_p).strip()
            vpl = vp.lower()
            if vpl in {"p.?", "p.?*", "p.", "p.=", "p.(=)"}:
                protein_uninformative = True
            else:
                # Heuristic: must contain a position to be useful (e.g., p.Arg273His, p.R273H, p.*577Qext*)
                if not re.search(r"p\.[A-Za-z\*][A-Za-z]{0,2}\d+", vp):
                    protein_uninformative = True

        # Use protein HGVS when available and informative; otherwise fall back to genomic string for reporting
        variant = (variant_p if (variant_p and not protein_uninformative) else hgvs_g)

        # 🚨 EARLY SKIP: skip ONLY if we have neither a usable protein HGVS
        # NOR a genomic HGVS. When hgvs_g is present (e.g. VCF-derived input
        # via vcf_to_discovery.py), the downstream genomic→protein converter
        # will try to resolve hgvs_p — don't short-circuit before that runs.
        if (not variant_p or protein_uninformative) and not hgvs_g:
            return {
                'skipped': True,
                'skip_reason': 'skipped_no_protein_consequence'
            }

        # Parse genomic coordinates and alleles from hgvs_g like: NC_000017.11:g.50201411T>C
        chrom_num = None
        pos = None
        ref = None
        alt = None
        if hgvs_g:
            m = re.search(r'NC_0*(\d+)\.\d+:g\.(\d+)([ACGT])>([ACGT])', hgvs_g)
            if m:
                chrom_num = m.group(1)  # e.g., '17'
                pos = int(m.group(2))
                ref = m.group(3)
                alt = m.group(4)

        # Fetch conservation if available
        conservation_score = None
        if self.conservation_fetcher and chrom_num and pos:
            chrom_label = f"chr{chrom_num}"
            conservation_score = self.conservation_fetcher.get_conservation_score(chrom_label, pos)
            if conservation_score is not None:
                print(f"✅ Real conservation score for {chrom_label}:{pos}: {conservation_score:.4f}")

        # Prefer frequency provided by input (e.g., ClinVar AF) if present
        gnomad_freq = 0.0
        possible_freq_keys = [
            'clinvar_af', 'ClinVar_AF', 'AF', 'AlleleFrequency', 'allele_frequency',
            'gnomad_freq', 'gnomAD frequency', 'gnomAD_frequency', 'global_af'
        ]
        for key in possible_freq_keys:
            if key in row and row.get(key) not in (None, ''):
                try:
                    val = str(row.get(key)).strip()
                    val = val.replace('%', '')
                    freq_val = float(val)
                    if freq_val > 1.0:
                        freq_val = freq_val / 100.0  # normalize percent to fraction
                    gnomad_freq = freq_val
                    print(f"🧾 Using provided frequency {gnomad_freq:.6f} from column '{key}'")
                    break
                except Exception:
                    pass

        # If not provided, early frequency lookup using coordinates (prefer coordinates we already have)
        # 🚨 DISABLED: API calls are broken and slow - TSV already has frequencies!
        # if gnomad_freq == 0.0 and chrom_num and pos and ref and alt:
        #     try:
        #         freq_result = self.frequency_fetcher.get_variant_frequency(chrom_num, pos, ref, alt)
        #         gnomad_freq = float(freq_result.get('frequency', 0.0) or 0.0)
        #         print(f"🌍 Pre-analysis frequency for {gene} {variant} at chr{chrom_num}:{pos} {ref}>{alt}: {gnomad_freq:.6f} (source: {freq_result.get('source','?')})")
        #     except Exception as fe:
        #         print(f"⚠️ Frequency lookup failed for {gene} {variant} at chr{chrom_num}:{pos} {ref}>{alt}: {fe}")

        # Derive/collect variant_type from ClinVar-style columns when present
        variant_type = None
        raw_molecular_consequence = ''
        possible_type_keys = [
            'molecular_consequence', 'MolecularConsequence', 'Consequence', 'consequence',
            'MC', 'SequenceOntology', 'sequence_ontology', 'variant_type'
        ]
        for key in possible_type_keys:
            val = row.get(key)
            if val:
                raw_molecular_consequence = str(val)
                v = raw_molecular_consequence.lower()
                # 🚨 CRITICAL: Check for nonsense/frameshift/start_loss FIRST (auto-pathogenic!)
                if 'nonsense' in v or 'stop_gained' in v or 'stop gain' in v:
                    variant_type = 'nonsense'
                    break
                if 'start_lost' in v or 'start loss' in v or 'start_loss' in v or 'initiator_codon' in v:
                    variant_type = 'start_loss'
                    break
                if 'frameshift' in v or 'frame_shift' in v or 'frame shift' in v:
                    variant_type = 'frameshift'
                    break
                if 'stop_lost' in v or 'stop loss' in v or 'stop_loss' in v:
                    variant_type = 'stop_loss'
                    break
                if 'splice' in v:
                    variant_type = 'splice'
                    break
                if 'synonymous' in v or '=)' in v:
                    variant_type = 'synonymous'
                    break
                if 'intron' in v or 'intronic' in v:
                    variant_type = 'intronic'
                    break
                if 'missense' in v:
                    variant_type = 'missense'
                    break

        # 🚨 HEURISTIC: Check variant string itself for nonsense/frameshift/start_loss indicators
        if not variant_type and variant:
            v_lower = str(variant).lower()
            # Start loss: p.Met1?, p.Met1Leu, p.M1?, p.M1L
            if variant.startswith('p.Met1') or variant.startswith('p.M1'):
                variant_type = 'start_loss'
            # Nonsense: p.Gly71*, p.Arg123Ter, p.R123*
            elif '*' in variant or 'ter' in v_lower or 'x' == variant[-1:]:
                variant_type = 'nonsense'
            # Frameshift: p.Arg123fs, p.Arg123Lysfs*10
            elif 'fs' in v_lower:
                variant_type = 'frameshift'
            # Indels: del, ins, dup, delins
            elif any(x in v_lower for x in ['del', 'ins', 'dup']):
                variant_type = 'indel'

        # If we lack a parseable protein HGVS, try converting to protein using the existing converter
        # Gate conversion to coding consequences to avoid noisy network calls for intronic/UTR
        coding_types = {'missense','stop_gained','stop_lost','frameshift','splice','splice_acceptor','splice_donor'}
        should_convert = hgvs_g and (protein_uninformative or (not variant_p)) and (variant_type in coding_types or variant_type is None) and (not os.getenv("CASCADE_DISABLE_CONVERSION"))
        if should_convert:
            try:
                converter = GenomicToProteinConverter()
                prot = converter.convert_genomic_to_protein(hgvs_g, gene_name=gene)
                if prot:
                    variant = prot
                    if not variant_type or variant_type == 'unknown':
                        variant_type = 'missense'
                    print(f"✅ Genomic→protein conversion succeeded for {gene}: {hgvs_g} → {variant}")
                else:
                    # Avoid mislabeling genomic/intractable rows as missense
                    if not variant_type:
                        variant_type = 'unknown'
            except Exception as e:
                print(f"⚠️ Genomic→protein conversion failed for {gene} {hgvs_g}: {e}")
                if not variant_type:
                    variant_type = 'unknown'

        # 🚨 SKIP variants with no protein consequence (intronic, UTR, etc.)
        # These can't be analyzed by cascade and will just waste time
        # Check if variant is still genomic (starts with NC_ or chr) or has no protein info
        is_still_genomic = variant and (variant.startswith('NC_') or variant.startswith('chr'))
        if is_still_genomic or variant_type in {'intronic', 'unknown'}:
            return {
                'skipped': True,
                'skip_reason': 'skipped_no_protein_consequence'
            }

        # Collect ClinVar-relevant passthroughs
        # Support both 'clinical_sig' and 'clinvar_sig' column names
        clinical_sig = row.get('clinical_sig', row.get('clinvar_sig', ''))
        review_status = row.get('review_status', row.get('ReviewStatus', ''))
        rcv = row.get('RCVaccession', row.get('RCV', ''))
        variation_id = row.get('VariationID', row.get('variation_id', ''))

        return {
            'gene': gene,
            'variant': variant,
            'gnomad_freq': gnomad_freq,
            'hgvs': hgvs_g,
            'clinical_sig': clinical_sig,
            'review_status': review_status,
            'molecular_consequence': raw_molecular_consequence,
            'rcv': rcv,
            'variation_id': variation_id,
            'conservation_score': conservation_score,
            'chrom': chrom_num,
            'pos': pos,
            'ref': ref,
            'alt': alt,
            'variant_type': variant_type,
        }

    def get_initial_stats(self):
        return {key: 0 for key in ['total_variants', 'processed', 'skipped_frequency', 'skipped_frameshift', 'skipped_synonymous', 'skipped_intronic', 'skipped_unparseable', 'skipped_no_protein_consequence', 'cascade_triggered', 'dn_only', 'failed', 'success', 'agreements', 'better_data_benign', 'better_data_pathogenic', 'clinical_correlation', 'disagreements', 'unclear']}

    def finalize_processing(self, results, output_path, stats):
        self.write_results_tsv(results, output_path)
        self.frequency_fetcher.save_cache()

    def write_results_tsv(self, results: List[Dict], output_path: str):
        if not results:
            # Write header-only file so downstream consumers don't fail
            # on "file not found" when every row happened to be skipped.
            Path(output_path).write_text("gene\tvariant\tstatus\n", encoding="utf-8")
            print(f"⚠️  No rows scored — wrote header-only file: {output_path}")
            return

        def map_clinvar_bucket(sig: str) -> str:
            if not sig: return 'NA'
            s = sig.strip().replace('_',' ').replace('-',' ').lower()
            # IMPORTANT: Check conflicting FIRST before pathogenic/benign keywords
            if 'conflicting' in s:
                return 'CONFLICTING'
            if 'likely pathogenic' in s:
                return 'LP'
            if 'pathogenic' in s:
                return 'P'
            if 'likely benign' in s:
                return 'LB'
            if 'benign' in s:
                return 'B'
            if 'uncertain' in s or 'vus' in s:
                return 'VUS'
            return 'NA'

        def map_ai_bucket(cls: str) -> str:
            if not cls: return 'NA'
            c = str(cls).upper()
            if c.startswith('VUS-P'): return 'VUS'  # bucket VUS-P with VUS for display bucket
            if c in ('P', 'LP'): return c
            if c in ('B', 'LB'): return c
            if c.startswith('VUS'): return 'VUS'
            return 'NA'

        def agreement_label(clin: str, ai: str) -> str:
            # Align with agreement_analysis.py semantics
            if clin == 'NA' or ai == 'NA':
                return 'NA'

            # CONFLICTING in ClinVar = experts can't agree, so ANY call is reasonable
            if clin == 'CONFLICTING':
                return 'CONFLICTING_OK'

            # VUS <-> VUS agreement
            if clin == 'VUS' and ai == 'VUS':
                return 'AGREE'

            # If ClinVar VUS and AI moves to a side, that's better data
            if clin == 'VUS' and ai in {'B','LB'}:
                return 'BETTER_DATA_to_benign'
            if clin == 'VUS' and ai in {'P','LP'}:
                return 'BETTER_DATA_to_pathogenic'

            # If WE say VUS, we're being appropriately cautious - not disagreeing
            if ai == 'VUS':
                return 'CAUTIOUS'

            # Bucket agreement across sides
            if (clin in {'B','LB'} and ai in {'B','LB'}) or (clin in {'P','LP'} and ai in {'P','LP'}):
                return 'AGREE'

            # Only TRUE disagreement is P/LP <-> B/LB flips
            return 'DISAGREE'

        columns = [
            # Identification
            'gene', 'variant', 'variant_type', 'hgvs', 'gnomad_freq', 'status',
            # ClinVar context
            'clinical_sig', 'review_status', 'molecular_consequence', 'rcv', 'variation_id',
            'clinvar_bucket',
            # Franklin-native context (downstream cumbursum adapter needs these)
            'zygosity', 'confidence',
            # Track 1: RAW MECHANISM (what the molecular analysis actually found)
            'raw_dn', 'raw_lof', 'raw_gof',
            'raw_score', 'raw_classification', 'raw_bucket', 'raw_vs_clinvar',
            # Track 2: ADJUSTED (after plausibility weighting + conservation nudge)
            'adj_dn', 'adj_lof', 'adj_gof',
            'adj_score', 'adj_classification', 'adj_bucket', 'adj_vs_clinvar',
            # The bridge: WHY the tracks differ
            'gene_family',
            'dn_plausibility', 'lof_plausibility', 'gof_plausibility',
            'phylop_score', 'conservation_nudge', 'pre_nudge_classification',
            'atypical_flags',
            'synergy_used',
            # LOF diagnostics
            'base_lof_score', 'stability_impact', 'structural_impact',
            'functional_impact', 'conservation_impact',
            'smart_multiplier', 'domain_multiplier', 'ml_proline_multiplier',
            # Evidence-quality confidence (added 2026-04-27 — see CONFIDENCE_SCORE_DESIGN.md)
            'evidence_confidence', 'evidence_confidence_factors',
            # Explanations
            'explanation', 'review_flags',
        ]

        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for r in results:
                scores = r.get('scores', {}) or {}
                filtered = r.get('plausibility_filtered_scores', {}) or {}

                # Raw mechanism scores (what the analyzers actually found)
                raw_dn = scores.get('DN', 0.0)
                raw_lof = scores.get('LOF', 0.0)
                raw_gof = scores.get('GOF', 0.0)

                # Adjusted scores (after plausibility weighting)
                adj_dn = filtered.get('DN', raw_dn)
                adj_lof = filtered.get('LOF', raw_lof)
                adj_gof = filtered.get('GOF', raw_gof)

                # Per-mechanism plausibility weights from the filter
                plaus_detail = r.get('plausibility_filter_detail', {})
                dn_plaus = plaus_detail.get('DN', {}).get('family_weight', '')
                lof_plaus = plaus_detail.get('LOF', {}).get('family_weight', '')
                gof_plaus = plaus_detail.get('GOF', {}).get('family_weight', '')

                # ClinVar bucketing
                clin_bucket = map_clinvar_bucket(r.get('clinical_sig', ''))

                # Track 1: Raw classification + ClinVar comparison
                raw_bucket = map_ai_bucket(r.get('raw_classification', ''))
                raw_vs_clinvar = agreement_label(clin_bucket, raw_bucket)

                # Track 2: Adjusted classification + ClinVar comparison
                adj_bucket = map_ai_bucket(r.get('final_classification', ''))
                adj_vs_clinvar = agreement_label(clin_bucket, adj_bucket)

                # Atypical mechanism flags (the "future science" signals)
                atypical = r.get('atypical_mechanisms', [])
                atypical_str = '; '.join(
                    f"{a['mechanism']}(weight={a.get('family_weight', '?')},raw={a.get('raw_score', '?'):.2f})"
                    for a in atypical
                ) if atypical else ''

                # LOF diagnostics
                lof_details = r.get('lof_details') or {}

                row = {
                    # Identification
                    'gene': r.get('gene', ''),
                    'variant': r.get('variant', ''),
                    'variant_type': r.get('variant_type', ''),
                    'hgvs': r.get('hgvs', ''),
                    'gnomad_freq': r.get('gnomad_freq', ''),
                    'status': r.get('status', ''),
                    # ClinVar
                    'clinical_sig': r.get('clinical_sig', ''),
                    'review_status': r.get('review_status', ''),
                    'molecular_consequence': r.get('molecular_consequence', ''),
                    'rcv': r.get('rcv', ''),
                    'variation_id': r.get('variation_id', ''),
                    'clinvar_bucket': clin_bucket,
                    # Track 1: Raw
                    'raw_dn': raw_dn,
                    'raw_lof': raw_lof,
                    'raw_gof': raw_gof,
                    'raw_score': r.get('raw_score', ''),
                    'raw_classification': r.get('raw_classification', ''),
                    'raw_bucket': raw_bucket,
                    'raw_vs_clinvar': raw_vs_clinvar,
                    # Track 2: Adjusted
                    'adj_dn': adj_dn,
                    'adj_lof': adj_lof,
                    'adj_gof': adj_gof,
                    'adj_score': r.get('final_score', ''),
                    'adj_classification': r.get('final_classification', ''),
                    'adj_bucket': adj_bucket,
                    'adj_vs_clinvar': adj_vs_clinvar,
                    # Bridge
                    'gene_family': r.get('gene_family', ''),
                    'dn_plausibility': dn_plaus,
                    'lof_plausibility': lof_plaus,
                    'gof_plausibility': gof_plaus,
                    'phylop_score': r.get('phylop_score', ''),
                    'conservation_nudge': r.get('conservation_nudge_applied', ''),
                    'pre_nudge_classification': r.get('pre_nudge_classification', ''),
                    'atypical_flags': atypical_str,
                    'synergy_used': r.get('synergy_used', False),
                    # LOF diagnostics
                    'base_lof_score': lof_details.get('base_lof_score', ''),
                    'stability_impact': lof_details.get('stability_impact', ''),
                    'structural_impact': lof_details.get('structural_impact', ''),
                    'functional_impact': lof_details.get('functional_impact', ''),
                    'conservation_impact': lof_details.get('conservation_impact', ''),
                    'smart_multiplier': lof_details.get('smart_multiplier', ''),
                    'domain_multiplier': lof_details.get('domain_multiplier', ''),
                    'ml_proline_multiplier': lof_details.get('ml_proline_multiplier', ''),
                    # Evidence-quality confidence (computed in cascade_analyzer)
                    'evidence_confidence': r.get('evidence_confidence', ''),
                    'evidence_confidence_factors': r.get('evidence_confidence_factors', ''),
                    # Explanations
                    'explanation': r.get('explanation', ''),
                    'review_flags': r.get('review_flags', ''),
                }
                writer.writerow(row)
        print(f"💾 CASCADE results saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="🌊 CASCADE BATCH PROCESSOR")
    parser.add_argument('--input', required=True, help='Input CSV/TSV file path')
    parser.add_argument('--output', required=True, help='Output TSV file path')
    parser.add_argument('--ignore-quality-floor', action='store_true',
                        help='Keep variants with Quality < 20 (normally dropped as likely sequencing '
                             'noise). Use for research paths that downweight rather than drop Low-conf '
                             'variants (e.g. BSZ v2 option-D weighting). Default off — status quo.')
    args = parser.parse_args()
    output_dir = os.path.dirname(args.output)
    if output_dir: os.makedirs(output_dir, exist_ok=True)
    processor = CascadeBatchProcessor(ignore_quality_floor=args.ignore_quality_floor)
    processor.process_csv(args.input, args.output)

if __name__ == "__main__":
    main()
