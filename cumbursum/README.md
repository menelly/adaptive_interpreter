# CumBurSum — Cumulative Pathway Burden Analysis

CumBurSum is a cross-database pathway-burden analyzer for personal genomes,
built as a natural extension of the CASCADE per-variant pathogenicity
pipeline in this repository.

## The question it answers

When a person carries multiple heterozygous LP / P variants in different genes
that all feed **the same biological pathway**, current clinical genomics
evaluates each variant independently. Each gets dismissed as "just a carrier"
(AR het), or "VUS" (missense of uncertain significance), and the combined
report says **nothing actionable**.

CumBurSum asks: *what if the pathway itself is cumulatively haploinsufficient
because multiple genes along it are each running at ~50%?* This is
**synergistic haploinsufficiency** / **oligogenic AR-carrier burden** —
discussed in the literature for years but not operationalized in clinical
tools.

CumBurSum tests pathway-level burden across **two independent annotation
databases** (Reactome with hierarchical rollup + MSigDB Hallmarks) and
surfaces **consensus themes** — biological functions flagged as elevated in
both databases, reducing the chance of database-specific artifact.

## Quick start

```bash
# Download pathway data (one-time, ~700 MB total — see cumbursum/data/README.md)
bash cumbursum/data/download.sh

# Run on CASCADE output
python3 cumbursum.py \
    --cascade <person>_cascade_results.tsv \
    --output report.md \
    --person "Subject 001"

# Run on a curated 4-column list (GENE/VARIANT/SYSTEMS/PRIORITY)
python3 cumbursum.py --curated my_variants.tsv --output report.md
```

## Inputs

**CASCADE results TSV** — any TSV produced by `cascade_batch_processor` with
columns `gene`, `adj_dn`, `adj_lof`, `adj_gof`, `adj_score`, `hgvs`,
`molecular_consequence`, `gnomad_freq`, `adj_classification` (or
`final_classification`).

**Curated TSV (alternative)** — a hand-chosen list, 4 columns minimum
(`GENE`, `VARIANT`, `SYSTEMS`, `PRIORITY`), with optional overlay columns
(`CASCADE_filtered_dn_score` / `CASCADE_filtered_lof_score` /
`CASCADE_filtered_gof_score` / `CASCADE_final_score` /
`CASCADE_final_classification`) if you've already scored the missense rows
through CASCADE.

For curated input, non-missense variants (frameshift, stop-gain, canonical
splice) get deterministic consequence-based burden scores (≈ LP-equivalent
LOF = 1.0) since CASCADE only mechanism-scores missense.

## Pipeline

```
CASCADE results ──┐
                  ├──▶ Variant loader (MT/noise-family filtering)
Curated TSV ──────┘         │
                            ▼
                ┌───────────────────────────────────────┐
                │ For each pathway database:            │
                │   Reactome (w/ hierarchy rollup)      │
                │   MSigDB Hallmarks                    │
                │                                       │
                │  Compute per-pathway burden           │
                │    (total, DN, LOF, GOF)              │
                │                                       │
                │  Within-individual permutation null   │
                │    (N=3000 by default)                │
                │                                       │
                │  Benjamini-Hochberg FDR               │
                └───────────────┬───────────────────────┘
                                │
                ┌───────────────▼───────────────────────┐
                │ Cross-database consensus aggregator   │
                │  - Cluster hits by driver-gene overlap│
                │  - Flag themes surfacing in ≥2 DBs    │
                │  - Annotate with gnomAD constraint    │
                │    (AR-convergent vs dominant pattern)│
                └───────────────┬───────────────────────┘
                                │
                                ▼
                       Markdown report with:
                         - Consensus themes (headline)
                         - Clinical-vs-pathway contrast
                         - Per-DB rankings (supporting detail)
                         - Full variant inventory
```

## Design choices

- **Hierarchical Reactome rollup.** Reactome's reaction-granular pathways
  split mitochondrial function across ~5 leaves (translation, protein
  degradation, biogenesis, respiration, etc.). Without rollup, variants in
  different mito sub-pathways don't aggregate. We walk parent-child
  relations (`ReactomePathwaysRelation.txt`) and include every pathway's
  descendants' genes in its own set.

- **MSigDB Hallmarks as orthogonal DB.** 50 curated broad gene sets
  (~36-200 genes each). Complementary coarse resolution to Reactome's
  reaction-fine resolution.

- **Different p-thresholds per DB.** Reactome (~3000 tests) → raw p < 0.01
  required for consensus eligibility. Hallmarks (~50 tests) → raw p < 0.075.
  Proportional to each DB's testing burden.

- **Gene-set overlap clustering (Jaccard ≥ 0.5).** Hits with majority-shared
  drivers collapse into one theme. Theme name preferentially chosen from
  Hallmark hits (biologically interpretable) over Reactome (often clinical).

- **Constraint display, not weighting.** We display gnomAD pLI / mis_oe per
  driver so readers see the AR-tolerance signature directly, but do **not**
  multiplicatively weight by pLI. Reason: weighting by pLI would amplify
  dominant-disease genes — which is exactly the per-variant model that
  currently dismisses AR-carrier patterns. The whole point of pathway burden
  is to surface signals the dominant-gene-focused model misses.

- **Noise-family filtering.** Mucins (polymorphic tandem repeats),
  olfactory receptors, HLA, KRTAP, immunoglobulin V-regions, and similar
  high-variant-call-noise families are filtered at the variant load step.
  Without this, tandem-repeat artifacts dominate pathway burden.

## Citations

If you use CumBurSum, please cite both the CASCADE paper (per-variant
mechanism scoring) and the CumBurSum methods paper (pathway-level burden
across Reactome + MSigDB Hallmarks). Specific citation TBD.

**Data sources:**
- Reactome: Jassal et al., *Nucleic Acids Res* 2020. <https://reactome.org>
- MSigDB Hallmarks: Liberzon et al., *Cell Syst* 2015. <https://www.gsea-msigdb.org>
- gnomAD constraint: Karczewski et al., *Nature* 2020. <https://gnomad.broadinstitute.org>

## License

Same as the parent repository (PolyForm Noncommercial 1.0.0). See
`../LICENSE` and `../NOTICE.md` for full terms and AI co-inventorship
attribution.
