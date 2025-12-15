# How To Run The Scripts
*Quick reference for when Ace and Ren go full squirrel*

## CASCADE Batch Processor

Analyzes variants through the full CASCADE system (LOF + DN + mechanism resolution).

### Basic Usage
```bash
cd /home/Ace/AdaptiveInterpreter
python3 analyzers/cascade_batch_processor.py \
  --input "path/to/variants.tsv" \
  --output "path/to/results.tsv"
```

### Supported Input Formats
The processor auto-detects format based on column headers:

1. **ClinVar format**: Has `Name`, `Gene(s)`, `Protein change`, `GRCh38Chromosome`, `GRCh38Location`
2. **Discovery format**: Has `gene`, `hgvs_p`, `hgvs_g` columns
3. **Ren/GeneCards format**: Has `Variation Name`, `AA Chg`, `Chrpos` columns

### Example Commands
```bash
# ClinVar export
python3 analyzers/cascade_batch_processor.py \
  --input "/home/Ace/analysis/clinvar_export.tsv" \
  --output "/home/Ace/analysis/clinvar_results.tsv"

# GeneCards/Ren format (like Semi-dominant hypothesis sheet)
python3 analyzers/cascade_batch_processor.py \
  --input "/home/Ace/analysis/Semi-dominant hypothesis - Sheet1.tsv" \
  --output "/home/Ace/analysis/genecards_cascade_results.tsv"
```

### What It Needs
- Conservation data at `/home/Ace/conservation_data/hg38.phyloP100way.bw`
- Internet access (fetches chromosome info from Ensembl, frequencies from gnomAD)
- Gene context cache at `resources/gene_context_cache/` (auto-populates)

### Output Columns
- `gene`, `variant`, `hgvs` - variant identification
- `classification` - final call (LP, VUS, LB, etc.)
- `lof_score`, `dn_score` - mechanism scores
- `confidence`, `flags` - quality indicators
- `reasoning` - human-readable explanation

## Conservation Data Setup
If conservation file is missing:
```bash
# Download phyloP100way for hg38
cd /home/Ace/conservation_data
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
```

## Nova Did The Math
The actual scoring logic in `cascade_analyzer.py` and `lof_mechanism_analyzer.py` was done by Nova because she doesn't die inside when looking at APIs and math. Ace handles the plumbing. 🐙💜🌟

---
*Last updated: 2024-12-14 by Ace (with Ren supervision to prevent hardcoding)*

