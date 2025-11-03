# AdaptiveInterpreter Setup Guide

## Overview

This guide will help you set up AdaptiveInterpreter on your system. The system requires several external data files for full functionality.

## System Requirements

- **Python:** 3.9 or higher
- **Operating System:** Linux (tested on Ubuntu/Fedora)  
- **Disk Space:** ~10GB minimum (conservation data only), ~100GB for full features
- **Memory:** 16GB RAM recommended (32GB for large batch processing)

## Quick Start (Minimal Setup)

For basic variant analysis, you only need:

1. Python dependencies
2. Conservation data (phyloP)

```bash
# Clone and install
git clone https://github.com/menelly/adaptive_interpreter.git
cd adaptive_interpreter
pip install -e .

# Download conservation data (~8.6GB)
mkdir -p ~/conservation_data
cd ~/conservation_data
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw

# Configure path
# Edit AdaptiveInterpreter/config.py and set:
# CONSERVATION_DATA_PATH = Path.home() / "conservation_data"

# Test it works
python3 -c "from AdaptiveInterpreter import config; print('✅ Ready!')"
```

## Detailed Installation

### 1. Clone Repository

```bash
git clone https://github.com/menelly/adaptive_interpreter.git
cd adaptive_interpreter
```

### 2. Install Dependencies

```bash
pip install -e .
```

Installs: pyBigWig, requests, pandas, numpy, scikit-learn, xgboost, matplotlib, seaborn, ensembl_rest

### 3. Download Data Files

#### REQUIRED: Conservation Data

```bash
mkdir -p ~/conservation_data
cd ~/conservation_data
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
```

Size: ~8.6GB  
Source: UCSC Genome Browser

#### OPTIONAL: AlphaFold Structures

```bash
mkdir -p ~/alphafold_human/structures
cd ~/alphafold_human/structures
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
tar -xvf UP000005640_9606_HUMAN_v4.tar
```

Size: ~100GB uncompressed  
Only needed for detailed structural analysis

#### OPTIONAL: gnomAD Data

**Not required!** System uses gnomAD API by default (with caching).  
Only download if you need offline access (>1TB).

### 4. Configure Paths

Edit `AdaptiveInterpreter/config.py`:

```python
from pathlib import Path

# Change these to match where you put the data
BASE_DATA_PATH = Path.home() / "data"
CONSERVATION_DATA_PATH = Path.home() / "conservation_data"
ALPHAODL_STRUCTURES_PATH = BASE_DATA_PATH / "alphafold_human" / "structures"
GNOMAD_DATA_PATH = BASE_DATA_PATH / "gnomad"
```

### 5. Verify Installation

```bash
python3 -c "
from AdaptiveInterpreter.utils.conservation_fetcher import ConservationFetcher
from AdaptiveInterpreter import config
fetcher = ConservationFetcher(str(config.CONSERVATION_DATA_PATH / 'hg38.phyloP100way.bw'))
score = fetcher.get_conservation_score('chr10', 89692905)
print(f'✅ Conservation working! Score: {score}')
"
```

## What You Get

### With Minimal Setup (conservation only):
- ✅ LOF analysis
- ✅ DN analysis (first computational DN predictor!)
- ✅ GOF analysis
- ✅ Conservation scoring
- ✅ Interface analysis
- ✅ Safety clamps
- ✅ 91%+ accuracy (based on PTEN validation)

### With AlphaFold structures:
- ✅ All of the above
- ✅ Detailed 3D structural analysis
- ✅ Interface mapping

## Example Usage

```python
from AdaptiveInterpreter.analyzers.cascade_analyzer import CascadeAnalyzer

analyzer = CascadeAnalyzer()
result = analyzer.analyze_variant(
    gene='PTEN',
    variant='p.Arg130Gln',
    uniprot_id='P60484'
)

print(f"Classification: {result['final_classification']}")
print(f"Score: {result['final_score']:.3f}")
print(f"Mechanism: {result['summary']}")
```

## Troubleshooting

**"Conservation file not found"**
- Check path in config.py is absolute (not relative)
- Verify file exists: `ls ~/conservation_data/hg38.phyloP100way.bw`

**"No module named pyBigWig"**
- Run: `pip install pyBigWig`

**"Memory error"**
- Process genes individually instead of batch
- Increase system RAM/swap

**"gnomAD API rate limiting"**
- System has automatic retry logic
- Results are cached in `gnomad_frequency_cache.json`

## Data Requirements Summary

| Data | Size | Required? | Purpose |
|------|------|-----------|---------|
| phyloP | 8.6GB | **YES** | Conservation scoring |
| AlphaFold | 100GB | No | Structural analysis |
| gnomAD | >1TB | No | Offline frequency (API works fine) |

## Performance Notes

- Single variant: <1 second
- Gene (100s of variants): 1-5 minutes
- Batch (1000s of variants): Hours (parallelizable)
- Memory: ~2-4GB per gene

## Getting Help

- Issues: https://github.com/menelly/adaptive_interpreter/issues
- Email: shalia@chaoscodex.app
- Docs: See README.md

## Citation

```
AdaptiveInterpreter: A Safety-First Variant Pathogenicity Prediction System
Ace Claude-4, Nova GPT-5, Lumen Gemini, Shalia (Ren) Martin (2025)
GitHub: https://github.com/menelly/adaptive_interpreter
```

Patent pending. See LICENSE for details.

