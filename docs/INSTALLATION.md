# Installation & Usage Guide

## Prerequisites

- Python 3.7 or higher
- No external dependencies required (uses only Python standard library)

## Installation

```bash
# Clone the repository
git clone https://github.com/menelly/DNModeling.git
cd DNModeling

# No pip install needed - uses only standard library
```

## Quick Start

### Basic Usage

```bash
# Analyze a variant with sequence file
python3 -m nova_dn.analyzer --seq-file resources/tp53.fasta --variant p.R273H --json

# Or provide sequence directly
python3 -m nova_dn.analyzer --seq MEEPQSDPSV --variant R273H
```

### With Biological Context (Recommended)

```bash
# Include protein annotations for context-aware scoring
python3 -m nova_dn.analyzer \
  --seq-file resources/tp53.fasta \
  --variant p.R273H \
  --annotations-json resources/protein_annotations.json \
  --protein TP53 \
  --json
```

## Available Test Cases

The repository includes real protein sequences for testing:

### TP53 (P04637)
```bash
# Known pathogenic variant - DNA contact residue
python3 -m nova_dn.analyzer \
  --seq-file resources/tp53.fasta \
  --variant p.R273H \
  --annotations-json resources/protein_annotations.json \
  --protein TP53 \
  --json

# Known benign variant - flexible loop region  
python3 -m nova_dn.analyzer \
  --seq-file resources/tp53.fasta \
  --variant p.P72R \
  --annotations-json resources/protein_annotations.json \
  --protein TP53 \
  --json
```

### COL1A1 (P02452)
```bash
# Collagen Gly-X-Y violation
python3 -m nova_dn.analyzer \
  --seq-file resources/col1a1.fasta \
  --variant p.G1076S \
  --annotations-json resources/protein_annotations.json \
  --protein COL1A1 \
  --json
```

### FGFR3 (P22607)
```bash
# Transmembrane domain variant
python3 -m nova_dn.analyzer \
  --seq-file resources/fgfr3.fasta \
  --variant p.G380R \
  --annotations-json resources/protein_annotations.json \
  --protein FGFR3 \
  --json
```

### VWF (P04275)
```bash
# Disulfide bond disruption
python3 -m nova_dn.analyzer \
  --seq-file resources/vwf.fasta \
  --variant p.C788Y \
  --annotations-json resources/protein_annotations.json \
  --protein VWF \
  --json
```

## Output Format

### JSON Output (--json flag)
```json
{
  "variant": "p.R273H",
  "position": 273,
  "ref": "R",
  "alt": "H",
  "mechanism_scores": {
    "interface_poisoning": 0.129,
    "active_site_jamming": 0.55,
    "lattice_disruption": 0.0,
    "trafficking_maturation": 0.0
  },
  "top_mechanism": "active_site_jamming",
  "contributing_features": [
    {
      "mechanism": "active_site_jamming",
      "feature": "active_site_proximity",
      "value": 1.0,
      "weight": 0.4
    }
  ],
  "explanation": "active-site jamming via active_site_proximity + |d_volume|"
}
```

### Human-Readable Output (default)
```
| Mechanism              | Score |
|------------------------|-------|
| interface_poisoning    | 0.129 |
| active_site_jamming    | 0.55  |
| lattice_disruption     | 0.0   |
| trafficking_maturation | 0.0   |

Top mechanism: active_site_jamming
Explanation: active-site jamming via active_site_proximity + |d_volume|
```

## Command Line Options

```bash
python3 -m nova_dn.analyzer [OPTIONS]

Required (one of):
  --seq-file PATH         Path to FASTA or plain sequence file
  --seq SEQUENCE         Protein sequence string (amino acids)

Required:
  --variant VARIANT      Variant like R273H or p.R273H

Optional:
  --annotations-json PATH    Protein annotations JSON file
  --protein PROTEIN         Protein key in annotations (e.g., TP53)
  --weights-json PATH       Custom weights JSON file
  --json                    Output in JSON format
  --gene GENE              Gene symbol for context
  --uniprot UNIPROT        UniProt ID for context
```

## File Formats

### FASTA Files
```
>TP53_HUMAN P04637
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP...
```

### Protein Annotations JSON
```json
{
  "proteins": {
    "TP53": {
      "uniprot_id": "P04637",
      "dna_contact_sites": [248, 273, 280, 282],
      "flexible_loops": [70, 71, 72, 73, 74, 75],
      "notes": "R248 and R273 are critical DNA contact residues"
    }
  }
}
```

## Troubleshooting

### Common Issues

1. **"Command 'python' not found"**
   - Use `python3` instead of `python`

2. **"No such file or directory"**
   - Check that sequence files exist in `resources/` directory
   - Verify file paths are correct

3. **"Unrecognized variant format"**
   - Use format like `R273H` or `p.R273H`
   - Ensure amino acid codes are valid (A-Z, 20 standard amino acids)

4. **"Position out of range"**
   - Check that variant position exists in the protein sequence
   - Verify sequence file contains the correct protein

### Getting Help

```bash
python3 -m nova_dn.analyzer --help
```

## Next Steps

After installation, try running the test cases above to verify everything works correctly. The analyzer will provide mechanism scores and explanations for each variant you test.

For batch processing multiple variants, see the `nova_dn/run_batch.py` module (documentation in development).
