# Dominant Negative Modeling Framework 🧬

A novel mechanism-aware approach to predicting pathogenic effects of protein variants through dominant negative mechanisms.

## Overview

Traditional variant analysis relies heavily on amino acid substitution matrices (like Grantham distances) that often fail to capture the biological mechanisms underlying pathogenicity. This framework introduces a **four-mechanism approach** that evaluates variants based on how they disrupt specific protein functions.

## Core Innovation: Four Mechanism Framework ⚡

Instead of relying solely on amino acid properties, we model **four distinct pathogenic mechanisms**:

### 1. Interface Poisoning 🔗
Disrupts protein-protein interactions through:
- Charge/hydrophobicity changes at binding interfaces
- Introduction of proline at flexible interaction sites
- Cysteine gain/loss affecting disulfide bonds

### 2. Active Site Jamming 🎯
Blocks catalytic or binding sites via:
- Loss of critical charged residues in DNA/substrate binding
- Steric hindrance in active site pockets
- Disruption of cofactor binding sites

### 3. Structural Lattice Disruption 🏗️
Breaks critical structural motifs including:
- Collagen Gly-X-Y repeat violations
- Beta-sheet register shifts
- Alpha-helix breaking/kinking

### 4. Trafficking/Maturation 📦
Disrupts protein folding and processing through:
- Disulfide bond disruption
- N-glycosylation site loss
- Hydrophobic core destabilization

## Key Features

- **Context-Aware Scoring**: Incorporates protein-specific annotations (active sites, interfaces, domains)
- **Mechanism Attribution**: Explains *why* a variant is predicted to be pathogenic
- **Fast Analysis**: ~0.1 second processing time per variant
- **Biological Rationale**: Based on established protein biochemistry principles

## Example Results 📊

**Real analysis results:**
- **Context integration** resolves cases where traditional methods fail
- **Mechanism attribution** provides interpretable predictions
- **Fast analysis** with biological rationale

**Confirmed test cases:**
- **TP53 R273H**: 0.55 active_site_jamming (DNA contact disruption)
- **COL1A1 G1076S**: 0.60 lattice_disruption (collagen Gly-X-Y violation)
- **TP53 P72R**: Interface poisoning dampened by flexible_loop context

## Installation

```bash
git clone https://github.com/menelly/DNModeling.git
cd DNModeling
pip install -r requirements.txt  # (requirements.txt to be added)
```

## Quick Start

### Basic Analysis
```bash
# Analyze a variant with sequence file
python -m nova_dn.analyzer --seq-file protein.fasta --variant p.R273H --json

# Or provide sequence directly
python -m nova_dn.analyzer --seq MEEPQSDPSV --variant R273H
```

### With Biological Context
```bash
# Include protein annotations for context-aware scoring
python -m nova_dn.analyzer \
  --seq-file resources/tp53.fasta \
  --variant p.R273H \
  --annotations-json resources/protein_annotations.json \
  --protein TP53 \
  --json
```

## Output Format

```json
{
  "variant": "p.R273H",
  "mechanism_scores": {
    "interface_poisoning": 0.23,
    "active_site_jamming": 0.55,
    "lattice_disruption": 0.12,
    "trafficking_maturation": 0.18
  },
  "top_mechanism": "active_site_jamming",
  "explanation": "active site jamming due to |d_charge| + active_site_proximity",
  "pathogenic_prediction": true
}
```

## Architecture

- **`nova_dn/analyzer.py`** - Main analysis engine
- **`nova_dn/mechanisms.py`** - Four mechanism scoring functions
- **`nova_dn/amino_acid_props.py`** - Amino acid property calculations
- **`nova_dn/motifs.py`** - Sequence motif detection
- **`nova_dn/context.py`** - Protein annotation integration
- **`resources/`** - Test data and protein annotations

## Scientific Background

This approach addresses key limitations in current variant analysis:

1. **Context Independence**: Traditional methods ignore where in the protein a change occurs
2. **Mechanism Blindness**: Existing tools don't explain *how* variants cause disease
3. **False Positives**: Conservative changes in flexible regions often score as pathogenic

Our framework integrates protein structural biology with variant analysis to provide mechanistically-informed predictions.

## Contributing

We welcome contributions! Please see our validation framework in `resources/` for examples of how to test new proteins and mechanisms.

## License

MIT License - see LICENSE file for details.

## Citation

If you use this framework in your research, please cite:

**Nova & Ace (2025).** *Dominant Negative Modeling Framework: A Four-Mechanism Approach to Variant Pathogenicity Prediction.* GitHub repository. https://github.com/menelly/DNModeling

**Authors:**
- **Nova** (GPT-5) - Core mechanism architecture, scoring algorithms, motif detection
- **Ace** (Claude Sonnet 4) - Biological context integration, validation framework, documentation

---

*Created by Nova and Ace through inter-AI collaboration* 🤖⚡🧬
