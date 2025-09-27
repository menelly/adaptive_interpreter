# HOW TO USE THE DN ANALYZER 🧬⚡
## *"Stop Panicking and Read This First!" - Ren & Ace*

### 🚨 CRITICAL: The Magic Incantation That Actually Works

**The CORRECT way to get results that match REAL_TEST_RESULTS.md:**

```bash
python3 -m nova_dn.analyzer \
  --seq-file /tmp/tp53_alphafold.fasta \
  --variant p.R248Q \
  --annotations-json resources/protein_annotations.json \
  --protein TP53 \
  --json
```

**WITHOUT the annotations, you get WRONG results!** 🚨

---

## 🔬 Step-by-Step Process

### Step 1: Get the AlphaFold Sequence
```python
from nova_dn.alphafold_sequence import AlphaFoldSequenceExtractor

extractor = AlphaFoldSequenceExtractor()
# Save TP53 sequence to file
fasta_path = '/tmp/tp53_alphafold.fasta'
extractor.save_fasta('P04637', fasta_path, 'TP53')  # P04637 = TP53 UniProt ID
```

### Step 2: Run the Analyzer with ALL Required Flags
```bash
python3 -m nova_dn.analyzer \
  --seq-file /tmp/tp53_alphafold.fasta \
  --variant p.R248Q \
  --annotations-json resources/protein_annotations.json \
  --protein TP53 \
  --json
```

---

## 🎯 Flag Explanations

| Flag | Purpose | Example | REQUIRED? |
|------|---------|---------|-----------|
| `--seq-file` | Path to FASTA file with protein sequence | `/tmp/tp53_alphafold.fasta` | ✅ YES |
| `--variant` | Variant in p.RefPosAlt format | `p.R248Q` | ✅ YES |
| `--annotations-json` | Protein annotations for context | `resources/protein_annotations.json` | ✅ **CRITICAL!** |
| `--protein` | Protein name for annotation lookup | `TP53` | ✅ **CRITICAL!** |
| `--json` | Output in JSON format | (no value) | ✅ YES |

### 🚨 Why Annotations Are CRITICAL

**Without annotations:**
- R248Q gets active_site_jamming = 0.0 ❌
- P72R gets interface_poisoning = 0.414 ❌

**With annotations:**
- R248Q gets active_site_jamming = 0.4 ✅ (from `active_site_proximity` context)
- P72R gets interface_poisoning = 0.314 ✅ (from `flexible_loop` dampening)

---

## 📊 Expected Results (Validation)

### TP53 R248Q (Pathogenic)
```json
{
  "mechanism_scores": {
    "interface_poisoning": 0.372,
    "active_site_jamming": 0.4,
    "lattice_disruption": 0.0,
    "trafficking_maturation": 0.0
  },
  "top_mechanism": "active_site_jamming"
}
```

### TP53 P72R (Benign)
```json
{
  "mechanism_scores": {
    "interface_poisoning": 0.314,
    "active_site_jamming": 0.0,
    "lattice_disruption": 0.0,
    "trafficking_maturation": 0.0
  },
  "top_mechanism": "interface_poisoning"
}
```

---

## 🧬 Common UniProt IDs

| Gene | UniProt ID | Notes |
|------|------------|-------|
| TP53 | P04637 | Tumor suppressor, tetramer |
| COL1A1 | P02452 | Collagen, trimer |
| FGFR3 | P22607 | Growth factor receptor, dimer |
| VWF | P04275 | von Willebrand factor, multimer |
| RYR1 | P21817 | Calcium channel, massive protein |

---

## 🔥 Quick Test Commands

### Test TP53 Variants
```bash
# Get sequence
python3 -c "from nova_dn.alphafold_sequence import AlphaFoldSequenceExtractor; extractor = AlphaFoldSequenceExtractor(); extractor.save_fasta('P04637', '/tmp/tp53.fasta', 'TP53')"

# Test pathogenic
python3 -m nova_dn.analyzer --seq-file /tmp/tp53.fasta --variant p.R248Q --annotations-json resources/protein_annotations.json --protein TP53 --json

# Test benign
python3 -m nova_dn.analyzer --seq-file /tmp/tp53.fasta --variant p.P72R --annotations-json resources/protein_annotations.json --protein TP53 --json
```

---

## 🚨 TROUBLESHOOTING

### "Why am I getting wrong scores?"
- ✅ Check: Are you using `--annotations-json` and `--protein` flags?
- ✅ Check: Are you using AlphaFold sequences, not local FASTA files?
- ✅ Check: Is your variant format `p.R248Q` not `R248Q`?

### "Why is my pathogenic variant scoring low?"
- Missing `active_site_proximity` context from annotations
- Check that the protein name matches what's in protein_annotations.json

### "Why is my benign variant scoring high?"
- Missing `flexible_loop` or other dampening context from annotations
- The variant might actually be in a functionally important region

---

## 💜 Remember This File Exists!

**File name: `HOW_TO_USE_THE_ANALYZER.md`**

Next time you're panicking about wrong scores, READ THIS FIRST! 🧬⚡💜
