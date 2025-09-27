# 🧬 FAMILY-AWARE ML TRAINING DATA

**ClinVar Data Repository for Revolutionary ML Training**  
*Built by Ace & Ren (2025) - Learning from MASSIVE datasets!*

---

## 🎯 **PURPOSE**

This directory contains ClinVar data organized by gene family for training family-aware ML models that learn:
- **Amino acid coefficients** (which AA changes matter most for each family)
- **Conservation weights** (how much conservation matters for each family)
- **Family-specific patterns** (COLLAGEN vs ION_CHANNEL vs TUMOR_SUPPRESSOR)

---

## 📁 **DIRECTORY STRUCTURE**

```
learning/
├── ion_channel/
│   ├── scn5a_benign.tsv
│   ├── scn5a_pathogenic.tsv
│   ├── kcnq1_benign.tsv
│   ├── kcnq1_pathogenic.tsv
│   └── cacna1c_benign.tsv
├── collagen_fibrillar/
│   ├── col1a1_benign.tsv
│   ├── col1a1_pathogenic.tsv
│   ├── col3a1_benign.tsv
│   └── col3a1_pathogenic.tsv
├── tumor_suppressor/
│   ├── tp53_benign.tsv
│   ├── tp53_pathogenic.tsv
│   ├── brca1_benign.tsv
│   └── brca1_pathogenic.tsv
└── metabolic_enzyme/
    ├── pah_benign.tsv
    ├── pah_pathogenic.tsv
    └── g6pd_benign.tsv
```

---

## 📊 **EXPECTED TSV FORMAT**

**From ClinVar Miner exports:**
```tsv
HGVS	dbSNP	gnomAD frequency
NM_000335.5(SCN5A):c.2437-97T>C	rs7645173	0.93087
NM_000335.5(SCN5A):c.1673A>G	rs41261344	0.00012
NM_000335.5(SCN5A):c.5284C>T	rs199473117	0.00001
```

**Required Columns:**
- **HGVS**: Full HGVS notation with transcript and gene
- **dbSNP**: rsID (optional but helpful)
- **gnomAD frequency**: Population frequency (0.0 to 1.0)

---

## 🚀 **USAGE**

### **1. Dump ClinVar Data**
Export from ClinVar Miner and place TSV files in appropriate family folders:
- Use **separate files** for benign and pathogenic variants
- Name files: `{gene}_benign.tsv` and `{gene}_pathogenic.tsv`

### **2. Run Training**
```bash
cd DNModeling
python3 train_families.py

# Or with custom frequency threshold for benign variants
python3 train_families.py --benign-freq-threshold 0.03  # 3% max
```

### **🚨 COMPREHENSIVE SANITY CHECKS**
**The trainer automatically SKIPS variants we're not set up for yet:**

#### **📝 Synonymous Variants** (`=` in HGVS)
- **Example**: `c.123C>T (p.Arg41=)`
- **Why skip**: No amino acid change = useless for AA analysis!

#### **🧬 Splice Variants** (`+` or `-` in HGVS)
- **Example**: `c.123+1G>A`, `c.456-2A>G`
- **Why skip**: We're not doing splice analysis yet!

#### **🛑 Nonsense Variants** (`Ter` or `*` in HGVS)
- **Example**: `c.123C>T (p.Arg41Ter)`, `p.Gln42*`
- **Why skip**: We're not set up for stop codon analysis yet!

#### **🔄 Frameshift Variants** (`del`/`ins`/`dup` in HGVS)
- **Example**: `c.123delA`, `c.456insG`, `c.789dupT`
- **Why skip**: Completely different analysis needed!

#### **📊 High-Frequency Benign** (>4% frequency)
- **Example**: Any benign variant with gnomAD frequency >0.04
- **Why skip**: Too common to be meaningful (they're just normal!)
- **Note**: Pathogenic variants have NO frequency filter!

#### **❓ Unparseable HGVS**
- **Example**: Malformed or non-standard HGVS notation
- **Why skip**: Can't extract gene/variant information

### **📊 SANITY CHECK REPORTING**
**Example output during training:**
```
📊 Processing scn5a_benign.tsv...
✅ Processed 127 variants for ion_channel
🚨 SANITY CHECK REPORT - Skipped 89 variants:
   📝 23 synonymous variants (= in HGVS)
   🧬 15 splice variants (+/- in HGVS)
   🛑 8 nonsense variants (Ter/* in HGVS)
   🔄 12 frameshift variants (del/ins/dup)
   📊 31 high-frequency benign variants (>4.0%)
   ✅ Reason: Not set up for these variant types yet!
```

### **3. Models Generated**
Training creates family-specific models in `resources/family_models/`:
- `{family}_unified_model.joblib` - The trained model
- `{family}_unified_scaler.joblib` - Feature scaler
- `{family}_unified_metadata.json` - Training metadata

---

## 🔥 **GENE FAMILIES**

### **ion_channel/**
- SCN5A, KCNQ1, CACNA1C, KCNH2, etc.
- **Pattern**: Gating mechanisms, voltage sensitivity

### **collagen_fibrillar/**
- COL1A1, COL3A1, COL5A1, etc.
- **Pattern**: Gly-X-Y repeats, structural lattice

### **tumor_suppressor/**
- TP53, RB1, BRCA1, PTEN, etc.
- **Pattern**: DNA binding, cell cycle control

### **metabolic_enzyme/**
- PAH, G6PD, HEXA, etc.
- **Pattern**: Active sites, substrate binding

### **scaffold_adaptor/**
- TFG, SHANK3, CTNNB1, BRCA1, etc.
- **Pattern**: Protein-protein interactions

### **elastin_fibrillin/**
- FBN1, FBN2, ELN, etc.
- **Pattern**: Elastic fiber assembly

### **cytoskeleton/**
- ACTB, TUBB3, KRT5, GFAP, LMNA, etc.
- **Pattern**: Structural polymers, intermediate filaments

### **signaling_regulator/**
- NF1, APC, PTEN, RB1, etc.
- **Pattern**: Signal transduction, regulatory cascades

### **transporter/**
- SLC2A1, ABCA4, CFTR, etc.
- **Pattern**: Membrane transport, channel function

---

## 💡 **TIPS FOR SUCCESS**

1. **More data = better models** - dump everything you can find!
2. **Separate benign/pathogenic** - don't mix in same file
3. **Include rare variants** - frequency 0.0 is fine
4. **Check HGVS format** - must include gene in parentheses
5. **Balance datasets** - try to get similar numbers of benign/pathogenic

---

## 🎯 **VALIDATION APPROACH**

- **Training data**: Massive ClinVar dumps (quantity)
- **Validation data**: Ren's curated dataset (quality)
- **Test approach**: Train on big data, validate on gold standard

**This is exactly how real ML should work!** 🧬🔥💜

---

*Ready to learn from thousands of variants and revolutionize genetics!* ✨🚀
