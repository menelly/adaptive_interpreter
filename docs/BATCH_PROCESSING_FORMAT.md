# 📊 BATCH PROCESSING FILE FORMATS

**CSV/TSV Input Format Guide for AdaptiveInterpreter Batch Processors**  
*Built by Ace & Nova (2025) - Making batch processing crystal clear!*

---

## 🎯 **OVERVIEW**

The AdaptiveInterpreter system supports multiple batch processors with different CSV/TSV input formats. This guide covers all supported formats and their requirements.

---

## 🌊 **CASCADE BATCH PROCESSOR** (Recommended)

**File:** `cascade_batch_processor.py`  
**Purpose:** Full DN → LOF → GOF cascade analysis with biological routing

### **Standard Format (Simple)**
```csv
gene,variant,hgvs,gnomad_freq
COL1A1,p.G1340S,NM_000088.4(COL1A1):c.4018G>A,0.00064
FBN1,p.R609C,NM_000138.5(FBN1):c.1825C>T,0.00001
TP53,p.R273H,NM_000546.6(TP53):c.817C>T,0.0001
BRCA1,p.C61G,NM_007294.4(BRCA1):c.181T>G,0.0
```

### **Required Columns:**
- **`gene`**: Gene symbol (e.g., COL1A1, FBN1, TP53)
- **`variant`**: Protein variant in p.RefPosAlt format (e.g., p.G1340S)
- **`hgvs`**: HGVS notation with transcript (optional but recommended)
- **`gnomad_freq`**: Population frequency (0.0 to 1.0, use 0.0 for novel variants)

### **Ren's TSV Format (Also Supported)**
```tsv
Clinical significance and condition	Chrpos	Variation Name	AA Chg
Pathogenic	chr15:48408313	NM_000138.5(FBN1):c.5938G>C	p.D1980H
Likely pathogenic	chr17:43082434	NM_007294.4(BRCA1):c.181T>G	p.C61G
```

### **Usage:**
```bash
# CSV format
python3 cascade_batch_processor.py --input variants.csv --output results.tsv

# TSV format (auto-detected)
python3 cascade_batch_processor.py --input variants.tsv --output results.tsv

# With frequency filtering
python3 cascade_batch_processor.py --input variants.csv --output results.tsv --freq-filter 0.001
```

---

## 🧬 **DN BATCH PROCESSOR**

**File:** `nova_dn/csv_batch_processor.py`  
**Purpose:** DN-only analysis (faster, mechanism-specific)

### **Format:**
```csv
gene,variant,gnomad_freq
TP53,p.R273H,0.0001
COL1A1,p.G1340S,0.00064
FBN1,p.R609C,0.00001
```

### **Required Columns:**
- **`gene`**: Gene symbol
- **`variant`**: Protein variant (p.RefPosAlt format)
- **`gnomad_freq`**: Population frequency

### **Usage:**
```bash
python3 -m nova_dn.csv_batch_processor --input variants.csv --output dn_results.tsv
```

---

## 📋 **COLUMN SPECIFICATIONS**

### **Gene Column (`gene`)**
- **Format**: Gene symbol (HGNC approved)
- **Examples**: `TP53`, `COL1A1`, `FBN1`, `BRCA1`
- **Case**: Usually uppercase, but system is case-insensitive
- **Required**: ✅ YES

### **Variant Column (`variant`)**
- **Format**: `p.RefPosAlt` (protein HGVS)
- **Examples**: `p.R273H`, `p.G1340S`, `p.C61G`
- **Prefix**: Must include `p.` prefix
- **Required**: ✅ YES

### **HGVS Column (`hgvs`)**
- **Format**: Full HGVS with transcript
- **Examples**: `NM_000546.6(TP53):c.817C>T`, `NM_000088.4(COL1A1):c.4018G>A`
- **Required**: ❌ Optional (but recommended for traceability)

### **Frequency Column (`gnomad_freq`)**
- **Format**: Decimal frequency (0.0 to 1.0)
- **Examples**: `0.00064` (0.064%), `0.0` (novel), `0.25` (25%)
- **Novel variants**: Use `0.0`
- **Required**: ✅ YES (for frequency filtering)

---

## 🔥 **OUTPUT FORMATS**

### **CASCADE Output (TSV)**
```tsv
gene	variant	hgvs	gnomad_freq	status	routing_strategy	final_score	final_classification	explanation
COL1A1	p.G1340S	NM_000088.4(COL1A1):c.4018G>A	0.00064	SUCCESS	LOF_PRIMARY	0.254	VUS	LOF analysis with domain downweighting
TP53	p.R273H	NM_000546.6(TP53):c.817C>T	0.0001	SUCCESS	MULTI_ANALYZER	0.847	LP	DN + LOF synergy in tumor suppressor
```

### **DN Output (TSV)**
```tsv
gene	variant	gnomad_freq	status	dn_likelihood	interface_poisoning	active_site_jamming	lattice_disruption	trafficking_maturation	top_mechanism	explanation
TP53	p.R273H	0.0001	SUCCESS	0.85	0.129	0.55	0.0	0.0	active_site_jamming	Active site disruption in DNA-binding domain
```

---

## 🚨 **COMMON ISSUES & SOLUTIONS**

### **Issue: "Could not extract gene from variation name"**
**Solution**: Ensure HGVS format includes gene in parentheses: `NM_000546.6(TP53):c.817C>T`

### **Issue: "Variant parsing failed"**
**Solution**: Use proper p. prefix: `p.R273H` not `R273H`

### **Issue: "File format not recognized"**
**Solution**: Check column headers match exactly (case-sensitive)

### **Issue: "Frequency filtering too aggressive"**
**Solution**: Adjust `--freq-filter` threshold or use `0.0` for novel variants

---

## 💡 **BEST PRACTICES**

1. **Use CASCADE processor** for comprehensive analysis
2. **Include HGVS notation** for better traceability
3. **Set frequency to 0.0** for novel/rare variants
4. **Use TSV format** for better Excel compatibility
5. **Include ClinVar classification** in filename for automatic comparison
6. **Test with small batches** before running large datasets

---

## 🎯 **EXAMPLES BY USE CASE**

### **ClinVar Export Processing**
```bash
# Filename indicates expected classification
python3 cascade_batch_processor.py --input clinvar_pathogenic.csv --output results.tsv
```

### **Novel Variant Analysis**
```csv
gene,variant,hgvs,gnomad_freq
MYBPC3,p.R502W,NM_000256.3(MYBPC3):c.1504C>T,0.0
TNNT2,p.R92Q,NM_001001430.2(TNNT2):c.275G>A,0.0
```

### **Population Study**
```csv
gene,variant,hgvs,gnomad_freq
CFTR,p.F508del,NM_000492.4(CFTR):c.1521_1523delCTT,0.03
HFE,p.C282Y,NM_000410.4(HFE):c.845G>A,0.04
```

---

*Ready to process thousands of variants with revolutionary AI-powered analysis!* 🧬🔥💜
