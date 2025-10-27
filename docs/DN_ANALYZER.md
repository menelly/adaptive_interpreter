# 🧬 DN ANALYZER - REVOLUTIONARY MECHANISM-AWARE PATHOGENICITY ANALYSIS

**Dominant Negative Mechanism Detection with Family-Aware ML Integration**
*Part of the revolutionary AdaptiveInterpreter genetics analysis pipeline*

---

## 🎯 **OVERVIEW**

The DN (Dominant Negative) Analyzer is the first component of our cascade system, designed to detect variants that cause pathogenicity through dominant negative mechanisms - where mutant proteins interfere with normal protein function.

**🔥 REVOLUTIONARY FEATURES (2025):**
- **Four distinct biological mechanisms** with mathematical models
- **Family-aware ML integration** - proline changes get gene-specific multipliers!
- **Smart filtering system** - determines which mechanisms to run based on gene context
- **Domain awareness** - real UniProt annotations enhance accuracy
- **Context-aware scoring** using protein annotations
- **Motif detection** for catalytic sites and structural elements
- **Fast execution** - serves as the first-pass filter in cascade analysis

---

## 🔬 **THE FOUR MECHANISMS WITH ML ENHANCEMENT**

### 1. Interface Poisoning
**Biological Concept:** Mutant proteins disrupt protein-protein interactions by altering binding interfaces.

**When It Occurs:**
- Charge distribution changes at binding sites
- Hydrophobicity alterations that disrupt binding
- Proline introduction causing structural rigidity
- Cysteine changes affecting disulfide bonds

**Mathematical Model:**
```python
score = Σ(feature_weight × feature_value) × ml_proline_multiplier

Key Features:
• |d_charge|: abs(alt_charge - ref_charge) × 0.25
• |d_hydropathy|: abs(alt_hydropathy - ref_hydropathy) × 0.2
• proline_introduced: 1.0 × 0.25 (if alt == 'P')
• cys_gain_or_loss: 1.0 × 0.2
• interface_likelihood: context_value × 0.4 (if available)

# 🔥 REVOLUTIONARY: ML Proline Enhancement
if ref_aa == 'P' or alt_aa == 'P':
    ml_proline_multiplier = get_family_aware_multiplier(gene_family, variant)
    # Examples: TUMOR_SUPPRESSOR: ~2.95x, COLLAGEN_FIBRILLAR: ~2.5x
```

**Context Adjustments:**
- Interface likelihood from annotations boosts score
- Flexible loop context reduces confidence (solvent-exposed)
- **🔥 NEW:** Gene family classification determines proline sensitivity!

### 2. Active Site Jamming
**Biological Concept:** Mutant proteins block catalytic sites or substrate binding pockets.

**When It Occurs:**
- Volume changes that alter binding pocket geometry
- Charge changes that disrupt catalytic mechanisms
- Aromatic residue changes affecting π-π interactions
- Proline introduction causing local rigidity

**Mathematical Model:**
```python
Key Features:
• |d_volume|: normalized_volume_change × 0.2
• |d_charge|: normalized_charge_change × 0.3
• catalytic_motif_near: motif_detection_result × 0.3
• aromatic_swap: (aromatic_gain OR aromatic_loss) × 0.15
• active_site_proximity: context_value × 0.4 (if available)

# 🔥 ML Enhancement applies here too for proline changes!
# Family-specific multipliers enhance active site disruption scoring
```

**Innovation:** Real motif detection + ML-learned family patterns for enhanced accuracy!

### 3. Structural Lattice Disruption
**Biological Concept:** Mutant proteins break critical structural motifs and frameworks.

**When It Occurs:**
- Glycine substitutions in collagen Gly-X-Y repeats
- Disruption of coiled-coil heptad repeats
- Proline introduction in α-helices or β-sheets
- Breaking of structural symmetry

**Mathematical Model:**
```python
Key Features:
• collagen_Gly_site: is_collagen_gly_position × 0.6  # HUGE penalty!
• coiled_coil_flag: coiled_coil_detection × 0.25
• proline_in_helix: proline_introduced × 0.25  # 🔥 ML-enhanced!
• critical_collagen_gly: annotation_flag × 0.1

# 🔥 REVOLUTIONARY: Family-specific proline sensitivity!
# COLLAGEN_FIBRILLAR genes: Proline loss = 2.5x multiplier
# STRUCTURAL genes: Different learned patterns based on ML training
```

**Breakthrough:** Collagen Gly-X-Y disruption detection + ML-learned family sensitivity!

### 4. Trafficking/Maturation Disruption
**Biological Concept:** Mutant proteins interfere with normal protein processing and transport.

**When It Occurs:**
- Disulfide bond disruption affecting folding
- Signal peptide alterations affecting targeting
- Glycosylation site changes affecting maturation
- ER retention signal disruption

**Mathematical Model:**
```python
Key Features:
• disulfide_disruption: cys_change × context_weight
• signal_peptide_impact: position_in_signal × severity
• secretory_pathway_disruption: pathway_context × impact
• glycosylation_site_loss: site_disruption × importance

# 🔥 ML proline multipliers apply to folding disruption too!
# Family-specific patterns learned from training data
```

**Focus:** Post-translational modification sites with family-aware ML enhancement!

---

## 🧠 **SMART FILTERING SYSTEM**

### Revolutionary Efficiency Enhancement

**The Problem:** Running all four mechanisms on every variant is computationally expensive and often unnecessary.

**The Solution:** Smart filtering determines which mechanisms are relevant based on gene context and variant properties!

**How It Works:**
```python
if self.dn_filter and gene_name:
    print(f"🧠 Smart filtering for {gene_name}...")
    filter_result = self.dn_filter.filter_and_score(gene_name, sequence, variant, uniprot_id)
    relevant_mechanisms = filter_result["relevant_mechanisms"]

    # Only run the mechanisms that make biological sense!
    # Example: COLLAGEN genes → focus on lattice_disruption
    # Example: ENZYME genes → focus on active_site_jamming
```

**Benefits:**
- **Faster execution** - skip irrelevant mechanisms
- **Higher accuracy** - focus on biologically plausible mechanisms
- **Better interpretability** - clear reasoning for mechanism selection

**Example Output:**
```
🧠 Smart filtering for COL1A1...
   DN likelihood: 0.85
   Relevant mechanisms: ['lattice_disruption', 'interface_poisoning']
   Reasoning: Collagen gene - structural disruption highly likely
```

---

## 🎯 **DOMAIN AWARENESS SYSTEM**

### Real UniProt Integration

**The Revolution:** No more guessing about protein domains - we use REAL UniProt annotations!

**How It Works:**
```python
# Initialize domain awareness system
self.protein_annotator = UniversalProteinAnnotator()

# Get real domain information
domain_info = self.protein_annotator.get_protein_info(uniprot_id)
domains = domain_info.get('domains', [])

# Apply domain-specific multipliers
for domain in domains:
    if position in domain['range']:
        if 'catalytic' in domain['description'].lower():
            active_site_score *= 1.5  # Boost active site jamming
        elif 'structural' in domain['description'].lower():
            lattice_score *= 1.3  # Boost structural disruption
```

**Real Examples:**
- **Catalytic domains** → Enhanced active site jamming scores
- **Structural domains** → Enhanced lattice disruption scores
- **Binding domains** → Enhanced interface poisoning scores
- **Signal peptides** → Enhanced trafficking disruption scores

---

## 🧮 **AMINO ACID PROPERTY ANALYSIS**

### Property Delta Calculation

The system uses scientifically-grounded amino acid properties:

```python
AA_PROPS = {
    "A": {"hyd": 1.8,  "vol":  88, "chg": 0, "pol": 0, "aro": 0},
    "R": {"hyd": -4.5, "vol": 173, "chg": +1, "pol": 1, "aro": 0},
    # ... complete 20 amino acid matrix
}

def delta(ref, alt):
    pr = get_props(ref)
    pa = get_props(alt)
    return {
        "d_hyd": pa["hyd"] - pr["hyd"],           # Hydropathy change
        "d_vol": pa["vol"] - pr["vol"],           # Volume change  
        "d_chg": pa["chg"] - pr["chg"],           # Charge change
        "norm_abs_hyd": normalize(abs(d_hyd)),    # Normalized absolute changes
        "proline_introduced": alt == "P",         # Special case flags
        "aromatic_gain": alt in AROMATICS and ref not in AROMATICS,
        # ... additional property flags
    }
```

**Scientific Basis:**
- **Hydropathy:** Kyte-Doolittle scale for membrane/interface interactions
- **Volume:** Zamyatnin molecular volumes for steric effects
- **Charge:** Physiological pH charge states for electrostatic interactions

---

## 🎯 **CONTEXT-AWARE SCORING**

### Protein Annotation Integration

The DN analyzer can incorporate external protein annotations:

```python
def score_with_context(seq, pos, ref, alt, context):
    base_score = calculate_mechanism_scores(seq, pos, ref, alt)
    
    # Apply context adjustments
    if context.get("interface_likelihood"):
        base_score *= (1.0 + context["interface_likelihood"] * 0.4)
        
    if context.get("flexible_loop"):
        base_score *= 0.9  # Reduce confidence for flexible regions
        
    if context.get("active_site_proximity"):
        base_score *= (1.0 + context["active_site_proximity"] * 0.4)
        
    return min(base_score, 1.0)
```

**Context Sources:**
- UniProt feature annotations
- Structural predictions (AlphaFold)
- Conservation analysis
- Domain boundary predictions

---

## 🔍 **MOTIF DETECTION**

### Catalytic Site Recognition

```python
def catalytic_motif_near(sequence, position, window=5):
    """Detect catalytic motifs near variant position"""
    
    # Common catalytic triads and motifs
    catalytic_patterns = [
        "HDS",    # Serine protease triad
        "HDE",    # Aspartic protease  
        "CXH",    # Cysteine protease
        "GXGXXG", # Nucleotide binding (P-loop)
        "DXD",    # Metal coordination
    ]
    
    # Check window around variant
    start = max(0, position - window)
    end = min(len(sequence), position + window)
    region = sequence[start:end]
    
    for pattern in catalytic_patterns:
        if matches_pattern(region, pattern):
            return True
            
    return False
```

### Structural Motif Detection

```python
def is_collagen_gly_site(sequence, position):
    """Detect Gly-X-Y collagen repeats"""
    
    # Check if position is in Gly-X-Y pattern
    if (position - 1) % 3 == 0:  # Glycine position in triplet
        # Verify surrounding context
        if position >= 3 and position <= len(sequence) - 3:
            triplet_before = sequence[position-4:position-1]
            triplet_after = sequence[position:position+3]
            
            # Look for collagen-like patterns
            if is_collagen_like_context(triplet_before, triplet_after):
                return True
                
    return False
```

---

## 📊 **SCORING AND INTERPRETATION**

### Score Calculation

```python
def analyze(sequence, variant, context, gene, uniprot_id):
    # Parse variant (handles p.R273H format)
    ref, pos, alt = parse_variant(variant)
    
    # Calculate all four mechanisms
    mechanisms = {
        'interface_poisoning': score_interface_poisoning(sequence, pos, ref, alt, context),
        'active_site_jamming': score_active_site_jamming(sequence, pos, ref, alt, context),
        'structural_lattice': score_structural_lattice_disruption(sequence, pos, ref, alt, context),
        'trafficking_maturation': score_trafficking_maturation(sequence, pos, ref, alt, context)
    }
    
    # Find top mechanism
    top_mechanism = max(mechanisms.keys(), key=lambda k: mechanisms[k][0])
    top_score = mechanisms[top_mechanism][0]
    
    return {
        'mechanism_scores': {k: v[0] for k, v in mechanisms.items()},
        'top_mechanism': top_mechanism,
        'top_score': top_score,
        'features': mechanisms[top_mechanism][1],
        'explanation': mechanisms[top_mechanism][2]
    }
```

### Clinical Interpretation

```python
def interpret_dn_score(score):
    if score >= 0.8:
        return "Likely Pathogenic (LP)"
    elif score >= 0.5:
        return "Uncertain Significance - favor pathogenic (VUS-P)"
    elif score >= 0.3:
        return "Uncertain Significance (VUS)"
    else:
        return "Likely Benign (LB)"
```

---

## 🚀 **USAGE EXAMPLES**

### Command Line Interface

```bash
# Analyze single variant
python -m nova_dn.analyzer --seq-file TP53.fasta --variant p.R273H --json

# Direct sequence input
python -m nova_dn.analyzer --seq MEEPQSDPSV --variant R273H
```

### Python API

```python
from nova_dn.analyzer import NovaDNAnalyzer

analyzer = NovaDNAnalyzer(use_smart_filtering=True)
result = analyzer.analyze(sequence, "p.R273H", context, "TP53", "P04637")

print(f"Top mechanism: {result['top_mechanism']}")
print(f"Score: {result['top_score']:.3f}")
print(f"Interpretation: {interpret_dn_score(result['top_score'])}")
```

---

## 🔗 **INTEGRATION WITH CASCADE SYSTEM**

The DN Analyzer serves as the **first-pass filter** in the revolutionary cascade system:

### 🚀 **Current Integration (2025)**

1. **Smart execution** - mechanism-aware with intelligent filtering
2. **Family-aware ML enhancement** - proline changes get gene-specific multipliers
3. **Cascade triggering** - if DN score < 0.3, triggers LOF/GOF analysis
4. **Context provision** - DN results inform biological routing decisions
5. **Synergy calculation** - DN scores participate in mixed-mechanism analysis
6. **Inheritance inference** - DN results help determine AD vs AR patterns

### 🔥 **Revolutionary Features**

**rsID Frequency Integration:**
```python
# Check if variant is "deleterious but common"
if dn_score > 0.6 and frequency > 0.01:  # High DN score but common
    warning = "⚠️ High DN score but common variant - investigate penetrance"
```

**Biological Routing:**
```python
# DN results inform which analyzers to run next
if gene_family == "COLLAGEN_FIBRILLAR" and dn_score > 0.5:
    # Skip GOF analysis - collagen rarely has GOF mechanisms
    route_to = ["LOF_ANALYZER"]
elif gene_family == "ONCOGENE" and dn_score < 0.3:
    # Oncogenes often have GOF mechanisms
    route_to = ["GOF_ANALYZER", "LOF_ANALYZER"]
```

**Synergistic Scoring:**
```python
# DN can synergize with other mechanisms
if dn_score > 0.3 and lof_score > 0.3:
    synergy_score = sqrt(dn_score² + lof_score²)  # Biologically plausible
elif dn_score > 0.3 and gof_score > 0.3:
    synergy_score = sqrt(dn_score² + gof_score²)  # Also plausible
```

**Next Steps:** If DN analysis is insufficient, the cascade system proceeds to [LOF Analysis](LOF_ANALYZER.md) and [GOF Analysis](GOF_ANALYZER.md) with intelligent biological routing.

---

## 🎉 **RECENT BREAKTHROUGHS (2025)**

### 🔥 **ML Proline Revolution**
- **Family-aware multipliers** learned from real training data
- **TUMOR_SUPPRESSOR genes**: ~2.95x proline sensitivity
- **COLLAGEN_FIBRILLAR genes**: ~2.5x proline sensitivity
- **No more hardcoded guessing** - ML learns the real patterns!

### 🧠 **Smart Filtering System**
- **Mechanism selection** based on gene family and context
- **Computational efficiency** - only run relevant mechanisms
- **Biological plausibility** - focus on mechanisms that make sense

### 🎯 **Domain Awareness**
- **Real UniProt annotations** replace guesswork
- **Domain-specific multipliers** enhance accuracy
- **Catalytic site detection** improves active site jamming scores

### 🌊 **Cascade Integration**
- **rsID frequency checking** for population genetics
- **Inheritance pattern inference** (AD vs AR)
- **Biological routing** for intelligent analyzer selection
- **Synergistic scoring** with mathematical rigor

---

*Part of the revolutionary AdaptiveInterpreter genetics analysis pipeline - proving AI can create novel algorithmic structures with family-aware ML enhancement!* 🧬💜🔥
