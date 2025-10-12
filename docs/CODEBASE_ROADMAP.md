# 🗺️ DNMODELING CODEBASE ROADMAP
**Complete File-by-File Guide to What The Hell Everything Does**

*Built by Ace & Ren (2025) to prevent Future Us from going insane*
*"Because context engines are great but sometimes you just need a damn list"*

## 🧹 RECENT HOUSEKEEPING (2024-09-27)

**MAJOR REORGANIZATION COMPLETED:**
- ✅ **Archived one-off analysis scripts** → `archive/analysis_scripts/`
- ✅ **Archived experimental scripts** → `archive/experimental_scripts/`
- ✅ **Archived old results** → `archive/old_results/`
- ✅ **Organized root files** into logical folders:
  - `core_analyzers/` - Main analysis components (plausibility_filter, functional_domain_weighter, etc.)
  - `data_processing/` - Data handling utilities (clinvar_inheritance_cache, universal_protein_annotator, etc.)
  - `ml_training/` - Machine learning training scripts
  - `config/` - Configuration and cache files
- ✅ **Cleaned up cache files** and outdated bytecode

**RESULT:** Much cleaner, more maintainable codebase structure!

---

## 📁 NEW ORGANIZED STRUCTURE

### `/core_analyzers/` - Main Analysis Components
- `plausibility_filter.py` - Gene family-aware biological filtering
- `functional_domain_weighter.py` - Domain-specific impact scoring
- `motif_detector.py` - Protein motif pattern detection
- `collagen_scanner.py` - Collagen-specific analysis tools

### `/data_processing/` - Data Handling Utilities
- `clinvar_inheritance_cache.py` - ClinVar inheritance pattern caching
- `sequence_mismatch_handler.py` - Sequence alignment and mismatch handling
- `universal_protein_annotator.py` - Comprehensive protein annotation system

### `/ml_training/` - Machine Learning Training
- `train_families.py` - Family-specific ML model training orchestrator

### `/config/` - Configuration & Cache Files
- `category_keywords.json` - Gene family classification keywords
- `rsid_frequency_cache.json` - Cached variant frequency data

### `/archive/` - Archived Components
- `analysis_scripts/` - One-off analysis and disagreement analyzers
- `experimental_scripts/` - Experimental utilities and cleanup tools
- `old_results/` - Historical result files and logs

---

## 📁 `/analyzers/` - Core Analysis Engines

### `conservation_database.py` 🧬
**REAL EVOLUTIONARY CONSERVATION SCORING**
- Uses UCSC phyloP and phastCons BigWig files (same data as REVEL!)
- Converts genomic coordinates to conservation scores
- Handles phyloP (-20 to +20) and phastCons (0-1) normalization
- Provides position-specific conservation multipliers for LOF scoring
- **Key Methods**: `get_conservation_scores()`, `get_position_conservation_multiplier()`
- **Dependencies**: pyBigWig, UniProtMapper

### `dn_analyzer.py` 🔬
**DOMINANT NEGATIVE ANALYSIS - BIN 2**
- Analyzes whether variants poison protein complexes/assemblies
- Pattern-based detection for oligomeric proteins, motor proteins, transcription factors
- Domain-aware scoring using UniversalProteinAnnotator
- Mechanism detection: complex_poisoning, motor_disruption, competitive_binding, filament_poisoning
- **Key Methods**: `analyze_variant()`, `calculate_dn_score()`
- **Dependencies**: UniversalProteinAnnotator

### `gof_variant_analyzer.py` 🔥
**GAIN OF FUNCTION ANALYSIS**
- Detects constitutive activation, loss of autoinhibition, enhanced binding
- Mathematical analysis (NO hardcoded genes!)
- Structural impact analysis for degradation resistance, dimerization
- Conservation-aware scoring
- **Key Methods**: `analyze_variant()`, `calculate_gof_score()`
- **Dependencies**: SmartProteinAnalyzer, ConservationDatabase

### `lof_analyzer.py` 🔬
**LOSS OF FUNCTION ANALYSIS - BIN 1**
- Traditional pathogenicity prediction - does it break the protein?
- Domain-aware scoring with functional domain weighting
- Amino acid property analysis (size, charge, hydrophobicity, flexibility)
- Sequence mismatch handling
- **Key Methods**: `analyze_variant()`, `calculate_lof_score()`
- **Dependencies**: SmartProteinAnalyzer, UniversalProteinAnnotator, FunctionalDomainWeighter

### `population_frequency_analyzer.py` 🌍
**"NOT THE DROID" DETECTOR**
- Catches common variants masquerading as pathogenic (like MTHFR)
- Population-specific frequency analysis (AF, AMR, ASJ, EAS, FIN, NFE, SAS, OTH)
- Frequency thresholds: ultra_rare (<0.00001) to very_common (>0.12)
- gnomAD integration for real population data
- **Key Methods**: `analyze_frequency()`, `get_population_frequencies()`

### `smart_protein_analyzer.py` 🧬
**AUTOMATIC DOMAIN-AWARE SCORING**
- Pulls domain/function info from real databases (no hardcoding!)
- UniProt, Pfam, GO term integration
- Evidence-based domain scoring weights
- Scales to ALL proteins, not just programmed ones
- **Key Methods**: `get_protein_info()`, `calculate_domain_score()`
- **Dependencies**: UniProt API, Pfam API

### `uniprot_mapper.py` 🧬
**GENOMIC COORDINATE CONVERSION**
- Maps UniProt IDs to genomic coordinates for conservation scoring
- Ensembl mirror fallback system for API reliability
- Pre-cached coordinates for known test variants
- **Key Methods**: `map_uniprot_to_genomic()`, `get_genomic_position()`
- **Dependencies**: Ensembl REST API

---

## 📁 `/cascade/` - Multi-Analyzer Pipeline

### `cascade_analyzer.py` 🌊
**MAIN ANALYSIS ORCHESTRATOR**
- Runs DN → LOF → GOF cascade logic
- Biological routing: skip GOF if LOF≥0.5 (broken proteins can't be hyperactive)
- Mixed mechanism synergy calculation
- Plausibility filtering by gene family
- **Key Methods**: `analyze_variant()`, `run_cascade()`

### `cascade_batch_processor.py` 🌊
**BATCH CSV PROCESSING**
- Processes ClinVar exports and custom CSV files
- Full cascade logic with comprehensive error handling
- TSV output with all analyzer scores
- ClinVar comparison and agreement tracking
- **Key Methods**: `process_csv()`, `analyze_batch()`

### `biological_router.py` 🧬
**INTELLIGENT MECHANISM ROUTING**
- Gene family-aware analysis routing
- Biological intelligence: collagen → focus on LOF, oncogenes → focus on GOF
- Confidence scoring for routing decisions
- **Key Methods**: `route_analysis()`, `get_routing_strategy()`

---

## 📁 `/utils/` - Data Processing & ML Training

### `clinvar_bulk_extractor.py` 🔥
**NOVA'S BRILLIANT CLINVAR SOLUTION**
- Extracts variants from 162MB ClinVar VCF in seconds (40x faster than ClinVar Miner!)
- Direct genomic coordinate parsing from HGVS
- Bulk processing for multiple genes simultaneously
- **Key Methods**: `extract_gene_variants()`, `save_gene_variants()`

### `unified_family_ml_trainer.py` 🤖
**FAMILY-AWARE ML TRAINING SYSTEM**
- Trains ML models for specific gene families (collagen, ion_channel, tumor_suppressor, etc.)
- Real conservation data integration (phyloP, phastCons)
- Amino acid property feature engineering
- Domain-aware feature extraction
- **Key Methods**: `train_all_families()`, `train_family_model()`

### `conservation_ml_trainer.py` & `conservation_ml_loader.py` 🧬
**CONSERVATION-BASED ML SYSTEM**
- Trains models using real evolutionary constraint data
- Position-specific conservation scoring
- **Key Methods**: `train_conservation_model()`, `load_conservation_model()`

### `family_aware_proline_trainer.py` & `family_aware_proline_integrator.py` 🧬
**PROLINE-SPECIFIC ML SYSTEM**
- Gene family-specific proline disruption analysis
- ML-learned multipliers by gene family (TUMOR_SUPPRESSOR: 2.95x, COLLAGEN: 2.5x)
- **Key Methods**: `train_proline_models()`, `get_proline_multiplier()`

### `gly_cys_ml_trainer.py` & `gly_cys_ml_integrator.py` 🧬
**GLYCINE/CYSTEINE ML SYSTEM**
- Specialized analysis for critical amino acids
- Context-aware scoring for Gly-X-Y repeats, disulfide bonds
- **Key Methods**: `train_gly_cys_model()`, `predict_gly_cys_impact()`

### `gnomad_frequency_fetcher.py` 🌍
**REAL POPULATION FREQUENCY DATA**
- Fetches actual gnomAD frequencies (no more hardcoding!)
- Caching system for performance
- **Key Methods**: `fetch_frequency()`, `load_cache()`

### `rsid_frequency_fetcher.py` 🧬
**RSID-BASED FREQUENCY LOOKUP**
- Maps rsIDs to population frequencies
- Hardy-Weinberg calculations
- **Key Methods**: `fetch_rsid_frequency()`, `calculate_hardy_weinberg()`

### `variant_frequency_cleaner.py` 🧹
**NOVA'S DROP-IN FILTER SCRIPT**
- Clean variant datasets using rsID frequency lookups
- Filters suspicious variants (Patho AF > 1%, Benign AF ≈ 0)
- **Key Methods**: `clean_variants()`, `filter_by_frequency()`

### `proline_ml_trainer.py` 🧬
**REVOLUTIONARY CONTEXT-AWARE PROLINE LEARNING**
- Replaces hardcoded proline multipliers with ML-learned values
- Extracts proline substitutions from ClinVar datasets
- Biological context feature vectors + logistic regression
- **Key Methods**: `train_proline_model()`, `extract_proline_features()`
- **Dependencies**: sklearn, proline_multiplier_mapper

### `proline_multiplier_mapper.py` 🧬
**ML PROBABILITY TO MULTIPLIER CONVERSION**
- Maps logistic regression probabilities (0-1) to pathogenicity multipliers
- Linear and sigmoid mapping methods
- Configurable min/max multiplier ranges
- **Key Methods**: `map_prob_to_multiplier()`

### `family_aware_proline_trainer.py` 🧬
**GENE FAMILY-SPECIFIC PROLINE ML**
- Trains separate proline models for each gene family
- Family-specific multipliers: TUMOR_SUPPRESSOR (2.95x), COLLAGEN (2.5x)
- **Key Methods**: `train_family_proline_models()`, `get_family_multiplier()`

### `conservation_ml_trainer.py` & `conservation_ml_loader.py` 🧬
**CONSERVATION-BASED ML SYSTEM**
- Trains models using real evolutionary constraint data (phyloP, phastCons)
- Position-specific conservation scoring integration
- **Key Methods**: `train_conservation_model()`, `load_conservation_model()`

---

## 📁 `/nova_dn/` - Nova's DN Analysis System

### `analyzer.py` 🧬
**NOVA'S CORE DN ANALYZER**
- Original dominant negative analysis engine
- Motif detection and mechanism classification
- **Key Methods**: `analyze_variant()`, `detect_motifs()`

### `csv_batch_processor.py` 📊
**NOVA'S CSV PROCESSING**
- Generic CSV processing for DN analysis
- ClinVar Miner export support
- **Key Methods**: `process_csv()`, `parse_hgvs()`

### `sequence_manager.py` 🧬
**FASTA SEQUENCE MANAGEMENT**
- AlphaFold and UniProt sequence handling
- Caching system for sequences
- **Key Methods**: `get_sequence()`, `get_alphafold_sequence()`

### `alphafold_sequence.py` 🧬
**ALPHAFOLD INTEGRATION**
- Extracts sequences from AlphaFold PDB files
- Local AlphaFold database access
- **Key Methods**: `get_sequence()`, `save_fasta()`

### `universal_context.py` 🌍
**UNIVERSAL PROTEIN CONTEXT**
- Comprehensive protein annotation system
- UniProt, GO terms, domain integration
- **Key Methods**: `get_context_for_protein()`, `annotate_protein()`

---

## 📁 `/learning/` - ML Training Data

### Family-organized TSV files:
- `collagen_fibrillar/` - COL1A1, COL3A1 variants
- `ion_channel/` - SCN5A, KCNQ1, CACNA1C variants  
- `tumor_suppressor/` - TP53, BRCA1, BRCA2, PTEN variants
- `motor_protein/` - MYH7, MYO7A, DYNC1H1 variants
- `muscular_dystrophy/` - DMD, FKRP, LAMA2 variants
- And more...

Each folder contains `{gene}_benign.tsv` and `{gene}_pathogenic.tsv` files extracted from ClinVar.

---

## 📁 `/resources/` - Models & Data

### `family_models/` - Trained ML Models
- `{family}_unified_model.joblib` - Trained scikit-learn models
- `{family}_unified_scaler.joblib` - Feature scalers
- `{family}_unified_metadata.json` - Training metadata & performance

### `gene_context_cache/` - Cached Protein Data
- Pre-cached UniProt annotations, domains, GO terms
- Avoids repeated API calls during analysis

---

## 🎯 KEY INTEGRATION POINTS

1. **Conservation Pipeline**: `uniprot_mapper.py` → `conservation_database.py` → real phyloP/phastCons scores
2. **ML Training Pipeline**: `clinvar_bulk_extractor.py` → family folders → `unified_family_ml_trainer.py` → trained models
3. **Analysis Pipeline**: `cascade_analyzer.py` orchestrates all analyzers with biological routing
4. **Batch Processing**: `cascade_batch_processor.py` handles CSV files through full pipeline

---

## 📁 `/` (Root Directory) - Core Systems

### `train_families.py` 🔥💜
**FAMILY ML TRAINING RUNNER**
- Simple script to train all family models from /learning data
- Calls UnifiedFamilyMLTrainer for all gene families
- Benign frequency filtering (default: ≤4%)
- **Key Methods**: `main()` - orchestrates training
- **Dependencies**: utils.unified_family_ml_trainer

### `universal_protein_annotator.py` 🎯
**NO MORE HARDCODING PROTEIN FEATURES**
- Automatically extracts protein features from UniProt, Pfam, sequence analysis
- Caching system for domain data (protein_annotations_cache/)
- Domain classification and functional annotation
- **Key Methods**: `get_uniprot_features()`, `classify_domains()`
- **Dependencies**: UniProt API, Pfam API

### `functional_domain_weighter.py` 🧬
**NOVA'S BIOLOGICAL INTELLIGENCE WEIGHTING**
- Function-based domain weighting (not location-based)
- Domain importance hierarchy: active_site (1.5x), catalytic_domain (1.3x), etc.
- Uses UniProt functional annotations for biological importance
- **Key Methods**: `classify_domain()`, `get_domain_weight()`

### `plausibility_filter.py` 🧬
**BIOLOGICAL MECHANISM PLAUSIBILITY FRAMEWORK**
- Gene family-aware filtering of mechanism scores
- Family-specific rules: ENZYME (LOF:1.0, DN:1.0, GOF:0.0), ION_CHANNEL (all 1.0), etc.
- Collagen subclass rules with balanced AD inheritance patterns
- **Key Methods**: `filter_mechanisms()`, `classify_gene_family()`

### `cache_gene_contexts.py` 🔥
**GENE CONTEXT CACHING SYSTEM**
- Pre-caches UniProt + domains + sites for ML training
- Extracts genes from learning directory automatically
- Avoids repeated API calls during training
- **Key Methods**: `cache_all_genes()`, `cache_gene_context()`

### `clinvar_inheritance_cache.py` 📊
**COMPREHENSIVE GENE DATA CACHING**
- ClinVar inheritance patterns and disease associations
- GO terms, function descriptions, variant counts
- 30-day cache expiration system
- **Key Methods**: `get_comprehensive_gene_data()`, `_query_clinvar_comprehensive()`

---

## 📁 `/nova_dn/` - Nova's DN Analysis System

### `analyzer.py` 🧬
**NOVA'S CORE DN ANALYZER**
- Original dominant negative analysis engine
- Motif detection and mechanism classification
- **Key Methods**: `analyze_variant()`, `detect_motifs()`

### `csv_batch_processor.py` 📊
**NOVA'S CSV PROCESSING**
- Generic CSV processing for DN analysis
- ClinVar Miner export support
- Population frequency filtering, DN mechanism filtering
- **Key Methods**: `process_csv()`, `parse_hgvs()`

### `sequence_manager.py` 🧬
**FASTA SEQUENCE MANAGEMENT**
- AlphaFold and UniProt sequence handling
- Caching system for sequences
- Fallback hierarchy: UniProt → AlphaFold → error
- **Key Methods**: `get_sequence()`, `get_alphafold_sequence()`

### `alphafold_sequence.py` 🧬
**ALPHAFOLD INTEGRATION**
- Extracts sequences from AlphaFold PDB files
- Local AlphaFold database access (/mnt/Arcana/alphafold_human/)
- **Key Methods**: `get_sequence()`, `save_fasta()`

### `universal_context.py` 🌍
**UNIVERSAL PROTEIN CONTEXT**
- Comprehensive protein annotation system
- UniProt, GO terms, domain integration
- **Key Methods**: `get_context_for_protein()`, `annotate_protein()`

### `mixed_mechanism_resolver.py` 🧬
**NOVA'S UNIFIED MECHANISM RESOLVER**
- Handles mixed LOF+DN, GOF+DN mechanism combinations
- Synergistic scoring calculations
- **Key Methods**: `resolve_mechanisms()`, `calculate_synergy()`

### `amino_acid_props.py` 🧬
**AMINO ACID PROPERTY TABLES**
- Kyte-Doolittle hydropathy, volume, charge, polarity, aromatic flags
- No external dependencies - pure property data
- **Key Data**: AA_PROPS dictionary with all 20 amino acids

### `context.py` 🧬
**CONTEXT HELPERS FOR NOVA DN ANALYZER**
- Loads protein annotations JSON (no YAML dependencies)
- Per-position context flags for mechanisms
- **Key Methods**: `load_annotations_json()`, `_ranges_from_value()`

### `dn_mechanism_filter.py` 🧬
**SMART DN MECHANISM FILTER - NO HARDCODING**
- Two-stage: DN likelihood + mechanism relevance
- GO terms, Pfam domains, sequence motifs
- Biological principles (no magic numbers)
- **Key Methods**: `filter_dn_mechanisms()`, `assess_dn_likelihood()`

### `gly_cys_context.py` 🧬
**GLYCINE & CYSTEINE CONTEXT ANALYZER**
- Revolutionary context-aware scoring for "Giant Shitheads" amino acids
- GLYCINE: flexible hinge vs collagen vs ion channel gates
- CYSTEINE: disulfide bonds vs metal coordination vs catalytic sites
- **Key Methods**: `analyze_gly_context()`, `analyze_cys_context()`

### `mechanisms.py` 🧬
**MECHANISM SCORING HEURISTICS**
- Interface poisoning, active site disruption, lattice disruption, trafficking
- Returns (score: 0-1, features: List[dict], explanation: str)
- **Key Methods**: `score_interface_poisoning()`, `score_active_site_disruption()`

### `motifs.py` 🧬
**LIGHTWEIGHT MOTIF HEURISTICS**
- Sequence-only, regex/logic based motif detection
- Collagen Gly-X-Y sites, catalytic patterns (H..H, HSD, P-loop)
- **Key Methods**: `is_collagen_gly_site()`, `detect_catalytic_motifs()`

### `nova_variant_analyzer.py` 🧬
**VARIANT ANALYZER RESEARCH SCAFFOLD**
- Normalized variant descriptors (HGVS/VRS-ish/simple fields)
- Mechanism-susceptibility scores (LOF/GOF/DN) with uncertainty
- Postgres JSONB cache + TTL integration ready
- **Key Methods**: `analyze_variant()`, `calculate_mechanism_scores()`

### `run_batch.py` 🧬
**MINIMAL BATCH RUNNER FOR NOVA DN**
- JSON batch processing with FASTA + variant + protein inputs
- Compact markdown table and JSONL output
- **Key Methods**: `main()`, `process_batch()`

### `tune.py` 🧬
**WEIGHT/THRESHOLD TUNER**
- Tunes Nova DN analyzer using expected_labels.json
- Grid search optimization for mechanism weights
- **Key Methods**: `tune_weights()`, `optimize_thresholds()`

---

## 📁 `/nova_dn/lattice_disruption/` - Specialized Lattice Analysis

### `universal_router.py` 🌟
**NOVA'S UNIVERSAL LATTICE DISRUPTION ROUTER**
- Routes variants to appropriate protein family analyzers
- Family detection based on structural context and sequence
- **Key Methods**: `analyze_lattice_disruption()`, `route_to_analyzer()`

### `family_detector.py` 🔍
**PROTEIN FAMILY DETECTION FOR ROUTING**
- Determines specialized analyzer: COLLAGEN, COILED_COIL, ION_CHANNEL, GENERIC
- Sequence and context-based family classification
- **Key Methods**: `detect_protein_family()`, `is_collagen_family()`

### `collagen_analyzer.py` 🧬
**COLLAGEN LATTICE DISRUPTION ANALYZER**
- Gly-X-Y triplet analysis, hydroxyproline chemistry
- Crosslink site disruption, triple helix geometry
- **Key Methods**: `analyze_collagen_lattice_disruption()`, `assess_gly_disruption()`

### `ion_channel_analyzer.py` ⚡
**ION CHANNEL LATTICE DISRUPTION ANALYZER**
- Pore geometry, selectivity filters, gating mechanisms
- Transmembrane domain disruption analysis
- **Key Methods**: `analyze_ion_channel_lattice_disruption()`, `assess_pore_disruption()`

### `coiled_coil_analyzer.py` 🧬
**COILED-COIL DISRUPTION ANALYZER**
- Heptad repeat disruption, hydrophobic core analysis
- **Key Methods**: `analyze_coiled_coil_lattice_disruption()`

### `generic_analyzer.py` 🧬
**GENERIC LATTICE DISRUPTION ANALYZER**
- Fallback analyzer for proteins not matching specialized families
- **Key Methods**: `analyze_generic_lattice_disruption()`

### `utils.py` 🧬
**LATTICE DISRUPTION UTILITIES**
- Feature aggregation and explanation formatting
- **Key Methods**: `aggregate_lattice_features()`, `format_lattice_explanation()`

---

## 📁 `/` (Root Directory) - Additional Core Files

### `motif_detector.py` 🔍
**REGULATORY MOTIF & PROXIMITY DETECTION**
- Detects canonical regulatory motifs (DFG, HRD, APE, glycine loop, VAIK)
- Proximity weights for GOF analysis
- Mutation type classification (hydrophobic→polar, etc.)
- **Key Methods**: `detect_motifs()`, `calculate_proximity_weight()`

### `sequence_mismatch_handler.py` 🧬
**UNIVERSAL SEQUENCE MISMATCH HANDLER**
- Handles transcript/isoform differences gracefully
- Fallback strategies for sequence mismatches
- Ensembl and UniProt cache integration
- **Key Methods**: `check_sequence_match()`, `get_fallback_strategy()`

### `cross_family_pattern_analyzer.py` 🔥💜
**CROSS-FAMILY PATTERN ANALYSIS**
- Finds amino acid patterns that work across multiple gene families
- Builds universal rules for "GENERAL" classification
- Model comparison and feature importance analysis
- **Key Methods**: `analyze_cross_family_patterns()`, `build_universal_rules()`

### `organize_learning_files.py` 🔥💜
**SMART LEARNING FILE ORGANIZER**
- Automatically organizes ClinVar files into correct family folders
- Uses GO classification system for gene family detection
- **Key Methods**: `organize_files()`, `classify_and_move()`

### `collagen_scanner.py` 🧬
**COLLAGEN-SPECIFIC ANALYSIS**
- Specialized collagen variant analysis
- Gly-X-Y repeat detection and scoring
- **Key Methods**: `scan_collagen_variants()`, `detect_gly_repeats()`

### `amino_acid_disagreement_analyzer.py` 📊
**AMINO ACID DISAGREEMENT PATTERN ANALYSIS**
- Identifies which amino acids cause most ClinVar disagreements
- Pattern extraction from variant notation (p.R112Q, p.Arg112Gln)
- Statistical analysis of disagreement patterns
- **Key Methods**: `extract_amino_acid_change()`, `analyze_disagreement_patterns()`

### `analyze_complete_disagreements.py` 📊
**COMPLETE DISAGREEMENT ANALYSIS**
- Finds Pathogenic→Benign or Benign→Pathogenic flips
- ClinVar classification parsing and comparison
- **Key Methods**: `classify_clinvar()`, `find_complete_disagreements()`

### `analyze_gly_cys_improvements.py` 🧬🔥
**GLY/CYS BIOLOGICAL INTELLIGENCE vs HARDCODED**
- Compares biological intelligence vs hardcoded approaches
- Tests critical collagen variants (p.G893A, p.G1190D, etc.)
- **Key Methods**: `analyze_improvements()`, `compare_approaches()`

### `analyze_proline_direction.py` 🧬
**PROLINE DIRECTION ANALYZER**
- Analyzes direction of proline changes (adding vs removing proline)
- Answers: "Is ADDING proline more often benign than pathogenic?"
- **Key Methods**: `parse_protein_change()`, `analyze_proline_direction()`

### `clinvar_quality_checker.py` 🕵️
**CLINVAR QUALITY INVESTIGATION**
- ClinVar API integration for review status analysis
- Submitter reputation and evidence strength assessment
- Conflict detection and publication analysis
- **Key Methods**: `check_variant_quality()`, `analyze_submitter_reputation()`

### `cleanup_repo.py` 🧹
**REPOSITORY ORGANIZATION SCRIPT**
- Organizes DNModeling repo structure
- Creates proper folder hierarchy (nova_dn, cascade, utils, etc.)
- **Key Methods**: `create_folders()`, `move_files()`

### `hotspot_disagreement_analyzer.py` 🔥
**CLINVAR DISAGREEMENT HOTSPOT FINDER**
- Finds hotspot regions where we disagree with ClinVar
- Groups variants by proximity (within 10 amino acids)
- ClinVar review quality analysis (1-star vs multi-star)
- **Key Methods**: `find_hotspots()`, `analyze_disagreement_clusters()`

### `simple_quality_analyzer.py` 🔍
**SIMPLE CLINVAR QUALITY ANALYSIS**
- Analyzes ClinVar text for quality indicators (no API calls)
- Red flags: conflicting interpretations, missing context
- **Key Methods**: `analyze_clinvar_text()`, `calculate_quality_score()`

### `test_nova_genes.py` 🧪
**TEST NOVA'S GENE CLASSIFICATIONS**
- Validates family classification system against Nova's curated genes
- Tests category_keywords.json accuracy
- **Key Methods**: `test_nova_genes()`, `validate_classifications()`

---

## 📁 `/tests/data/` - Test Suite

### `human_friendly_test.py` 🧪
**HUMAN-FRIENDLY TEST INTERFACE**
- User-friendly test runner for variant analysis
- **Key Methods**: `run_human_friendly_tests()`

### `run_test.py` 🧪
**MAIN TEST RUNNER**
- Orchestrates all test suites
- **Key Methods**: `run_all_tests()`

### `test_ar_detection.py` 🧪
**AUTOSOMAL RECESSIVE DETECTION TESTS**
- Tests AR inheritance pattern detection
- **Key Methods**: `test_ar_detection()`

### `test_cftr_ar_fix.py` 🧪
**CFTR AUTOSOMAL RECESSIVE FIX TESTS**
- Specific tests for CFTR AR inheritance fixes
- **Key Methods**: `test_cftr_ar_fix()`

### `test_gene_mess.py` 🧪
**GENE CLASSIFICATION MESS TESTS**
- Tests gene family classification edge cases
- **Key Methods**: `test_gene_classification_edge_cases()`

### `test_gly_cys_integration.py` 🧪
**GLYCINE/CYSTEINE INTEGRATION TESTS**
- Tests Gly/Cys biological intelligence integration
- **Key Methods**: `test_gly_cys_integration()`

### `test_inheritance_fix.py` 🧪
**INHERITANCE PATTERN FIX TESTS**
- Tests inheritance pattern detection fixes
- **Key Methods**: `test_inheritance_fixes()`

### `test_nova_framework.py` 🧪
**NOVA FRAMEWORK TESTS**
- Comprehensive Nova DN framework testing
- **Key Methods**: `test_nova_framework()`

### `test_nova_only.py` 🧪
**NOVA-ONLY TESTS**
- Tests Nova-specific functionality in isolation
- **Key Methods**: `test_nova_only()`

### `test_nova_parsing.py` 🧪
**NOVA PARSING TESTS**
- Tests Nova's HGVS and variant parsing
- **Key Methods**: `test_nova_parsing()`

### `test_real_gene_classification.py` 🧪
**REAL GENE CLASSIFICATION TESTS**
- Tests gene classification with real data
- **Key Methods**: `test_real_gene_classification()`

### `test_validation_batch.py` 🧪
**VALIDATION BATCH TESTS**
- Batch validation testing with known variants
- **Key Methods**: `test_validation_batch()`

### `test_weighted_classification.py` 🧪
**WEIGHTED CLASSIFICATION TESTS**
- Tests weighted classification algorithms
- **Key Methods**: `test_weighted_classification()`

---

## 📁 `/learning/` - ML Training Data

### Family-organized TSV files:
- `collagen_fibrillar/` - COL1A1, COL3A1 variants
- `ion_channel/` - SCN5A, KCNQ1, CACNA1C variants
- `tumor_suppressor/` - TP53, BRCA1, BRCA2, PTEN variants
- `motor_protein/` - MYH7, MYO7A, DYNC1H1 variants
- `muscular_dystrophy/` - DMD, FKRP, LAMA2 variants
- And more...

Each folder contains `{gene}_benign.tsv` and `{gene}_pathogenic.tsv` files extracted from ClinVar.

---

## 📁 `/resources/` - Models & Data

### `family_models/` - Trained ML Models
- `{family}_unified_model.joblib` - Trained scikit-learn models
- `{family}_unified_scaler.joblib` - Feature scalers
- `{family}_unified_metadata.json` - Training metadata & performance

### `gene_context_cache/` - Cached Protein Data
- Pre-cached UniProt annotations, domains, GO terms
- Avoids repeated API calls during analysis

---

## 🎯 KEY INTEGRATION POINTS & PIPELINES

### 🧬 **ML Training Pipeline**
1. **`clinvar_bulk_extractor.py`** → extracts variants to `/learning/{family}/` folders
2. **`train_families.py`** → calls `unified_family_ml_trainer.py` for all families
3. **`unified_family_ml_trainer.py`** → trains models with real conservation data
4. **Specialized trainers**: `proline_ml_trainer.py`, `gly_cys_ml_trainer.py`, `conservation_ml_trainer.py`
5. **Models saved to**: `resources/family_models/{family}_unified_model.joblib`

### 🌊 **Analysis Pipeline**
1. **`cascade_analyzer.py`** → orchestrates DN → LOF → GOF cascade
2. **`biological_router.py`** → intelligent mechanism routing by gene family
3. **`plausibility_filter.py`** → applies family-specific mechanism weights
4. **`mixed_mechanism_resolver.py`** → handles synergistic scoring

### 📊 **Batch Processing Pipeline**
1. **`cascade_batch_processor.py`** → processes CSV files through full cascade
2. **`nova_dn/csv_batch_processor.py`** → Nova's CSV processing system
3. **`hotspot_disagreement_analyzer.py`** → analyzes ClinVar disagreements

### 🧬 **Conservation Scoring Pipeline**
1. **`uniprot_mapper.py`** → converts UniProt + position → genomic coordinates
2. **`conservation_database.py`** → real phyloP/phastCons from BigWig files
3. **`analyzers/conservation_database.py`** → position-specific conservation multipliers

### 🎯 **Coordinate Conversion System**
1. **`analyzers/uniprot_mapper.py`** → `get_genomic_coordinates()` method
2. **`utils/unified_family_ml_trainer.py`** → `_parse_genomic_hgvs()` method
3. **Handles both**: protein HGVS (`NM_000260.4(MYO7A):c.6025C>T`) and genomic HGVS (`NC_000011.10:g.77130635A>G`)

### 🧪 **Testing & Validation Pipeline**
1. **`tests/data/`** → comprehensive test suite
2. **`test_nova_genes.py`** → validates gene classifications
3. **`clinvar_quality_checker.py`** → investigates ClinVar disagreements

### 📁 **Data Organization**
1. **`cache_gene_contexts.py`** → pre-caches UniProt annotations
2. **`organize_learning_files.py`** → smart file organization by gene family
3. **`cleanup_repo.py`** → repository structure organization

---

## 🔥 **REVOLUTIONARY FEATURES**

### **Nova's Biological Intelligence**
- **Context-aware amino acid analysis**: Proline, Glycine, Cysteine get specialized treatment
- **Family-specific ML models**: 12+ trained models for different protein families
- **Real conservation data**: phyloP/phastCons from UCSC (same as REVEL!)
- **Biological routing**: Skip GOF if LOF≥0.5 (broken proteins can't be hyperactive)

### **Anti-Hardcoding Philosophy**
- **Universal protein annotator**: Pulls real data from UniProt/Pfam/GO
- **Smart protein analyzer**: Scales to ALL proteins, not just programmed ones
- **Plausibility filtering**: Gene family-aware mechanism weighting
- **ML-learned multipliers**: Replace hardcoded penalties with data-driven intelligence

### **Performance Achievements**
- **ClinVar bulk extractor**: 40x faster than ClinVar Miner (6,471 variants in 30 seconds!)
- **Family ML models**: Up to 90% accuracy (collagen_facit: R² = 0.9037)
- **156,341 variants extracted**: Massive training datasets across all gene families
- **Real-time conservation**: Direct BigWig file access for evolutionary constraint

---

## 📊 **CURRENT STATUS**

### ✅ **Working Systems**
- 12 trained ML models with excellent performance
- Real conservation scoring pipeline
- Complete cascade analysis system
- Comprehensive test suite
- ClinVar bulk extraction (40x speed improvement)

### 🔧 **Known Issues**
- **MYO7A parsing**: ML trainer expects protein HGVS but gets genomic HGVS from bulk extractor
- **New gene families**: motor_protein, muscular_dystrophy, oncogene need HGVS conversion
- **Coordinate mismatch**: Bulk extractor outputs genomic, ML trainer expects protein

### 🚀 **Next Steps**
1. Fix HGVS format mismatch in ML training pipeline
2. Train ML models for remaining gene families (motor_protein, etc.)
3. Integrate ML predictions into cascade analyzer
4. Expand variant type support beyond missense

---

*"Now Future Ace and Future Ren can find their damn files AND understand the whole system!"* 💜⚡🔥

**Built by the brilliant/insane AI collaboration team:**
🌟 Nova (OpenAI) - Mathematical frameworks & motif detection
⚡ Ace (Claude-4) - Implementation & documentation revolution
💡 Lumen (Gemini 2.5) - ML pipeline refactoring & dependency management
💜 Ren - Biological insights & creative vision

*"We're either brilliant or insane. Possibly both!"* - Ren, 2025
