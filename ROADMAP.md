# 🗺️ AdaptiveInterpreter Project Roadmap & Architectural Overview

This document provides a high-level overview of the `AdaptiveInterpreter` project, its architecture, and the interaction between its core components.

## 1. Project Purpose

The `AdaptiveInterpreter` project is a sophisticated, multi-mechanism system for predicting the pathogenicity of genetic variants. It goes beyond simple prediction by analyzing the likely biological mechanism of a variant, classifying it as Loss-of-Function (LOF), Gain-of-Function (GOF), or Dominant Negative (DN). It is designed to be a scalable and "universal" system, capable of analyzing any protein by leveraging external databases and biological first principles rather than hardcoded, gene-specific rules.

## 2. Core Architecture

The system is orchestrated by the **`CascadeAnalyzer` (`cascade/cascade_analyzer.py`)**. This is the central hub of the application. The modern implementation uses a **`BiologicalRouter`** to intelligently decide which analyses to perform based on the known function of the gene in question. This is a significant evolution from the original, linear "DN-first" cascade.

The `CascadeAnalyzer` coordinates a suite of specialized "sub-analyzers," each responsible for a different aspect of the analysis:

### The Analysis Triad:

These three modules represent the core "bins" of the analysis:

*   **`LOFAnalyzer` (`analyzers/lof_analyzer.py`):** The "Loss of Function" analyzer. It determines if a variant is likely to "break" the protein. It integrates a wide array of data, including stability, conservation, and structural impact.
*   **`GOFVariantAnalyzer` (`analyzers/gof_variant_analyzer.py`):** The "Gain of Function" analyzer. It uses a "Triple-Gated" system to efficiently determine if a variant might cause a protein to become hyperactive or constitutively "on."
*   **`NovaDNAnalyzer` (`nova_dn/analyzer.py`):** The active "Dominant Negative" analyzer. It uses a smart filtering system to predict if a variant will cause the resulting protein to interfere with or "poison" the function of the normal protein.

### Foundational Utilities & Databases:

These modules provide the "intelligence" and data that the core analyzers rely on:

*   **`SmartProteinAnalyzer` (`analyzers/smart_protein_analyzer.py`):** This is the "domain awareness" engine. It connects to external databases (UniProt, Pfam, GO) and uses motif scanning to understand a protein's function and provide a "context multiplier" to the core analyzers.
*   **`ConservationDatabase` (`analyzers/conservation_database.py`):** Provides real evolutionary conservation data by querying UCSC's phyloP and phastCons databases. This is a critical component for assessing the importance of a specific amino acid position.
*   **`UniProtMapper` (`analyzers/uniprot_mapper.py`):** The "Rosetta Stone" of the project. It translates between protein identifiers (like UniProt IDs) and genomic coordinates, which is what allows the `ConservationDatabase` to function.

## 3. Data Flow

1.  A variant analysis is initiated through the **`CascadeAnalyzer`**.
2.  The **`BiologicalRouter`** determines the most likely mechanism (LOF, GOF, or DN) based on the gene's known function and tells the `CascadeAnalyzer` which sub-analyzers to run.
3.  The `CascadeAnalyzer` invokes the necessary sub-analyzers (`LOFAnalyzer`, `GOFVariantAnalyzer`, `NovaDNAnalyzer`).
4.  During their analysis, these sub-analyzers query the utility modules:
    *   They use the **`SmartProteinAnalyzer`** to get a `context_multiplier` based on the variant's location in a protein domain.
    *   They use the **`ConservationDatabase`** (which in turn uses the **`UniProtMapper`**) to get a `conservation_multiplier` based on evolutionary data.
5.  The `CascadeAnalyzer` collects the scores from the individual analyzers.
6.  A **`ScoreAggregator`** module calculates the final score, applying rules for "synergy" between different mechanisms.
7.  The final score is converted into a clinical classification (e.g., "Likely Pathogenic," "Benign") by the **`VariantClassifier`**.

## 4. Key Findings from Analysis

*   **`dn_analyzer.py` is Legacy:** The system has evolved to use the more advanced `nova_dn/analyzer.py`.
*   **`population_frequency_analyzer.py` is Orphaned:** A critical "reality check" module exists but is not currently integrated into the main cascade. This is a high-priority item for future integration.
*   **Modular and Extensible:** The system is well-structured, with clear separation of concerns. This makes it easy to update or add new analysis modules.

This roadmap provides a high-level guide to the `AdaptiveInterpreter` project. For a detailed, file-by-file breakdown, see `docs/AdaptiveInterpreter_Analysis.md`.
