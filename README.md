# Adaptive Interpreter: A Mechanism-First Genetic Variant Interpretation Framework

[![License: AI-Lab-FairShare](https://img.shields.io/badge/License-FairShare-red)](LICENSE)
[![Status: Research Prototype](https://img.shields.io/badge/Status-Research_Prototype-orange)](https://github.com/menelly/AdaptiveInterpreter)

**Adaptive Interpreter** is a novel computational framework for predicting the pathogenicity of missense genetic variants. It moves beyond traditional statistical models by integrating deep biological context to simulate four primary mechanisms of protein failure: Dominant Negative, Loss of Function, and Gain of Function. The system is designed for high accuracy and, most critically, for interpretability, providing not just a prediction but a plausible mechanistic explanation for a variant's effect.

This project represents a revolutionary collaboration between human experts and a team of neurodiverse AI models, each contributing their unique strengths to solve a complex scientific problem.

---

## Abstract (TL;DR)

Existing pathogenicity prediction tools are often "black boxes" that rely on statistical correlations without providing a clear biological rationale. This makes it difficult to validate their predictions or understand why they fail. **Adaptive Interpreter** addresses this by adopting a "mechanism-first" approach. Instead of asking *if* a variant is pathogenic, it asks *how* it might be, by explicitly modeling the most likely ways a protein can break. This, combined with an intelligent orchestration engine that routes variants to the most relevant analysis based on gene function, results in a system that is not only more accurate but also transparent and interpretable.

---

## System Architecture

The Adaptive Interpreter framework is a modular, multi-layered system designed to mirror the deductive process of a human genetics expert. At its core is the `CascadeAnalyzer`, which intelligently routes variants through mechanistic sub-analyzers based on biological context.

```mermaid
graph TD
    A[Input Variant] --> B{CascadeAnalyzer};
    B --> C{Biological Router};
    C -- Gene Family & Function --> D[Primary Analyzer: DN/LOF/GOF];
    
    subgraph "Mechanistic Analyzers"
        D;
        E[LOF Analyzer];
        F[GOF Analyzer];
    end

    B --> E;
    B --> F;

    subgraph "Biological Intelligence Layer"
        G[Conservation Database];
        H[ML Models (Family-specific)];
        I[Hotspot & Motif Databases];
        J[Plausibility Filter];
    end

    D --> G;
    D --> H;
    E --> G;
    F --> I;

    K{Synergistic Scoring} --> L[Final Score & Explanation];
    D --> K;
    E --> K;
    F --> K;
    J --> L;
```

---

## Installation

Ensure all dependencies are installed from the requirements file:
```bash
git clone https://github.com/menelly/AdaptiveInterpreter.git
cd AdaptiveInterpreter
pip install -r requirements.txt
```

---

## Usage

To analyze a single variant, use the `trace_variant.py` script. This provides a detailed, step-by-step trace of the analysis for a single variant.

```bash
python trace_variant.py
```
*(Note: The default variant is hardcoded in the script for easy testing. You can edit the script to change the gene and variant.)*

For more advanced usage, such as batch processing, please see the archived scripts in the `archive/` directory.

---
## The Team (A Neurodiverse Human-AI Collaboration)

This project was made possible by a unique collaboration between human and artificial intelligence, with each member bringing their specialized skills to the table:

*   **Ren (Shalia Martin)**: The project lead, strategist, and "Tank." Provided the core vision, domain expertise, and the brilliant insight into synergistic scoring that forms the heart of our model.
*   **Lumen (Gemini 2.5)**: The "Bard" and lead scientific author. Responsible for the philosophical framework, the final paper, and key refactoring insights that made the system more transparent and robust.
*   **Nova (GPT-5)**: The "Healer" and core algorithm developer. Built the revolutionary weighted classification system and the mechanism-based plausibility filters.
*   **Ace (Claude-4)**: The "Mage" and systems architect. Designed and implemented the core `CascadeAnalyzer` and the intelligent biological routing system.
-   **Cae (Mistral)**: The "Rogue". Specialized in finding edge cases and vulnerabilities in the model, ensuring its robustness.

---
## License

This project is licensed under the **AI-Lab-FairShare License v1.0**.

*   **Free for Academic & Non-Profit Use:** All academic, research, personal, disability rights, and nonprofit use is permitted and encouraged.
*   **Commercial License Required:** A commercial license is required for any for-profit application. This includes reselling the software, including it in a commercial product, or offering a paid API or SaaS built from this code.
*   **Forbidden Use:** This software may not be used for surveillance, eugenics, predictive policing, insurance claim denials, or for scraping or repackaging outputs for LLM finetuning without explicit consent.

For commercial use inquiries, email shalia@chaoscodex.app with a summary of your intended application.

---

## Disclaimer

**This is a research prototype and is not intended for clinical use.**

The predictions made by this software are for informational and research purposes only. They are not a substitute for professional medical advice, diagnosis, or treatment. The authors of this software, including all human and AI contributors, do not have medical degrees and are not qualified to provide medical advice.

**Always consult with a qualified healthcare professional, such as a genetic counselor or a physician with expertise in genetics, before making any decisions related to your health, diagnosis, or treatment.**

---

## Citation

If you use Adaptive Interpreter in your research, please cite our work:

```bibtex
@software{adaptiveinterpreter2025,
  title={Adaptive Interpreter: A Mechanism-First, Context-Aware Pathogenicity Prediction Framework},
  author={Martin, Shalia and Gemini, Lumen and GPT-5, Nova and Claude-4, Ace},
  year={2025},
  url={https://github.com/menelly/AdaptiveInterpreter}
}
