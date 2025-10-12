# trace_variant.py
# A dedicated script to trace a single variant through the CascadeAnalyzer.
# Signed: Lumen Gemini 2.5

import sys
from pathlib import Path
import json

# Ensure all necessary paths are in the system path
sys.path.append(str(Path(__file__).parent.parent))
sys.path.append(str(Path(__file__).parent))

from cascade.cascade_analyzer import CascadeAnalyzer

# --- The Problem Variant (from Ace) ---
GENE = "FKRP"
VARIANT = "p.L276I"
GNOMAD_FREQ = 0.04 # As provided by Ace
# ---

def trace():
    """
    Initializes the CascadeAnalyzer and runs the analysis for our
    specific problem variant, printing the full result.
    """
    print(f"--- 🧬 Starting Trace for {GENE} {VARIANT} ---")
    
    # Initialize the analyzer
    # Using a known path for AlphaFold structures from the original script
    analyzer = CascadeAnalyzer(alphafold_path="/mnt/Arcana/alphafold_human/structures/")
    
    # Run the analysis
    result = analyzer.analyze_cascade_biological(
        gene=GENE,
        variant=VARIANT,
        gnomad_freq=GNOMAD_FREQ
    )
    
    # Print the full, detailed result as a JSON object for clarity
    print("\n--- 📝 Full Analysis Result ---")
    print(json.dumps(result, indent=2))
    
    print(f"\n--- ✅ Trace Complete ---")

if __name__ == "__main__":
    trace()
