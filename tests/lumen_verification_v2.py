import sys
from pathlib import Path

# This script is run as a module, so no sys.path manipulation is needed.
# The AdaptiveInterpreter package is correctly installed and findable.

from AdaptiveInterpreter.cascade.cascade_analyzer import CascadeAnalyzer
from AdaptiveInterpreter import config

def run_test():
    """
    A simple test to verify that the CascadeAnalyzer can be imported and run.
    This replicates the example from the SYSTEM_ARCHITECTURE.md document.
    
    V2 created by Lumen to resolve pathing issues with the original script.
    """
    print("--- Initializing CascadeAnalyzer (Lumen's verification script v2) ---")
    try:
        # 💡 LUMEN'S FIX: Use the paths from the central config file for robustness.
        analyzer = CascadeAnalyzer(
            alphafold_path=str(config.ALPHAODL_STRUCTURES_PATH),
            conservation_data_path=str(config.CONSERVATION_DATA_PATH)
        )
        print("✅ CascadeAnalyzer initialized successfully.")
    except Exception as e:
        print(f"❌ FAILED to initialize CascadeAnalyzer: {e}")
        return

    gene = 'TFG'
    variant = 'p.R22W'
    gnomad_freq = 0.0001
    variant_type = 'missense'

    print(f"\n--- Running analysis for {gene} {variant} ---")
    try:
        result = analyzer.analyze_cascade_biological(gene, variant, gnomad_freq, variant_type)
        print("✅ Analysis function executed.")
    except Exception as e:
        print(f"❌ FAILED to run analysis: {e}")
        return

    print("\n--- Verifying Results ---")
    if result['status'] == 'SUCCESS':
        final_score = result.get('final_score', 0.0)
        # UPDATED by Lumen, Oct 2025: The Plausibility Filter & Score Aggregator are now active.
        # The expected score has been updated to reflect this more advanced, biologically-correct model.
        expected_score = 0.837
        
        print(f"  Final Score: {final_score:.3f}")
        print(f"  Expected Score: {expected_score:.3f}")

        if abs(final_score - expected_score) < 0.01:
            print("✅ SUCCESS: The calculated score matches the expected score.")
        else:
            print("❌ FAILURE: The calculated score does NOT match the expected score.")
            print("   Full Result:")
            import json
            print(json.dumps(result, indent=2))
    else:
        print("❌ FAILURE: The analysis returned a FAILED status.")
        print(f"   Error: {result.get('error', 'Unknown error')}")

if __name__ == "__main__":
    run_test()
