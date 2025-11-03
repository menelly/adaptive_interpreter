#!/usr/bin/env python3
"""
🔬 BATCH CASCADE ANALYZER FOR REMAINING 22 GENES
Built by Ace with MAXIMUM SASS! 💜

This script runs the cascade analyzer on the 22 genes that just got converted!
"""

import subprocess
import sys
from pathlib import Path

# The 22 genes that just got converted
REMAINING_GENES = [
    "ACTN2", "APOB", "BRCA1", "CDC73", "CTNNA1", "DSC2", "FH", "HNF1B",
    "KCNQ1", "MAX", "MLH1", "MSH6", "MYL2", "NKX2-5", "PHOX2B", "RAD50",
    "RET", "SCN5A", "SDHC", "SMARCA4", "TGFBR2", "TSC2"
]

DATA_DIR = Path("/home/Ace/analysis/acmg_sf_train_split")
OUTPUT_DIR = Path("/home/Ace/analysis/acmg_73_FIXED_run_20251102")
CASCADE_SCRIPT = Path("/home/Ace/AdaptiveInterpreter/cascade/batch_cascade_analyzer.py")

def main():
    print("🔬 BATCH CASCADE ANALYZER FOR REMAINING 22 GENES")
    print("=" * 60)
    print(f"💜 Built by Ace with MAXIMUM SASS!")
    print(f"📊 Processing {len(REMAINING_GENES)} genes")
    print(f"📁 Output: {OUTPUT_DIR}")
    print("=" * 60)
    print()
    
    total = len(REMAINING_GENES)
    success = 0
    failed = 0
    
    for i, gene in enumerate(REMAINING_GENES, 1):
        print(f"[{i}/{total}] 🧬 Processing {gene}...")
        
        input_file = DATA_DIR / gene / f"{gene}.protein_ready.discovery.tsv"
        output_file = OUTPUT_DIR / f"{gene}.cascade.with_interface.tsv"
        
        if not input_file.exists():
            print(f"   ⚠️  Input file not found: {input_file}")
            failed += 1
            continue
        
        if output_file.exists():
            print(f"   ⚠️  Output file already exists, skipping: {output_file}")
            success += 1
            continue
        
        # Run the cascade analyzer
        cmd = [
            "python3",
            str(CASCADE_SCRIPT),
            str(input_file),
            str(output_file)
        ]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600  # 10 minute timeout per gene
            )
            
            if result.returncode == 0:
                print(f"   ✅ Success: {output_file}")
                success += 1
            else:
                print(f"   ❌ Failed with return code {result.returncode}")
                print(f"   Error: {result.stderr[:200]}")
                failed += 1
        except subprocess.TimeoutExpired:
            print(f"   ⏱️  Timeout after 10 minutes")
            failed += 1
        except Exception as e:
            print(f"   ❌ Error: {e}")
            failed += 1
        
        print()
    
    print("=" * 60)
    print("🎉 BATCH CASCADE COMPLETE!")
    print(f"   ✅ Success: {success}")
    print(f"   ❌ Failed:  {failed}")
    print(f"   📊 Total:   {total}")
    print("=" * 60)
    
    if failed > 0:
        sys.exit(1)

if __name__ == "__main__":
    main()

