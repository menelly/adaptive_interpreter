#!/bin/bash
# Batch convert the 22 missing genes from .other.tsv to .protein_ready.discovery.tsv

set -e  # Exit on error

SCRIPT_DIR="/home/Ace/AdaptiveInterpreter/utils"
DATA_DIR="/home/Ace/analysis/acmg_sf_train_split"

# List of genes that need conversion (missing protein_ready.discovery.tsv)
MISSING_GENES=(
    "ACTN2"
    "APOB"
    "BRCA1"
    "CDC73"
    "CTNNA1"
    "DSC2"
    "FH"
    "HNF1B"
    "KCNQ1"
    "MAX"
    "MLH1"
    "MSH6"
    "MYL2"
    "NKX2-5"
    "PHOX2B"
    "RAD50"
    "RET"
    "SCN5A"
    "SDHC"
    "SMARCA4"
    "TGFBR2"
    "TSC2"
)

echo "🔬 Starting batch conversion of 22 missing genes..."
echo "=================================================="
echo ""

TOTAL=${#MISSING_GENES[@]}
CURRENT=0
SUCCESS=0
FAILED=0

for gene in "${MISSING_GENES[@]}"; do
    CURRENT=$((CURRENT + 1))
    echo "[$CURRENT/$TOTAL] Processing $gene..."
    
    INPUT_FILE="$DATA_DIR/$gene/${gene}.other.tsv"
    OUTPUT_FILE="$DATA_DIR/$gene/${gene}.protein_ready.discovery.tsv"
    
    if [ ! -f "$INPUT_FILE" ]; then
        echo "   ⚠️  Input file not found: $INPUT_FILE"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    if [ -f "$OUTPUT_FILE" ]; then
        echo "   ⚠️  Output file already exists, skipping: $OUTPUT_FILE"
        SUCCESS=$((SUCCESS + 1))
        continue
    fi
    
    # Run genomic_to_protein.py as a module
    if (cd /home/Ace && python3 -m AdaptiveInterpreter.utils.genomic_to_protein "$INPUT_FILE" "$OUTPUT_FILE") > /dev/null 2>&1; then
        echo "   ✅ Success: $OUTPUT_FILE"
        SUCCESS=$((SUCCESS + 1))
    else
        echo "   ❌ Failed to convert $gene"
        FAILED=$((FAILED + 1))
    fi
    
    echo ""
done

echo "=================================================="
echo "🎉 Batch conversion complete!"
echo "   Success: $SUCCESS"
echo "   Failed:  $FAILED"
echo "   Total:   $TOTAL"
echo "=================================================="

