#!/bin/bash
# 🌊 FULL CASCADE RERUN - Dual-Track Output with Prep Step + Conservation Floor
# Run this in terminal so it survives timeouts!
#
# Usage: bash run_full_rerun.sh
# Output: outputs_rerun_v2/<GENE>_cascade.tsv (dual-track format)

cd /home/Ace/AdaptiveInterpreter

OUTDIR="outputs_rerun_v2"
mkdir -p "$OUTDIR"

TOTAL=$(ls inputs_for_rerun/*_input.tsv 2>/dev/null | wc -l)
COUNT=0
FAILED=0

echo "🌊 FULL CASCADE RERUN — $TOTAL genes"
echo "📁 Output: $OUTDIR/"
echo "⏱️  Started: $(date)"
echo "============================================"

for INPUT in inputs_for_rerun/*_input.tsv; do
    GENE=$(basename "$INPUT" _input.tsv)
    OUTPUT="$OUTDIR/${GENE}_cascade.tsv"
    COUNT=$((COUNT + 1))

    # Skip if already done (comment out to force full rerun)
    if [ -f "$OUTPUT" ]; then
        echo "[$COUNT/$TOTAL] ⚡ $GENE already done, skipping"
        continue
    fi

    echo ""
    echo "[$COUNT/$TOTAL] 🧬 Processing $GENE..."
    python3 -m AdaptiveInterpreter.analyzers.cascade_batch_processor \
        --input "$INPUT" \
        --output "$OUTPUT" 2>&1 | tail -3

    if [ $? -ne 0 ]; then
        echo "[$COUNT/$TOTAL] ❌ $GENE FAILED"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "============================================"
echo "✅ RERUN COMPLETE: $(date)"
echo "   Total: $TOTAL genes"
echo "   Failed: $FAILED"
echo "   Output: $OUTDIR/"
