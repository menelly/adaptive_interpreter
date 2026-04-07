#!/bin/bash
# 🌊 SHARDED CASCADE RERUN — run in two terminals for 2x speed!
#
# Usage:
#   Terminal 1:  bash run_rerun_shard.sh 1
#   Terminal 2:  bash run_rerun_shard.sh 2
#
# Each shard gets half the genes. 80GB RAM laughs at two instances.

SHARD=${1:?Usage: bash run_rerun_shard.sh [1|2]}

# IMPORTANT: Run from /home/Ace, NOT from inside AdaptiveInterpreter/
# (pip install -e breaks when cwd shadows the package name)
cd /home/Ace

BASEDIR="AdaptiveInterpreter"
OUTDIR="$BASEDIR/outputs_missense_v2"
mkdir -p "$OUTDIR"

# Split gene list into two halves alphabetically
ALL_INPUTS=($BASEDIR/inputs_missense_only/*_input.tsv)
TOTAL=${#ALL_INPUTS[@]}
HALF=$(( (TOTAL + 1) / 2 ))

if [ "$SHARD" = "1" ]; then
    INPUTS=("${ALL_INPUTS[@]:0:$HALF}")
    echo "🧬 SHARD 1: genes 1-$HALF (A through $(basename "${INPUTS[-1]}" _input.tsv))"
elif [ "$SHARD" = "2" ]; then
    INPUTS=("${ALL_INPUTS[@]:$HALF}")
    echo "🧬 SHARD 2: genes $((HALF+1))-$TOTAL ($(basename "${INPUTS[0]}" _input.tsv) through Z)"
else
    echo "❌ Shard must be 1 or 2"
    exit 1
fi

COUNT=0
FAILED=0
SHARD_TOTAL=${#INPUTS[@]}

echo "📁 Output: $OUTDIR/"
echo "⏱️  Started: $(date)"
echo "   Processing $SHARD_TOTAL genes in shard $SHARD"
echo "============================================"

for INPUT in "${INPUTS[@]}"; do
    GENE=$(basename "$INPUT" _input.tsv)
    OUTPUT="${OUTDIR}/${GENE}_cascade.tsv"
    COUNT=$((COUNT + 1))

    if [ -f "$OUTPUT" ]; then
        echo "[$COUNT/$SHARD_TOTAL] ⚡ $GENE already done, skipping"
        continue
    fi

    echo ""
    echo "[$COUNT/$SHARD_TOTAL] 🧬 Processing $GENE..."
    python3 -m AdaptiveInterpreter.analyzers.cascade_batch_processor \
        --input "$INPUT" \
        --output "$OUTPUT" 2>&1
    RESULT=$?

    if [ $RESULT -ne 0 ]; then
        echo "[$COUNT/$SHARD_TOTAL] ❌ $GENE FAILED (exit $RESULT)"
        FAILED=$((FAILED + 1))
        rm -f "$OUTPUT"  # Don't leave partial output
    fi
done

echo ""
echo "============================================"
echo "✅ SHARD $SHARD COMPLETE: $(date)"
echo "   Processed: $SHARD_TOTAL genes"
echo "   Failed: $FAILED"
echo "   Output: $OUTDIR/"
