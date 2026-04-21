#!/usr/bin/env bash
# Downloads all pathway annotation + constraint files CumBurSum needs.
# Run from the cumbursum/data/ directory. See README.md for source licenses.

set -euo pipefail

cd "$(dirname "$0")"

echo "Downloading Reactome files..."
curl -sSL -O https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt
curl -sSL -O https://reactome.org/download/current/ReactomePathways.txt
curl -sSL -O https://reactome.org/download/current/ReactomePathwaysRelation.txt

echo "Downloading MSigDB Hallmarks..."
curl -sSL -O 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt'

echo "Downloading gnomAD v4.1 constraint..."
curl -sSL -o gnomad_constraint.v4.1.tsv \
  'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv'

echo
echo "Done. Files in $(pwd):"
ls -la *.txt *.gmt *.tsv 2>/dev/null | awk '{print "  "$NF, "("$5" bytes)"}'
echo
echo "Run the CLI from the repo root, e.g.:"
echo "    python3 cumbursum.py --cascade my_variants.tsv --output report.md"
