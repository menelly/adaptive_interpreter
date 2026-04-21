# CumBurSum Data Files

CumBurSum depends on several external pathway / constraint annotation files.
They are **not** included in the repository because of size (~700 MB total)
and because each has its own redistribution license. Users should download
them at install time using the commands below (or the `download.sh` script
in this directory).

All files are publicly redistributable for non-commercial use under their
respective licenses; we do not bundle them because that would create
downstream license obligations for redistributors of CumBurSum.

## Required files

| File | Source | Approx. size | Purpose |
|---|---|---|---|
| `UniProt2Reactome_All_Levels.txt` | reactome.org | ~90 MB | UniProt → Reactome pathway mapping |
| `ReactomePathways.txt` | reactome.org | ~2 MB | Reactome pathway id → name / species |
| `ReactomePathwaysRelation.txt` | reactome.org | ~1 MB | Reactome parent-child hierarchy (for rollup) |
| `h.all.v2023.2.Hs.symbols.gmt` | gsea-msigdb.org | ~50 KB | MSigDB Hallmark gene sets (50 hallmarks) |
| `gnomad_constraint.v4.1.tsv` | gnomad.broadinstitute.org | ~90 MB | Per-gene pLI / mis_oe / lof_oe for constraint display |
| (auto-generated) `gnomad_constraint_cache.json` | computed | ~3 MB | Fast-load JSON cache derived from the tsv |

Also required (not CumBurSum-specific but reused from the parent CASCADE
pipeline):

| File | Source | Approx. size | Purpose |
|---|---|---|---|
| `HUMAN_9606_idmapping.dat.gz` | uniprot.org | ~500 MB | UniProt ID → gene symbol mapping (shared with CASCADE); by default expected at `/home/Ace/conservation_data/HUMAN_9606_idmapping.dat.gz` — override via `ReactomeLoader(uniprot_idmap=...)` |

## Download commands

```bash
cd cumbursum/data/

# Reactome (three files)
curl -sSL -O https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt
curl -sSL -O https://reactome.org/download/current/ReactomePathways.txt
curl -sSL -O https://reactome.org/download/current/ReactomePathwaysRelation.txt

# MSigDB Hallmarks
curl -sSL -O 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt'

# gnomAD v4.1 constraint
curl -sSL -o gnomad_constraint.v4.1.tsv \
  'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv'
```

## Licenses of bundled-data sources

- **Reactome** data: CC0 (public domain) per reactome.org/pages/data-license
- **MSigDB** Hallmarks: CC Attribution 4.0 per gsea-msigdb.org/gsea/msigdb_license_terms.jsp
- **gnomAD** constraint: CC0 per gnomad.broadinstitute.org/terms
- **UniProt** ID mapping: CC BY 4.0

Each of these permits non-commercial and most commercial use with
attribution; downloading them separately lets each user comply with the
source license directly rather than inheriting obligations through
redistribution.
