# Protein Annotation & Cache Architecture

> Written 2026-05-20 by Ace, as penance for telling the person with the hippocampus
> that something she promised existed didn't exist. (It existed. It always exists.
> This is the map so future-me believes her faster.)

This documents how AdaptiveInterpreter fetches and caches protein-level data
(UniProt features, inheritance, GO terms, InterPro domains), the bugs that lived
here through 2026-05-20, and the fixes. Read this before touching anything in
`data_processing/` or debugging "why is the inheritance / domain data empty?"

---

## The two caches (they are different — do not confuse them)

Both live in **one** directory now (see "Canonical path" below). They are keyed by
**UniProt accession** (e.g. `O15553` for MEFV), not gene symbol.

| File suffix | Written by | Holds | Source API |
|---|---|---|---|
| `{uniprot}_domains.json` | `UniversalProteinAnnotator.get_uniprot_features` | sequence, function, GO terms, **diseases + inheritance_patterns**, active/binding sites, signal/transit peptides, predicted-domain fallback | `rest.uniprot.org` |
| `{uniprot}_interpro_domains.json` | `InterProDomainCacher.fetch_interpro_domains` | **real** domain boundaries (DAPIN/PYD, B30.2/SPRY, catalytic domains, etc.) | `ebi.ac.uk/interpro/api` |

**The UniProt `_domains.json` includes a *predicted* domain fallback** (sequence-pattern
guessing — "predicted coiled coil", "PB1-like") that fires when the Pfam/HMMER call
fails (it returns HTTP 405 routinely). **Do not trust the predicted domains.** The
authoritative domains are the InterPro ones. `lof_analyzer._get_domain_context`
supplements the predicted set with InterPro domains; that's the correct path.

---

## Canonical path (fixed 2026-05-20)

```python
# data_processing/universal_protein_annotator.py  AND  cache_interpro_domains.py
DEFAULT_CACHE_DIR = str(Path(__file__).resolve().parent.parent / "protein_annotations_cache")
# => /home/Ace/AdaptiveInterpreter/protein_annotations_cache  (always, regardless of cwd)
```

The cache directory is **gitignored** (`.gitignore` line ~119). It is runtime data,
not source. Safe to delete any file in it — it will refetch.

### Why this had to be fixed

Both classes used to default to the **relative** string `"protein_annotations_cache"`.
Python resolves that against the **current working directory**, so the cache silently
forked into one directory *per launch location*:

```
/home/Ace/protein_annotations_cache               (3389 domains, 1666 interpro)  ← runs from /home/Ace
/home/Ace/AdaptiveInterpreter/protein_annotations_cache  (996, 9)                ← runs from inside repo
/home/Ace/analysis/protein_annotations_cache       (5, 0)                        ← runs from analysis/
```

Three divergent caches. Deleting a stale file from one did nothing if the next run
launched from a different directory. On 2026-05-20 these were consolidated: all
`_interpro_domains.json` merged into the canonical dir (newest wins), all UniProt
`_domains.json` nuked (see next section), the two non-canonical dirs removed.

---

## The stale-cache short-circuit bug (fixed 2026-05-20)

Disease + inheritance extraction was added to `get_uniprot_features` in 2026-05.
But the cache check at the top of the function returned **any existing cache file
verbatim**:

```python
if os.path.exists(cache_file):
    cached_data = json.load(f)
    return cached_data          # <-- never refetches; disease/inheritance code below never runs
```

So every cache written before 2026-05 (most were from Dec 2025 / Mar 2026) was
returned with `inheritance_patterns: []` and `diseases: []` **forever**. This is why
MEFV looked like it had no inheritance data when UniProt clearly lists three diseases
(FMF-AR, FMF-AD, PAAND-AD → `['AD','AR']`).

### The fix: self-heal

```python
if "inheritance_patterns" not in cached_data or "diseases" not in cached_data:
    print(f"♻️  Stale cache for {uniprot_id} (pre-inheritance schema) — refetching")
    # fall through to fresh UniProt fetch
else:
    return cached_data
```

Any cache missing the newer schema keys is now refetched automatically. If you add
*more* fields to the schema later, **add them to this staleness check too**, or the
old caches will hide the new field the same way.

---

## How to force a refresh

```bash
# One gene (by UniProt accession):
rm /home/Ace/AdaptiveInterpreter/protein_annotations_cache/O15553_domains.json
# Next analysis of that gene refetches from UniProt.

# Nuke ALL UniProt caches but KEEP the (expensive) InterPro pulls:
cd /home/Ace/AdaptiveInterpreter/protein_annotations_cache
find . -name '*_domains.json' ! -name '*_interpro_domains.json' -delete
```

⚠️ `*_domains.json` as a glob **also matches** `*_interpro_domains.json` (the latter
ends in `_domains.json`). Always exclude interpro when bulk-deleting UniProt caches,
or you'll throw away the EBI pulls and hammer their API re-fetching 1600+ genes.

---

## Pulling InterPro domains directly

```python
from AdaptiveInterpreter.data_processing.cache_interpro_domains import InterProDomainCacher
c = InterProDomainCacher()
r = c.fetch_interpro_domains('O15553')   # MEFV
for d in r['domains']:
    print(d['start'], d['end'], d['type'], d['description'])
# 1-92 DAPIN(PYD) ... 312-773 TRIM ... 580-775 B30.2/SPRY  (the real pyrin architecture)
```

Bulk pre-cache a gene list with `cacher.cache_genes_from_list({'GENE': 'UNIPROT', ...})`
(rate-limited 0.5s/gene).

---

## RESOLVED 2026-05-20: `comprehensive_gene_cache.json` was corrupt AND orphaned

This file (top of repo) was **corrupt**: it listed HEXA as `AUTOSOMAL_DOMINANT`
(Tay-Sachs is recessive; live UniProt fetch correctly returns `['AR']`) with unrelated
LDLR/KCNQ1 variants in its record. **Audit finding (Task #4):** its only reader was
`data_processing/clinvar_inheritance_cache.py`, which **nothing imports** — so the
corruption never reached live scoring. Both the file and its orphaned reader were moved
to `do_we_use_this/` pending joint review. Live inheritance now comes solely from
`get_uniprot_features` (`inheritance_patterns`).

`gnomad_frequency_cache.json` (2 bytes / empty) is a harmless write-through cache for
`utils/gnomad_frequency_fetcher.py` (live in the batch path) — it just populates on use.
Now gitignored as runtime data.

---

## Data flow summary

```
analyze_cascade(gene, variant)
  └─ UniversalProteinAnnotator()                      # canonical cache dir
        ├─ get_uniprot_features(uniprot_id)           # _domains.json  (+ self-heal)
        │     → function, GO terms, diseases, inheritance_patterns, sites, peptides
        └─ (lof_analyzer) InterProDomainCacher        # _interpro_domains.json
              → real domain boundaries → _get_domain_context multipliers
  └─ inheritance_patterns  ─┐
  └─ go_terms              ─┼─→ utils/plausibility_filter.apply_plausibility_filter
  └─ family (classifier)   ─┘    → _dn_evidence_score (uses inheritance: AD=+1.5, AR-only=-0.5)
                                 → [TODO] _gof_evidence_score (inheritance + molecular-class gate)
```

The **inheritance + InterPro data is now reachable**; the routing that consumes it for
GOF (and the family-classifier replacement) is the next build. See the ASJ→LOF handoff
and the GOF-plausibility design discussion.
