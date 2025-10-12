# 📊 DNModeling Data Sources

This document lists the external data files required for the DNModeling system to function correctly.

## 1. Conservation Data

-   **Source:** UCSC Genomics Institute
-   **Files:**
    -   `hg38.phyloP100way.bw`: Evolutionary conservation scores based on the phyloP algorithm.
    -   `hg38.phastCons100way.bw`: Evolutionary conservation scores based on the phastCons algorithm.
-   **Expected Location:** As defined in `config.py`, defaults to `~/genetics_data/conservation/`.
-   **Download:** These files can be downloaded from the [UCSC Genome Browser download page](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/) and [phastCons100way page](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/).

## 2. gnomAD Data

-   **Source:** Genome Aggregation Database (gnomAD)
-   **Files:** `gnomad.genomes.v4.0.sites.chr*.vcf.bgz`
-   **Expected Location:** As defined in `config.py`, defaults to `~/genetics_data/gnomad/`.
-   **Download:** These files can be downloaded from the [gnomAD downloads page](https://gnomad.broadinstitute.org/downloads).

## 3. AlphaFold Protein Structures

-   **Source:** AlphaFold Protein Structure Database
-   **Files:** PDB files for human proteins (e.g., `AF-P04637-F1-model_v4.pdb`).
-   **Expected Location:** As defined in `config.py`, defaults to `~/genetics_data/alphafold_human/structures/`.
-   **Download:** The complete set of human protein structures can be downloaded from the [AlphaFold database](https://alphafold.ebi.ac.uk/download).

## 4. UniProt ID Mapping and Domain Data

-   **Source:** UniProt
-   **Files:**
    -   `uniprot_sprot.dat.gz`: The main UniProt data file.
    -   `uniprot_id_mapping.dat.gz`: ID mapping data.
-   **Expected Location:** The `UniProtMapper` class in `analyzers/uniprot_mapper.py` will download and cache these files automatically to a subdirectory within the `conservation_data` path.
-   **Download:** No manual download is required.
