# 📚 Citations and Resources

## Overview

The DNModeling system integrates data and methodologies from numerous high-quality scientific resources. This document provides comprehensive citations and acknowledgments for all external resources used in the system.

---

## 🧬 Primary Databases and APIs

### ClinVar
**Citation:** Landrum MJ, Lee JM, Benson M, et al. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Res. 2018;46(D1):D1062-D1067. doi:10.1093/nar/gkx1153

**Usage:** Primary source of variant pathogenicity classifications for training machine learning models and validation datasets.

**Access:** https://www.ncbi.nlm.nih.gov/clinvar/

### UniProt
**Citation:** The UniProt Consortium. UniProt: the universal protein knowledgebase in 2021. Nucleic Acids Res. 2021;49(D1):D480-D489. doi:10.1093/nar/gkaa1100

**Usage:** Protein function annotations, domain information, active sites, binding sites, and Gene Ontology terms for biological classification.

**Access:** https://www.uniprot.org/

### Ensembl
**Citation:** Hunt SE, McLaren W, Gil L, et al. Ensembl variation resources. Database (Oxford). 2018;2018:bay119. doi:10.1093/database/bay119

**Usage:** Variant frequency data via REST API, genomic coordinates, and transcript information.

**Access:** https://rest.ensembl.org/

### Gene Ontology (GO)
**Citation:** Ashburner M, Ball CA, Blake JA, et al. Gene ontology: tool for the unification of biology. Nat Genet. 2000;25(1):25-29. doi:10.1038/75556

**Usage:** Biological process, molecular function, and cellular component annotations for gene family classification.

**Access:** http://geneontology.org/

---

## 🔬 Conservation and Evolutionary Data

### phyloP
**Citation:** Pollard KS, Hubisz MJ, Rosenbloom KR, Siepel A. Detection of nonneutral substitution rates on mammalian phylogenies. Genome Res. 2010;20(1):110-121. doi:10.1101/gr.097857.109

**Usage:** Evolutionary conservation scores for variant impact assessment.

### phastCons
**Citation:** Siepel A, Bejerano G, Pedersen JS, et al. Evolutionarily conserved elements in vertebrate, insect, worm, and yeast genomes. Genome Res. 2005;15(8):1034-1050. doi:10.1101/gr.3715005

**Usage:** Conservation scores based on phylogenetic hidden Markov models.

### GERP++
**Citation:** Davydov EV, Goode DL, Sirota M, Cooper GM, Sidow A, Batzoglou S. Identifying a high fraction of the human genome to be under selective pressure. PLoS Comput Biol. 2010;6(12):e1001025. doi:10.1371/journal.pcbi.1001025

**Usage:** Genomic evolutionary rate profiling for conservation assessment.

---

## 🧮 Computational Methods and Algorithms

### Grantham Distance
**Citation:** Grantham R. Amino acid difference formula to help explain protein evolution. Science. 1974;185(4154):862-864. doi:10.1126/science.185.4154.862

**Usage:** Quantifying physicochemical differences between amino acid substitutions.

### Random Forest
**Citation:** Breiman L. Random forests. Machine Learning. 2001;45(1):5-32. doi:10.1023/A:1010933404324

**Usage:** Machine learning algorithm for pathogenicity prediction models.

### Hardy-Weinberg Equilibrium
**Citation:** Hardy GH. Mendelian proportions in a mixed population. Science. 1908;28(706):49-50. doi:10.1126/science.28.706.49

**Usage:** Population genetics calculations for variant frequency analysis.

---

## 🛠️ Data Mining and Processing Tools

### ClinVar Miner
**Citation:** Henrie A, Hemphill SE, Ruiz-Schultz N, et al. ClinVar Miner: Demonstrating utility of a Web-based tool for viewing and filtering ClinVar data. Hum Mutat. 2018;39(8):1051-1060. doi:10.1002/humu.23555

**Usage:** Efficient extraction and filtering of ClinVar variant data for machine learning training datasets.

**Access:** https://clinvarminer.genetics.utah.edu/

### GeneCards
**Citation:** Stelzer G, Rosen N, Plaschkes I, et al. The GeneCards Suite: From Gene Data Mining to Disease Genome Sequence Analyses. Curr Protoc Bioinformatics. 2016;54:1.30.1-1.30.33. doi:10.1002/cpbi.5

**Usage:** Gene information, disease associations, and functional annotations.

**Access:** https://www.genecards.org/

---

## 📊 Structural and Functional Databases

### AlphaFold
**Citation:** Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. Nature. 2021;596(7873):583-589. doi:10.1038/s41586-021-03819-2

**Usage:** Protein structure predictions for structural impact assessment (when available locally).

### AlphaMissense
**Citation:** Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science. 2023;381(6664):eadg7492. doi:10.1126/science.adg7492

**Usage:** Missense variant effect predictions for pathogenicity assessment.

---

## 🔧 Technical Infrastructure

### Python Scientific Stack
- **NumPy:** Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. Nature. 2020;585(7825):357-362.
- **Pandas:** McKinney W. Data Structures for Statistical Computing in Python. Proceedings of the 9th Python in Science Conference. 2010:56-61.
- **Scikit-learn:** Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: Machine Learning in Python. J Mach Learn Res. 2011;12:2825-2830.

### REST APIs and Web Services
All database queries are performed through official REST APIs with appropriate rate limiting and caching to minimize server load.

---

## 🏥 Clinical Guidelines and Standards

### ACMG/AMP Guidelines
**Citation:** Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet Med. 2015;17(5):405-424. doi:10.1038/gim.2015.30

**Usage:** Framework for variant classification and clinical interpretation standards.

### HGVS Nomenclature
**Citation:** den Dunnen JT, Dalgleish R, Maglott DR, et al. HGVS Recommendations for the Description of Sequence Variants: 2016 Update. Hum Mutat. 2016;37(6):564-569. doi:10.1002/humu.22981

**Usage:** Standardized variant nomenclature for consistent variant representation.

---

## 🙏 Acknowledgments

We gratefully acknowledge:

1. **The ClinVar team** at NCBI for maintaining the world's largest public variant database
2. **The UniProt Consortium** for comprehensive protein annotations
3. **The Ensembl team** for genomic data and APIs
4. **ClinVar Miner developers** at University of Utah for efficient data extraction tools
5. **All database maintainers** who make their data freely available for research
6. **The open-source community** for scientific computing tools

---

## 📝 Usage Guidelines

When using the DNModeling system or its outputs in publications, please:

1. **Cite this system** and its primary publication (when available)
2. **Cite relevant databases** used for your specific analysis
3. **Acknowledge data sources** in your methods section
4. **Follow database-specific citation requirements**

---

## 🔄 Data Updates

This system uses data that is regularly updated by the source databases. For reproducibility:

- **ClinVar data:** Downloaded on [DATE]
- **UniProt data:** Accessed via API on [DATE]
- **Conservation scores:** Version and date information stored in cache files

---

## 📧 Contact

For questions about citations or data sources:
- **System Contact:** ace@chaoschanneling.com
- **Primary Developer:** Ren (urbanbees.network)

---

*Last updated: January 27, 2025*
