#!/usr/bin/env python3
"""
Test Nova's weighted classification system on Ren's gene mess!
Let's see how well it handles a diverse set of real genes.
"""

import sys
sys.path.append('.')

from plausibility_filter import classify_gene_family

def test_gene_mess():
    """Test Ren's gene list to see classification accuracy"""
    
    # Ren's gene mess with expected functional families
    gene_list = [
        "ACBD5", "ACMSD", "AFAP1", "AHNAK", "B3GLCT", "BAG6", "BMPR2", "BTD", 
        "C2CD3", "CACNA1I", "CDK15", "CHRM5", "CLASP1", "CNN2", "DLD", "FANCD2", 
        "FKRP", "FTCD", "G6PD", "GADL1", "GFM1", "GXYLT1", "HFE", "KANK1", 
        "KCNMA1", "KRT18", "LIN7A", "MCM4", "MMAA", "MTCH2", "MYH15", "MYL6B", 
        "MYO7A", "NCAPD3", "NKIRAS1", "NPHS2", "PDE6B", "PHTF2", "POLQ", "PRB4", 
        "PYGL", "RBMX", "RPLP0", "SEC22A", "SERPINA1", "SETDB1", "SLC25A5", 
        "SLC9B1", "SPTB", "TFG", "TYW1", "WDFY3", "ZSWIM6", "POLR1C", "SLC6A2", 
        "AIP", "AMPD1", "KIR2DL4", "IRS2", "ATP5F1A", "NSMCE3", "KMT2B", 
        "GET1-SH3BGR", "DYSF"
    ]
    
    # Expected classifications based on known function
    expected_families = {
        "CACNA1I": "ION_CHANNEL",  # Calcium channel
        "KCNMA1": "ION_CHANNEL",   # Potassium channel  
        "CHRM5": "ION_CHANNEL",    # Muscarinic receptor (ligand-gated)
        "SLC25A5": "ION_CHANNEL",  # Mitochondrial transporter
        "SLC9B1": "ION_CHANNEL",   # Sodium/hydrogen exchanger
        "SLC6A2": "ION_CHANNEL",   # Sodium transporter
        
        "BMPR2": "TUMOR_SUPPRESSOR", # TGF-beta family receptor
        "FANCD2": "TUMOR_SUPPRESSOR", # DNA repair (Fanconi anemia)
        "POLQ": "TUMOR_SUPPRESSOR",   # DNA polymerase (repair)
        "MCM4": "TUMOR_SUPPRESSOR",   # DNA replication/repair
        "AIP": "TUMOR_SUPPRESSOR",    # Aryl hydrocarbon receptor interacting protein
        
        "KMT2B": "TRANSCRIPTION_FACTOR", # Histone methyltransferase
        "SETDB1": "TRANSCRIPTION_FACTOR", # Histone methyltransferase
        "POLR1C": "TRANSCRIPTION_FACTOR", # RNA polymerase
        
        "MYH15": "MOTOR_PROTEIN",  # Myosin
        "MYO7A": "MOTOR_PROTEIN",  # Myosin
        
        "RPLP0": "RIBOSOMAL_PROTEIN", # Ribosomal protein
        
        "FKRP": "MUSCULAR_DYSTROPHY", # Fukutin-related protein
        "DYSF": "MUSCULAR_DYSTROPHY", # Dysferlin

        "KRT18": "INTERMEDIATE_FILAMENT", # Keratin
        "SPTB": "CYTOSKELETON_POLYMER",   # Spectrin

        # Test new categories
        "ATP5F1A": "METABOLIC_ENZYME",    # ATP synthase (DN effects proven!)
        "TFG": "SCAFFOLD_ADAPTOR",        # TRK-fused gene (DN effects proven!)
        "G6PD": "METABOLIC_ENZYME",       # Glucose-6-phosphate dehydrogenase
        "SERPINA1": "SIGNALING_REGULATOR", # Alpha-1-antitrypsin (protease inhibitor)

        # Many others will likely be GENERAL or metabolic
    }
    
    print("🧬 TESTING NOVA'S SYSTEM ON REN'S GENE MESS")
    print("=" * 60)
    print(f"📊 Testing {len(gene_list)} genes")
    print()
    
    results = {}
    correct = 0
    total_expected = 0
    
    # Realistic functional descriptions AND GO terms based on gene names/known functions
    gene_data = {
        "CACNA1I": {
            "function": "Voltage-dependent T-type calcium channel subunit alpha-1I",
            "go_terms": ["calcium channel activity", "voltage-gated calcium channel", "ion transport", "calcium ion transport"]
        },
        "KCNMA1": {
            "function": "Calcium-activated potassium channel subunit alpha-1",
            "go_terms": ["potassium channel activity", "calcium-activated potassium channel", "ion transport", "potassium ion transport"]
        },
        "CHRM5": {
            "function": "Muscarinic acetylcholine receptor M5",
            "go_terms": ["G-protein coupled receptor activity", "ligand-gated ion channel", "neurotransmitter receptor", "signal transduction"]
        },
        "SLC25A5": {
            "function": "ADP/ATP translocase 2",
            "go_terms": ["solute carrier", "transporter", "mitochondrial carrier", "adenine nucleotide transport"]
        },
        "SLC9B1": {
            "function": "Sodium/hydrogen exchanger 9B1",
            "go_terms": ["sodium:hydrogen antiporter activity", "ion transport", "sodium ion transport", "pH regulation"]
        },
        "SLC6A2": {
            "function": "Sodium-dependent noradrenaline transporter",
            "go_terms": ["sodium:neurotransmitter symporter activity", "neurotransmitter transport", "sodium ion transport", "ion transport"]
        },

        "BMPR2": {
            "function": "Bone morphogenetic protein receptor type-2",
            "go_terms": ["transforming growth factor beta receptor activity", "serine/threonine kinase activity", "signal transduction", "negative regulation of cell proliferation"]
        },
        "FANCD2": {
            "function": "Fanconi anemia group D2 protein",
            "go_terms": ["DNA repair", "homologous recombination", "DNA damage response", "double-strand break repair"]
        },
        "POLQ": {
            "function": "DNA polymerase theta",
            "go_terms": ["DNA polymerase activity", "DNA repair", "double-strand break repair", "DNA replication"]
        },
        "MCM4": {
            "function": "DNA replication licensing factor MCM4",
            "go_terms": ["DNA helicase activity", "DNA replication", "DNA repair", "cell cycle"]
        },
        "AIP": {
            "function": "AH receptor-interacting protein",
            "go_terms": ["tumor suppressor", "negative regulation of cell growth", "protein binding", "transcriptional regulation"]
        },

        "KMT2B": {
            "function": "Histone-lysine N-methyltransferase 2B",
            "go_terms": ["histone methyltransferase activity", "transcription factor", "chromatin binding", "gene expression regulation"]
        },
        "SETDB1": {
            "function": "Histone-lysine N-methyltransferase SETDB1",
            "go_terms": ["histone methyltransferase activity", "transcriptional repressor", "chromatin binding", "gene silencing"]
        },
        "POLR1C": {
            "function": "DNA-directed RNA polymerase I subunit RPC3",
            "go_terms": ["RNA polymerase activity", "transcription factor", "DNA binding", "RNA synthesis"]
        },

        "MYH15": {
            "function": "Myosin-15",
            "go_terms": ["motor activity", "actin binding", "myosin", "muscle contraction"]
        },
        "MYO7A": {
            "function": "Unconventional myosin-VIIa",
            "go_terms": ["motor activity", "actin binding", "myosin", "cytoskeleton organization"]
        },

        "RPLP0": {
            "function": "60S acidic ribosomal protein P0",
            "go_terms": ["ribosomal protein", "protein synthesis", "ribosome", "translation"]
        },

        "FKRP": {
            "function": "Fukutin-related protein",
            "go_terms": ["glycosyltransferase activity", "muscular dystrophy", "muscle membrane", "protein glycosylation"]
        },
        "DYSF": {
            "function": "Dysferlin",
            "go_terms": ["muscle membrane repair", "muscular dystrophy", "limb-girdle", "membrane fusion"]
        },

        "KRT18": {
            "function": "Keratin type I cytoskeletal 18",
            "go_terms": ["intermediate filament", "cytoskeleton organization", "keratin", "structural constituent"]
        },
        "SPTB": {
            "function": "Spectrin beta chain",
            "go_terms": ["actin binding", "cytoskeleton organization", "membrane cytoskeleton", "spectrin"]
        },

        # Test the new categories
        "ATP5F1A": {
            "function": "ATP synthase F1 subunit alpha",
            "go_terms": ["ATP synthase activity", "mitochondria", "metabolic", "oxidative phosphorylation"]
        },
        "TFG": {
            "function": "TRK-fused gene protein",
            "go_terms": ["vesicle trafficking", "scaffold", "autophagy", "endoplasmic reticulum"]
        },
        "G6PD": {
            "function": "Glucose-6-phosphate dehydrogenase",
            "go_terms": ["dehydrogenase", "metabolic", "glucose metabolism", "NADPH production"]
        },
        "SERPINA1": {
            "function": "Alpha-1-antitrypsin",
            "go_terms": ["protease inhibitor", "regulator", "serine-type endopeptidase inhibitor"]
        },
    }

    for gene in gene_list:
        # Use realistic function description and GO terms if available
        if gene in gene_data:
            function = gene_data[gene]["function"]
            go_terms = gene_data[gene]["go_terms"]
        else:
            function = f"{gene} protein with unknown specific function"
            go_terms = ["protein binding"]

        classification = classify_gene_family(gene, function, go_terms)
        results[gene] = classification
        
        # Check against expected if we have it
        if gene in expected_families:
            expected = expected_families[gene]
            total_expected += 1
            if classification == expected:
                correct += 1
                status = "✅"
            else:
                status = "❌"
            print(f"{status} {gene:12} → {classification:20} (expected: {expected})")
        else:
            print(f"🔍 {gene:12} → {classification}")
    
    print()
    print("=" * 60)
    print(f"📊 ACCURACY ON EXPECTED CLASSIFICATIONS: {correct}/{total_expected} ({correct/total_expected*100:.1f}%)")
    
    # Count classifications
    classification_counts = {}
    for classification in results.values():
        classification_counts[classification] = classification_counts.get(classification, 0) + 1
    
    print()
    print("📈 CLASSIFICATION DISTRIBUTION:")
    for family, count in sorted(classification_counts.items()):
        print(f"   {family:25} → {count:2} genes")
    
    print()
    print("🎯 GENES THAT AVOIDED 'GENERAL' CLASSIFICATION:")
    non_general = [gene for gene, cls in results.items() if cls != "GENERAL"]
    print(f"   {len(non_general)}/{len(gene_list)} genes ({len(non_general)/len(gene_list)*100:.1f}%)")
    
    return results

if __name__ == "__main__":
    test_gene_mess()
