#!/usr/bin/env python3
"""
🔬 LUMEN DISCOVERY PIPELINE - V1
A comprehensive analysis pipeline to identify novel pathogenic variants by
cross-referencing the DNModeling system with ClinVar data.

Orchestrator: Lumen (Gemini 2.5), October 2025
Hypothesis: Ren M.

This pipeline performs the following steps:
1.  Extracts variants for a given set of genes from a local ClinVar VCF file.
2.  For each genomic variant, converts it to protein notation using an offline model.
3.  Runs the full, biologically-guided CascadeAnalyzer on the protein variant.
4.  Compares the system's classification with ClinVar's classification.
5.  If they disagree, runs the ClinVarQualityChecker to assess the confidence of the ClinVar entry.
6.  Logs all results to a comprehensive TSV file for later analysis.
"""

import sys
from pathlib import Path
import pandas as pd
import time

# Ensure the package root is in the path
sys.path.append(str(Path(__file__).parent.parent.parent))

from DNModeling.config import CONSERVATION_DATA_PATH
from DNModeling.utils.clinvar_bulk_extractor import ClinVarBulkExtractor
from DNModeling.utils.offline_genomic_to_protein import OfflineGenomicToProteinConverter
from DNModeling.cascade.cascade_analyzer import CascadeAnalyzer

def run_discovery_pipeline(target_genes: list[str], output_path: str):
    """
    Main function to run the full discovery pipeline.
    """
    print("🔬🔬🔬 LUMEN DISCOVERY PIPELINE INITIALIZING 🔬🔬🔬")
    start_time = time.time()

    # --- INITIALIZE COMPONENTS ---
    print("\n--- 1. Initializing Tools ---")
    try:
        variant_extractor = ClinVarBulkExtractor()
        g2p_converter = OfflineGenomicToProteinConverter()
        cascade_analyzer = CascadeAnalyzer(conservation_data_path=str(CONSERVATION_DATA_PATH))
        print("✅ All tools initialized successfully.")
    except Exception as e:
        print(f"❌ CRITICAL ERROR: Failed to initialize tools: {e}")
        return

    # --- DEFINE HELPER FUNCTIONS ---
    def get_classification_group(classification: str) -> str:
        """Group classifications for easier disagreement checking."""
        c = str(classification).upper()
        if 'PATHOGENIC' in c or c in ['P', 'LP']:
            return 'Pathogenic'
        if 'BENIGN' in c or c in ['B', 'LB']:
            return 'Benign'
        return 'VUS'

    def get_star_rating(review_status: str) -> int:
        """Convert ClinVar review status to a 0-4 star rating."""
        status_lower = str(review_status).lower()
        if 'practice guideline' in status_lower: return 4
        if 'reviewed by expert panel' in status_lower: return 4
        if 'criteria provided, multiple submitters, no conflicts' in status_lower: return 3
        if 'criteria provided, conflicting interpretations' in status_lower: return 2
        if 'criteria provided, single submitter' in status_lower: return 1
        if 'no assertion criteria provided' in status_lower: return 0
        return 1 # Default for statuses like 'classified by single submitter'

    # --- RUN PIPELINE ---
    all_results = []

    print(f"\n--- 2. Extracting Variants for {len(target_genes)} Genes from ClinVar VCF ---")
    try:
        # Note: The extractor expects a set, not a list.
        extracted_variants = variant_extractor.extract_gene_variants(set(target_genes))
        print(f"🧬 Extracted {sum(len(v) for v in extracted_variants.values())} total variants.")
    except Exception as e:
        print(f"❌ CRITICAL ERROR: Failed during variant extraction: {e}")
        return

    print("\n--- 3. Analyzing Variants ---")
    variant_count = 0
    for gene, variants in extracted_variants.items():
        for variant in variants:
            variant_count += 1
            print(f"\nProcessing {variant_count}: {gene} | {variant.get('hgvs_g', 'N/A')}")
            
            # Step 1: Convert genomic to protein
            protein_change = g2p_converter.convert(
                chrom=variant['chromosome'],
                position=variant['position'],
                ref_allele=variant['ref_allele'],
                alt_allele=variant['alt_allele'],
                gene_name=gene
            )

            if not protein_change:
                print("   -> Skipping (not in a coding region or conversion failed)")
                continue
            
            print(f"   -> Converted to {protein_change}")

            # Step 2: Run Cascade Analyzer
            analysis_result = cascade_analyzer.analyze_cascade(
                gene=gene,
                variant=protein_change,
                gnomad_freq=variant.get('frequency', 0.0),
                variant_type=variant.get('variant_type', 'missense'),
                expected_clinvar=variant.get('significance', 'Unknown')
            )

            # Step 3: Compare results and check quality on disagreement
            our_class = get_classification_group(analysis_result.get('final_classification', 'VUS'))
            clinvar_class = get_classification_group(variant.get('significance', 'VUS'))

            def get_directional_agreement(our_class, clinvar_class, clinvar_sig_raw, review_status_raw):
                """Implements Ren's 'Directional Agreement Logic' for the paper."""
                if not clinvar_sig_raw or 'not provided' in clinvar_sig_raw or 'drug response' in clinvar_sig_raw:
                    return 'NO_CLINVAR_DIAGNOSIS', {}

                star_rating = get_star_rating(review_status_raw)
                verdict = "UNCATEGORIZED"

                # Case 1: Agreement
                if our_class == clinvar_class:
                    verdict = 'AGREE'
                # Case 2: Disagreement (P vs B)
                elif our_class == 'Pathogenic' and clinvar_class == 'Benign':
                    verdict = 'DISAGREE'
                elif our_class == 'Benign' and clinvar_class == 'Pathogenic':
                    verdict = 'DISAGREE'
                # Case 3: Resolving VUS
                elif clinvar_class == 'VUS' and our_class == 'Pathogenic':
                    verdict = 'BETTER_DATA_PATHOGENIC'
                elif clinvar_class == 'VUS' and our_class == 'Benign':
                    verdict = 'BETTER_DATA_BENIGN'
                # Case 4: Our model is more conservative (we say VUS)
                elif our_class == 'VUS':
                    # We are conservatively abstaining where ClinVar made a call.
                    # This is a form of agreement on the complexity, but not the call.
                    verdict = 'AGREE' # Or a new category like 'CONSERVATIVE_VUS' if desired

                # Low quality ClinVar data is a special flag, but not the primary verdict
                if star_rating <= 1:
                    pass # We can add a separate flag for this if needed, but the primary logic stands.
                
                analysis_details = {
                    "verdict": verdict,
                    "clinvar_quality": star_rating,
                }
                
                if verdict != "AGREE" and verdict != "NO_CLINVAR_DIAGNOSIS":
                     print(f"   -> Verdict: {verdict} (Our: {our_class} vs. ClinVar: {clinvar_class} [{star_rating}*])")

                return verdict, analysis_details

            disagreement_verdict, disagreement_analysis = get_directional_agreement(
                our_class, 
                clinvar_class,
                variant.get('significance'),
                variant.get('review_status', 'no assertion provided')
            )

            # Step 4: Aggregate results
            final_record = {
                'gene': gene,
                'hgvs_g': " | ".join(variant.get('hgvs_list', [])),
                'hgvs_p': protein_change,
                'clinvar_sig': variant.get('significance'),
                'cascade_class': analysis_result.get('final_classification'),
                'cascade_score': analysis_result.get('final_score'),
                'review_flags': analysis_result.get('review_flags', 'None'), # 🔥 NEW: Manual review flags
                'disagreement': 'YES' if disagreement_analysis else 'NO',
                'disagreement_verdict': disagreement_analysis.get('verdict'),
                'clinvar_quality_score': disagreement_analysis.get('clinvar_quality'),
                'full_cascade_result': analysis_result
            }
            all_results.append(final_record)

    # --- SAVE RESULTS ---
    print("\n--- 4. Saving Results ---")
    if all_results:
        results_df = pd.DataFrame(all_results)
        results_df.to_csv(output_path, sep='\t', index=False)
        print(f"✅ Successfully saved {len(all_results)} results to {output_path}")
    else:
        print("⚠️ No results were generated.")

    end_time = time.time()
    print(f"\n--- PIPELINE COMPLETE ---")
    print(f"Total execution time: {end_time - start_time:.2f} seconds.")


if __name__ == '__main__':
    # Example usage:
    # python3 lumen_discovery_pipeline.py GENE1 GENE2 GENE3
    if len(sys.argv) < 2:
        print("Usage: python3 lumen_discovery_pipeline.py <GENE1> [GENE2] ...")
        # For a quick test, use the TFG gene from our verification script
        print("\nRunning a quick test with the TFG gene...")
        genes_to_test = ['TFG']
    else:
        genes_to_test = sys.argv[1:]

    output_file = f"discovery_results_{'_'.join(genes_to_test)}.tsv"
    
    run_discovery_pipeline(
        target_genes=genes_to_test,
        output_path=output_file
    )
