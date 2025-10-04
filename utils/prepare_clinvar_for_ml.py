#!/usr/bin/env python3
"""
🧬 PREPARE CLINVAR DATA FOR ML TRAINING
Convert ClinVar genomic variants to protein notation for ML training

Takes the output from clinvar_bulk_extractor.py and converts genomic
coordinates to protein changes using our offline converter.

Usage:
    python prepare_clinvar_for_ml.py --gene FKRP --input tests/clinvar_extracted/ --output learning/
"""

import argparse
import pandas as pd
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))
from utils.offline_genomic_to_protein import OfflineGenomicToProteinConverter


def parse_genomic_hgvs(hgvs: str):
    """Parse genomic HGVS to extract chromosome, position, ref, alt"""
    import re
    
    # Pattern: NC_000019.10:g.46755451A>C
    match = re.search(r'NC_0+(\d+)(?:\.\d+)?:g\.(\d+)([ATCG])>([ATCG])', hgvs)
    if match:
        chrom_num = match.group(1)
        position = int(match.group(2))
        ref_allele = match.group(3)
        alt_allele = match.group(4)
        
        # Convert NC number to chr
        if chrom_num == "23":
            chrom = "chrX"
        elif chrom_num == "24":
            chrom = "chrY"
        elif chrom_num == "12920":
            chrom = "chrM"
        else:
            chrom = f"chr{chrom_num}"
        
        return chrom, position, ref_allele, alt_allele
    
    return None, None, None, None


def convert_clinvar_to_ml_format(gene: str, input_dir: Path, output_dir: Path, converter: OfflineGenomicToProteinConverter):
    """Convert ClinVar TSV files to ML training format"""
    
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"🧬 Converting ClinVar data for {gene}")
    print(f"   Input: {input_dir}")
    print(f"   Output: {output_dir}")
    
    for significance in ['pathogenic', 'benign']:
        input_file = input_dir / f"{gene.lower()}_{significance}.tsv"
        
        if not input_file.exists():
            print(f"   ⚠️ File not found: {input_file}")
            continue
        
        print(f"\n📖 Processing {significance} variants...")
        df = pd.read_csv(input_file, sep='\t')
        print(f"   Found {len(df)} variants")
        
        converted_variants = []
        failed_conversions = 0
        
        for idx, row in df.iterrows():
            if idx % 50 == 0 and idx > 0:
                print(f"      Processed {idx}/{len(df)} variants...")
            
            genomic_hgvs = row['HGVS']
            
            # Parse genomic coordinates
            chrom, position, ref_allele, alt_allele = parse_genomic_hgvs(genomic_hgvs)
            
            if not chrom:
                failed_conversions += 1
                continue
            
            # Convert to protein notation
            protein_change = converter.convert(chrom, position, ref_allele, alt_allele, gene)
            
            if not protein_change:
                failed_conversions += 1
                continue
            
            # Create HGVS-like format for ML trainer
            # Format: NM_XXXXXX.X(GENE):c.XXX (p.RefPosAlt)
            ml_hgvs = f"NM_000000.0({gene}):c.{position}A>C ({protein_change})"
            
            converted_variants.append({
                'HGVS': ml_hgvs,
                'dbSNP': row.get('dbSNP', ''),
                'gnomAD frequency': row.get('gnomAD frequency', 0.0),
                'ClinVar_significance': row.get('ClinVar_significance', significance),
                'Review_status': row.get('Review_status', 'unknown'),
                'genomic_hgvs': genomic_hgvs,
                'protein_change': protein_change
            })
        
        print(f"   ✅ Converted {len(converted_variants)} variants")
        print(f"   ⚠️ Failed to convert {failed_conversions} variants")
        
        if converted_variants:
            # Save to output
            output_file = output_dir / f"{gene.lower()}_{significance}.tsv"
            output_df = pd.DataFrame(converted_variants)
            
            # Keep only columns needed for ML training
            ml_columns = ['HGVS', 'dbSNP', 'gnomAD frequency', 'ClinVar_significance', 'Review_status']
            output_df[ml_columns].to_csv(output_file, sep='\t', index=False)
            
            print(f"   💾 Saved to {output_file}")
            
            # Also save a detailed version with genomic info
            detail_file = output_dir / f"{gene.lower()}_{significance}_detailed.tsv"
            output_df.to_csv(detail_file, sep='\t', index=False)
            print(f"   💾 Saved detailed version to {detail_file}")


def main():
    parser = argparse.ArgumentParser(description="Prepare ClinVar data for ML training")
    parser.add_argument('--gene', required=True, help='Gene symbol (e.g., FKRP)')
    parser.add_argument('--input', required=True, help='Input directory with ClinVar TSV files')
    parser.add_argument('--output', default='learning/', help='Output directory for ML training files')
    
    args = parser.parse_args()
    
    print("🔥💜 CLINVAR TO ML CONVERTER")
    print("=" * 60)
    
    # Initialize converter
    print("\n🔨 Initializing offline genomic-to-protein converter...")
    converter = OfflineGenomicToProteinConverter()
    
    # Convert files
    print("\n🧬 Converting ClinVar data...")
    convert_clinvar_to_ml_format(
        gene=args.gene,
        input_dir=Path(args.input),
        output_dir=Path(args.output),
        converter=converter
    )
    
    print("\n🎉 CONVERSION COMPLETE!")
    print(f"   Ready for ML training with: python ml_training/train_families.py")


if __name__ == '__main__':
    main()

