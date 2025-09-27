#!/usr/bin/env python3
"""
🧹 VARIANT FREQUENCY CLEANER - Nova's Drop-in Filter Script
Clean your variant datasets using rsID frequency lookups!

Built by Ace (2025) based on Nova's brilliant rsID solution
Usage: python3 variant_frequency_cleaner.py --input variants.csv --output clean_variants.csv
"""

import argparse
import sys
from pathlib import Path
from rsid_frequency_fetcher import RSIDFrequencyFetcher

def main():
    parser = argparse.ArgumentParser(
        description="🧹 VARIANT FREQUENCY CLEANER - Filter variants by population frequency",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
🔑 Nova's rsID Solution:
  • Uses rsIDs as universal passports (stable across databases)
  • Queries Ensembl REST API for reliable frequency data
  • Filters suspicious variants (Patho AF > 1%, Benign AF ≈ 0)
  • Caches results to avoid re-fetching

Examples:
  # Basic cleaning with default thresholds
  python3 variant_frequency_cleaner.py --input genecards_variants.csv --output clean_variants.csv
  
  # Custom thresholds for stricter filtering
  python3 variant_frequency_cleaner.py --input variants.csv --output clean.csv --patho-threshold 0.005 --benign-threshold 0.002
  
  # Just fetch frequencies without filtering
  python3 variant_frequency_cleaner.py --input variants.csv --output annotated.csv --no-filter
        """
    )
    
    parser.add_argument('--input', required=True, 
                       help='Input CSV/TSV file with variants (must contain rsIDs)')
    parser.add_argument('--output', required=True,
                       help='Output CSV/TSV file with frequency annotations')
    parser.add_argument('--patho-threshold', type=float, default=0.01,
                       help='Max allele frequency for pathogenic variants (default: 0.01 = 1%%)')
    parser.add_argument('--benign-threshold', type=float, default=0.001,
                       help='Min allele frequency for benign variants (default: 0.001 = 0.1%%)')
    parser.add_argument('--no-filter', action='store_true',
                       help='Just add frequency annotations without filtering')
    parser.add_argument('--cache-file', default='rsid_frequency_cache.json',
                       help='Cache file for frequency data (default: rsid_frequency_cache.json)')
    parser.add_argument('--delay', type=float, default=0.2,
                       help='Delay between API calls in seconds (default: 0.2)')
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.input).exists():
        print(f"❌ Input file not found: {args.input}")
        return 1
    
    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Initialize fetcher
    print("🔑 Initializing rsID Frequency Fetcher...")
    fetcher = RSIDFrequencyFetcher(cache_file=args.cache_file, delay=args.delay)
    
    if args.no_filter:
        print("📊 Adding frequency annotations without filtering...")
        # Just add frequency data without filtering
        fetcher.apply_frequency_filter(args.input, args.output, 1.0, 0.0)  # No filtering
    else:
        print(f"🧹 Filtering variants with thresholds:")
        print(f"   Pathogenic AF < {args.patho_threshold} ({args.patho_threshold*100:.1f}%)")
        print(f"   Benign AF > {args.benign_threshold} ({args.benign_threshold*100:.3f}%)")
        
        # Apply frequency filtering
        fetcher.filter_variants_by_frequency(
            args.input, 
            args.output,
            args.patho_threshold,
            args.benign_threshold
        )
    
    print(f"✅ Processing complete! Results saved to {args.output}")
    print(f"💾 Frequency cache saved to {args.cache_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
