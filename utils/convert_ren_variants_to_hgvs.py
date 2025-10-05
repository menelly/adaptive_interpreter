#!/usr/bin/env python3
"""
Quick converter for Ren's personal variants TSV to HGVS format
Converts: GENE | TRANSCRIPT | VARIANT
To: HGVS | gnomAD frequency
"""

import sys

def convert_to_hgvs(gene, transcript, variant):
    """Convert Ren's format to HGVS"""
    import re

    # Clean up variant
    variant = variant.strip()

    # If it's in format like "c.1613C>T (p.Ser538Leu)", extract protein change
    if '(' in variant and 'p.' in variant:
        match = re.search(r'\(p\.([^)]+)\)', variant)
        if match:
            protein_change = f"p.{match.group(1)}"
        else:
            return None
    # If it's in format like "c.64C>T p.Arg22Trp" (space instead of parentheses)
    elif 'p.' in variant:
        match = re.search(r'p\.([A-Za-z0-9*]+)', variant)
        if match:
            protein_change = f"p.{match.group(1)}"
        else:
            return None
    # If variant already has p. prefix only
    elif variant.startswith('p.'):
        protein_change = variant
    # If it's in format like "c.386G>A" with no protein change, skip
    elif variant.startswith('c.') and 'p.' not in variant:
        return None
    else:
        # Assume it's already in protein format
        protein_change = f"p.{variant}" if not variant.startswith('p.') else variant

    # Build HGVS string - ALWAYS use parentheses format for parser
    if transcript and transcript.strip():
        hgvs = f"{transcript}({gene}):{protein_change}"
    else:
        # No transcript, use fake transcript format so parser can find gene in parentheses
        hgvs = f"NM_UNKNOWN({gene}):{protein_change}"

    return hgvs

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 convert_ren_variants_to_hgvs.py input.tsv output.tsv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write header
        outfile.write("HGVS\tgnomAD frequency\n")
        
        for line in infile:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 2:
                print(f"⚠️  Skipping malformed line: {line}")
                continue
            
            gene = parts[0].strip()
            
            # Handle different formats
            if len(parts) == 3:
                # GENE | TRANSCRIPT | VARIANT
                transcript = parts[1].strip()
                variant = parts[2].strip()
            elif len(parts) == 2:
                # GENE | VARIANT (no transcript)
                transcript = ""
                variant = parts[1].strip()
            else:
                print(f"⚠️  Unexpected format: {line}")
                continue
            
            hgvs = convert_to_hgvs(gene, transcript, variant)
            if hgvs:
                # Default frequency to 0 (will be looked up)
                outfile.write(f"{hgvs}\t0\n")
                print(f"✅ {gene} → {hgvs}")
            else:
                print(f"⏭️  Skipped {gene} {variant} (couldn't convert)")
    
    print(f"\n✅ Conversion complete! Output: {output_file}")

if __name__ == '__main__':
    main()

