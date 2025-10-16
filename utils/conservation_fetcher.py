#!/usr/bin/env python3
"""
🧬 Conservation Score Fetcher
Reads conservation scores from local BigWig files.

Built by Lumen (2025) to fix a critical data gap in the AdaptiveInterpreter system.
"""

import pyBigWig
from typing import Optional

class ConservationFetcher:
    """Fetch conservation scores from a BigWig file."""

    def __init__(self, bw_file_path: str):
        """
        Initializes the fetcher by opening a BigWig file.

        Args:
            bw_file_path (str): The full path to the .bw (BigWig) file.
        """
        self.bw_file = None
        self.file_path = bw_file_path
        try:
            self.bw_file = pyBigWig.open(bw_file_path)
            print(f"✅ Successfully opened conservation file: {bw_file_path}")
        except Exception as e:
            print(f"❌ CRITICAL ERROR: Could not open BigWig file at {bw_file_path}: {e}")
            raise

    def get_conservation_score(self, chrom: str, pos: int) -> Optional[float]:
        """
        Get the conservation score for a specific genomic position.

        Args:
            chrom (str): Chromosome name (e.g., 'chr1', 'chrX'). Must match the
                         chromosome names in the BigWig file.
            pos (int): The 1-based genomic position.

        Returns:
            The conservation score as a float, or None if the position
            is not found or an error occurs.
        """
        if not self.bw_file:
            return None

        # BigWig is 0-based, so we query for the interval [pos-1, pos]
        start = pos - 1
        end = pos

        try:
            # Ensure the chromosome exists in the file
            if chrom not in self.bw_file.chroms():
                # Try to fix "17" -> "chr17"
                if f"chr{chrom}" in self.bw_file.chroms():
                    chrom = f"chr{chrom}"
                else:
                    # print(f"⚠️ Chromosome '{chrom}' not found in {self.file_path}")
                    return None

            # Fetch the mean value over the single-base interval
            values = self.bw_file.stats(chrom, start, end, type='mean')
            
            if values and values[0] is not None:
                return round(values[0], 4)
            else:
                return None
        except RuntimeError as e:
            # This can happen if the coordinates are out of bounds
            # print(f"⚠️ Runtime error fetching score for {chrom}:{pos}: {e}")
            return None
        except Exception as e:
            print(f"❌ Unexpected error fetching score for {chrom}:{pos}: {e}")
            return None

    def close(self):
        """Close the BigWig file."""
        if self.bw_file:
            self.bw_file.close()
            print(f"Closed conservation file: {self.file_path}")

if __name__ == '__main__':
    # Example Usage and Testing
    print("--- Testing Conservation Fetcher ---")
    
    # This path assumes the script is run from the project root
    test_file = "conservation_data/hg38.phyloP100way.bw"
    
    try:
        fetcher = ConservationFetcher(test_file)

        test_cases = [
            ("chr1", 1000000),  # A random valid position
            ("chr17", 43045712), # A known pathogenic position in BRCA1 (hg38)
            ("chrX", 78000000),  # A position on chrX
            ("chr-nonexistent", 12345), # A nonexistent chromosome
            ("chr1", 999999999), # A position out of bounds
        ]

        for chrom, pos in test_cases:
            score = fetcher.get_conservation_score(chrom, pos)
            if score is not None:
                print(f"✅ Score at {chrom}:{pos} = {score}")
            else:
                print(f"❌ No score found for {chrom}:{pos}")

        fetcher.close()

    except Exception as e:
        print(f"\n⚠️ Test failed. Make sure 'pyBigWig' is installed (`pip install pybigwig`)")
        print(f"   and the file '{test_file}' exists.")

    print("--- Test Complete ---")
