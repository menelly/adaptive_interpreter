#!/usr/bin/env python3
"""
AlphaFold Sequence Extractor for DN Analysis

Extracts protein sequences from local AlphaFold structure files
for use with the Nova DN analyzer.

Authors: Nova & Ace (2025)
"""

import gzip
import os
from typing import Optional, Dict


class AlphaFoldSequenceExtractor:
    """Extract protein sequences from local AlphaFold PDB files"""
    
    def __init__(self, alphafold_path: str = "/mnt/Arcana/alphafold_human/structures/"):
        self.alphafold_path = alphafold_path
        
        # 3-letter to 1-letter amino acid code mapping
        self.aa_map = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
    
    def get_sequence(self, uniprot_id: str) -> Optional[str]:
        """
        Extract protein sequence from AlphaFold PDB file
        
        Args:
            uniprot_id: UniProt ID (e.g., 'P04637' for TP53)
            
        Returns:
            Protein sequence string, or None if not found
        """
        pdb_path = os.path.join(self.alphafold_path, f"AF-{uniprot_id}-F1-model_v4.pdb.gz")
        
        if not os.path.exists(pdb_path):
            print(f"AlphaFold structure not found: {pdb_path}")
            return None
        
        try:
            with gzip.open(pdb_path, 'rt') as f:
                lines = f.readlines()
            
            sequence = ''
            prev_resnum = -1
            
            # Extract sequence from ATOM records (CA atoms only)
            for line in lines:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    resnum = int(line[22:26])
                    if resnum != prev_resnum:
                        aa_code = line[17:20].strip()
                        if aa_code in self.aa_map:
                            sequence += self.aa_map[aa_code]
                        prev_resnum = resnum
            
            return sequence if sequence else None
            
        except Exception as e:
            print(f"Error extracting sequence from {pdb_path}: {e}")
            return None
    
    def save_fasta(self, uniprot_id: str, output_path: str, gene_name: str = None) -> bool:
        """
        Extract sequence and save as FASTA file
        
        Args:
            uniprot_id: UniProt ID
            output_path: Path to save FASTA file
            gene_name: Optional gene name for FASTA header
            
        Returns:
            True if successful, False otherwise
        """
        sequence = self.get_sequence(uniprot_id)
        if not sequence:
            return False
        
        try:
            with open(output_path, 'w') as f:
                header = f">{gene_name or uniprot_id} {uniprot_id}"
                f.write(f"{header}\n")
                f.write(f"{sequence}\n")
            
            print(f"Saved {len(sequence)} residue sequence to {output_path}")
            return True
            
        except Exception as e:
            print(f"Error saving FASTA file: {e}")
            return False
    
    def verify_variant(self, uniprot_id: str, position: int, expected_aa: str) -> Dict:
        """
        Verify that a variant position matches the expected amino acid
        
        Args:
            uniprot_id: UniProt ID
            position: 1-based position
            expected_aa: Expected amino acid (single letter)
            
        Returns:
            Dict with verification results
        """
        sequence = self.get_sequence(uniprot_id)
        if not sequence:
            return {"error": "Could not extract sequence"}
        
        if position < 1 or position > len(sequence):
            return {"error": f"Position {position} out of range (1-{len(sequence)})"}
        
        actual_aa = sequence[position - 1]  # Convert to 0-based
        
        return {
            "position": position,
            "expected": expected_aa,
            "actual": actual_aa,
            "match": actual_aa == expected_aa,
            "sequence_length": len(sequence),
            "context": sequence[max(0, position-6):position+5]  # 5 residues each side
        }


def main():
    """Test the AlphaFold sequence extractor"""
    extractor = AlphaFoldSequenceExtractor()
    
    # Test proteins
    test_cases = [
        ("TP53", "P04637", 273, "R"),
        ("ATP5F1A", "P25705", 130, "I"),
        ("TFG", "Q92734", 22, "R"),
        ("PYGL", "P06737", 634, "D"),
        ("MYO7A", "Q13402", 220, "H"),
        ("FKRP", "Q9H9S5", 276, "L"),
        ("DLD", "P09622", 34, "T"),
        ("ACMSD", "Q8TDX5", 175, "P")
    ]
    
    print("🧬 Testing AlphaFold Sequence Extractor 🧬")
    print("=" * 50)
    
    for gene, uniprot_id, pos, expected_aa in test_cases:
        print(f"\n{gene} ({uniprot_id}):")
        
        # Verify variant
        result = extractor.verify_variant(uniprot_id, pos, expected_aa)
        if "error" in result:
            print(f"  ❌ {result['error']}")
        else:
            match_symbol = "✅" if result["match"] else "❌"
            print(f"  {match_symbol} Position {pos}: {result['actual']} (expected {expected_aa})")
            print(f"     Context: {result['context']}")
            print(f"     Length: {result['sequence_length']} residues")
        
        # Save FASTA
        fasta_path = f"/tmp/{gene.lower()}_alphafold.fasta"
        if extractor.save_fasta(uniprot_id, fasta_path, gene):
            print(f"  💾 Saved to {fasta_path}")


if __name__ == "__main__":
    main()
