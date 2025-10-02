#!/usr/bin/env python3
"""
Sequence Manager for DN Analysis 🧬📁
Manages both UniProt and AlphaFold sequences with persistent FASTA caching

Features:
- Downloads UniProt sequences once, caches forever
- Falls back to full UniProt when AlphaFold is truncated
- Organized folder structure for GitHub sharing
- No repeated API calls (UniProt-friendly!)

Authors: Ace & Nova (2025)
"""

import os
import requests
import time
from pathlib import Path
from typing import Optional, Dict, Tuple
from .alphafold_sequence import AlphaFoldSequenceExtractor
from .universal_context import UniversalContext


class SequenceManager:
    """Manages protein sequences with persistent FASTA caching"""
    
    def __init__(self, base_path: str = "sequences"):
        self.base_path = Path(base_path)
        self.uniprot_path = self.base_path / "uniprot"
        self.alphafold_path = self.base_path / "alphafold"
        
        # Create directories
        self.uniprot_path.mkdir(parents=True, exist_ok=True)
        self.alphafold_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize extractors
        self.alphafold_extractor = AlphaFoldSequenceExtractor()
        self.universal_context = UniversalContext()
        
        # Rate limiting for UniProt (be nice!)
        self.last_uniprot_request = 0
        self.min_request_interval = 1.0  # 1 second between requests
    
    def get_uniprot_fasta_path(self, gene_name: str, uniprot_id: str) -> Path:
        """Get path for UniProt FASTA file"""
        return self.uniprot_path / f"{gene_name}_{uniprot_id}.fasta"
    
    def get_alphafold_fasta_path(self, gene_name: str, uniprot_id: str) -> Path:
        """Get path for AlphaFold FASTA file"""
        return self.alphafold_path / f"{gene_name}_{uniprot_id}_AF.fasta"
    
    def load_fasta_sequence(self, fasta_path: Path) -> Optional[str]:
        """Load sequence from FASTA file"""
        if not fasta_path.exists():
            return None
        
        try:
            with open(fasta_path, 'r') as f:
                lines = f.readlines()
            
            # Skip header line, join sequence lines
            sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
            return sequence if sequence else None
            
        except Exception as e:
            print(f"⚠️ Error reading FASTA {fasta_path}: {e}")
            return None
    
    def save_fasta_sequence(self, sequence: str, fasta_path: Path, 
                           gene_name: str, uniprot_id: str, source: str = "UniProt") -> bool:
        """Save sequence to FASTA file"""
        try:
            with open(fasta_path, 'w') as f:
                f.write(f">{gene_name} {uniprot_id} | {source} | Length: {len(sequence)}\n")
                # Write sequence in 80-character lines (standard FASTA format)
                for i in range(0, len(sequence), 80):
                    f.write(f"{sequence[i:i+80]}\n")
            
            print(f"💾 Saved {len(sequence)} residue {source} sequence to {fasta_path}")
            return True
            
        except Exception as e:
            print(f"❌ Error saving FASTA {fasta_path}: {e}")
            return False
    
    def download_uniprot_sequence(self, uniprot_id: str) -> Optional[str]:
        """Download sequence from UniProt API with rate limiting"""
        
        # Rate limiting - be nice to UniProt!
        current_time = time.time()
        time_since_last = current_time - self.last_uniprot_request
        if time_since_last < self.min_request_interval:
            sleep_time = self.min_request_interval - time_since_last
            print(f"⏱️ Rate limiting: sleeping {sleep_time:.1f}s")
            time.sleep(sleep_time)
        
        try:
            url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
            print(f"🌐 Downloading {uniprot_id} from UniProt...")
            
            response = requests.get(url, timeout=30)
            self.last_uniprot_request = time.time()
            
            if response.status_code == 200:
                # Parse FASTA content
                lines = response.text.strip().split('\n')
                sequence = ''.join(line for line in lines if not line.startswith('>'))
                return sequence if sequence else None
            else:
                print(f"❌ UniProt request failed: {response.status_code}")
                return None
                
        except Exception as e:
            print(f"❌ Error downloading from UniProt: {e}")
            return None
    
    def get_uniprot_sequence(self, gene_name: str, uniprot_id: str) -> Optional[str]:
        """Get UniProt sequence (cached or downloaded)"""
        
        # Check if we already have it cached
        fasta_path = self.get_uniprot_fasta_path(gene_name, uniprot_id)
        sequence = self.load_fasta_sequence(fasta_path)
        
        if sequence:
            print(f"✅ Using cached UniProt sequence for {gene_name} ({len(sequence)} residues)")
            return sequence
        
        # Try to get from universal context (which might have it cached)
        context = self.universal_context.get_context_for_protein(gene_name, uniprot_id)
        if 'sequence' in context and context['sequence']:
            sequence = context['sequence']
            print(f"✅ Got UniProt sequence from context cache ({len(sequence)} residues)")
            
            # Save to FASTA for future use
            self.save_fasta_sequence(sequence, fasta_path, gene_name, uniprot_id, "UniProt")
            return sequence
        
        # Last resort: download directly from UniProt
        sequence = self.download_uniprot_sequence(uniprot_id)
        if sequence:
            print(f"✅ Downloaded UniProt sequence ({len(sequence)} residues)")
            self.save_fasta_sequence(sequence, fasta_path, gene_name, uniprot_id, "UniProt")
            return sequence
        
        print(f"❌ Could not get UniProt sequence for {gene_name} ({uniprot_id})")
        return None
    
    def get_alphafold_sequence(self, gene_name: str, uniprot_id: str) -> Optional[str]:
        """Get AlphaFold sequence (cached or extracted)"""
        
        # Check if we already have it cached
        fasta_path = self.get_alphafold_fasta_path(gene_name, uniprot_id)
        sequence = self.load_fasta_sequence(fasta_path)
        
        if sequence:
            print(f"✅ Using cached AlphaFold sequence for {gene_name} ({len(sequence)} residues)")
            return sequence
        
        # Extract from AlphaFold PDB
        sequence = self.alphafold_extractor.get_sequence(uniprot_id)
        if sequence:
            print(f"✅ Extracted AlphaFold sequence ({len(sequence)} residues)")
            self.save_fasta_sequence(sequence, fasta_path, gene_name, uniprot_id, "AlphaFold")
            return sequence
        
        print(f"❌ Could not get AlphaFold sequence for {gene_name} ({uniprot_id})")
        return None
    
    def get_best_sequence(self, gene_name: str, uniprot_id: str, 
                         variant_position: Optional[int] = None) -> Tuple[str, str, str]:
        """
        Get the best sequence for analysis
        
        Returns:
            (sequence, source, temp_fasta_path)
        """
        
        # Get both sequences
        uniprot_seq = self.get_uniprot_sequence(gene_name, uniprot_id)
        alphafold_seq = self.get_alphafold_sequence(gene_name, uniprot_id)
        
        # Decide which to use
        # IMPORTANT: Always use UniProt sequence for residue indexing to avoid numbering drift.
        # AlphaFold sequence is used for structure/context, not for AA-indexing checks.
        if uniprot_seq:
            chosen_seq = uniprot_seq
            source = "UniProt"
            if variant_position and alphafold_seq:
                print(f"🧬 Using UniProt sequence ({len(chosen_seq)} residues); AlphaFold available for structure")
            else:
                print(f"🧬 Using UniProt sequence ({len(chosen_seq)} residues)")
        elif alphafold_seq:
            # Fallback to AlphaFold if UniProt failed
            chosen_seq = alphafold_seq
            source = "AlphaFold"
            print(f"🧬 Fallback to AlphaFold sequence ({len(chosen_seq)} residues)")
        else:
            raise ValueError(f"Could not get any sequence for {gene_name} ({uniprot_id})")

        # Create temporary FASTA file for analyzer
        import tempfile
        temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        temp_fasta.write(f">{gene_name} {uniprot_id} | {source}\n")
        temp_fasta.write(f"{chosen_seq}\n")
        temp_fasta.close()
        
        return chosen_seq, source, temp_fasta.name
    
    def verify_variant_position(self, gene_name: str, uniprot_id: str, 
                               position: int, expected_aa: str) -> Dict:
        """Verify variant position against sequence"""
        
        uniprot_seq = self.get_uniprot_sequence(gene_name, uniprot_id)
        if not uniprot_seq:
            return {"error": "Could not get sequence"}
        
        if position < 1 or position > len(uniprot_seq):
            return {"error": f"Position {position} out of range (1-{len(uniprot_seq)})"}
        
        actual_aa = uniprot_seq[position - 1]  # Convert to 0-based
        
        return {
            "gene": gene_name,
            "uniprot_id": uniprot_id,
            "position": position,
            "expected": expected_aa,
            "actual": actual_aa,
            "match": actual_aa == expected_aa,
            "sequence_length": len(uniprot_seq),
            "context": uniprot_seq[max(0, position-6):position+5]  # 5 residues each side
        }


def test_sequence_manager():
    """Test the sequence manager"""
    manager = SequenceManager()
    
    # Test with RYR1 - should show the truncation issue
    print("🧬 Testing RYR1 sequence management")
    print("=" * 50)
    
    gene = "RYR1"
    uniprot_id = "P21817"
    variant_pos = 4842  # This should be beyond AlphaFold coverage
    
    try:
        sequence, source, temp_path = manager.get_best_sequence(gene, uniprot_id, variant_pos)
        print(f"✅ Got {len(sequence)} residue sequence from {source}")
        print(f"📁 Temp FASTA: {temp_path}")
        
        # Verify the problematic variant
        verification = manager.verify_variant_position(gene, uniprot_id, variant_pos, "V")
        if verification.get("match"):
            print(f"✅ Variant p.V{variant_pos}M verified: {verification['context']}")
        else:
            print(f"❌ Variant mismatch: expected V, got {verification.get('actual', 'unknown')}")
        
        # Cleanup
        os.unlink(temp_path)
        
    except Exception as e:
        print(f"❌ Test failed: {e}")


if __name__ == "__main__":
    test_sequence_manager()
