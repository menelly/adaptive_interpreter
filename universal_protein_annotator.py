#!/usr/bin/env python3
"""
Universal Protein Annotator - NO MORE HARDCODING!
Automatically extracts protein features from UniProt, Pfam, and sequence analysis
"""

import requests
import re
import json
from typing import Dict, List, Tuple, Optional
import time

class UniversalProteinAnnotator:
    def __init__(self):
        self.uniprot_base = "https://rest.uniprot.org"
        self.pfam_base = "https://pfam.xfam.org"
        
    def get_uniprot_features(self, uniprot_id: str) -> Dict:
        """Extract all features from UniProt API"""
        url = f"{self.uniprot_base}/uniprotkb/{uniprot_id}.json"
        
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            features = {
                "uniprot_id": uniprot_id,
                "sequence": data.get("sequence", {}).get("value", ""),
                "length": data.get("sequence", {}).get("length", 0),
                "gene_name": self._extract_gene_name(data),
                "function": self._extract_function(data),
                "domains": [],
                "active_sites": [],
                "binding_sites": [],
                "transmembrane": [],
                "signal_peptide": [],
                "coiled_coil": [],
                "disulfide_bonds": [],
                "glycosylation": [],
                "phosphorylation": []
            }

            # Debug: print function extraction
            print(f"Extracted function for {uniprot_id}: {features['function'][:100] if features['function'] else 'NONE'}...")
            
            # Extract features from UniProt annotations
            for feature in data.get("features", []):
                self._parse_uniprot_feature(feature, features)
                
            return features
            
        except Exception as e:
            print(f"UniProt API failed for {uniprot_id}: {e}")
            return {"uniprot_id": uniprot_id, "error": str(e)}
    
    def _extract_gene_name(self, data: Dict) -> str:
        """Extract primary gene name"""
        genes = data.get("genes", [])
        if genes and "geneName" in genes[0]:
            return genes[0]["geneName"]["value"]
        return ""
    
    def _extract_function(self, data: Dict) -> str:
        """Extract function description"""
        comments = data.get("comments", [])
        for comment in comments:
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    return texts[0].get("value", "")
        return ""
    
    def _parse_uniprot_feature(self, feature: Dict, features: Dict):
        """Parse individual UniProt feature"""
        feature_type = feature.get("type")
        location = feature.get("location", {})
        start = location.get("start", {}).get("value")
        end = location.get("end", {}).get("value")
        description = feature.get("description", "")
        
        if not start:
            return
            
        # Map UniProt feature types to our categories
        if feature_type == "DOMAIN":
            features["domains"].append({
                "start": start, "end": end or start,
                "type": "domain", "description": description
            })
        elif feature_type == "ACT_SITE":
            features["active_sites"].append(start)
        elif feature_type == "BINDING":
            features["binding_sites"].append({
                "position": start, "description": description
            })
        elif feature_type == "TRANSMEM":
            features["transmembrane"].append({
                "start": start, "end": end or start
            })
        elif feature_type == "SIGNAL":
            features["signal_peptide"].append({
                "start": start, "end": end or start
            })
        elif feature_type == "COILED":
            features["coiled_coil"].append({
                "start": start, "end": end or start
            })
        elif feature_type == "DISULFID":
            # Parse disulfide bond positions
            if ";" in description:
                positions = re.findall(r'\d+', description)
                if len(positions) >= 2:
                    features["disulfide_bonds"].append([int(positions[0]), int(positions[1])])
        elif feature_type == "CARBOHYD":
            features["glycosylation"].append(start)
        elif feature_type == "MOD_RES" and "phospho" in description.lower():
            features["phosphorylation"].append(start)
    
    def detect_sequence_motifs(self, sequence: str) -> Dict:
        """Detect common sequence motifs"""
        motifs = {
            "collagen_repeats": [],
            "walker_a": [],
            "walker_b": [],
            "leucine_zipper": [],
            "zinc_finger": [],
            "helix_turn_helix": []
        }
        
        # Collagen Gly-X-Y repeats - STRICT detection
        collagen_regions = []
        i = 0
        while i < len(sequence) - 20:  # Need at least 7 triplets
            if sequence[i] == 'G':
                # Look for consecutive Gly-X-Y pattern
                triplet_count = 0
                j = i
                while j < len(sequence) - 2:
                    if sequence[j] == 'G':
                        triplet_count += 1
                        j += 3
                    else:
                        break

                # Only count as collagen if we have 7+ consecutive triplets
                if triplet_count >= 7:
                    collagen_regions.append({
                        "start": i + 1,  # 1-based
                        "end": j
                    })
                    i = j
                else:
                    i += 1
            else:
                i += 1

        if collagen_regions:
            motifs["collagen_repeats"] = collagen_regions
        
        # Walker A motif: [AG]XXXXGK[ST]
        walker_a_pattern = r'[AG].{4}GK[ST]'
        for match in re.finditer(walker_a_pattern, sequence):
            motifs["walker_a"].append({
                "start": match.start() + 1,
                "end": match.end()
            })
        
        # Walker B motif: [RK]X{2}[DE]
        walker_b_pattern = r'[RK].{2}[DE]'
        for match in re.finditer(walker_b_pattern, sequence):
            motifs["walker_b"].append({
                "start": match.start() + 1,
                "end": match.end()
            })
        
        # Leucine zipper: L-X6-L-X6-L pattern
        leu_zip_pattern = r'L.{6}L.{6}L'
        for match in re.finditer(leu_zip_pattern, sequence):
            motifs["leucine_zipper"].append({
                "start": match.start() + 1,
                "end": match.end()
            })
        
        return motifs
    
    def predict_functional_sites(self, sequence: str, features: Dict) -> Dict:
        """Predict functional sites from sequence and known features"""
        predictions = {
            "dna_contact_sites": [],
            "interface_likelihood": [],
            "flexible_loops": [],
            "critical_residues": []
        }
        
        # DNA-binding prediction: look for basic residues in domains
        for domain in features.get("domains", []):
            if any(term in domain.get("description", "").lower() 
                   for term in ["dna", "bind", "helix", "zinc finger"]):
                # Find basic residues in this domain
                start, end = domain["start"], domain["end"]
                for i in range(start-1, min(end, len(sequence))):
                    if sequence[i] in "RKH":
                        predictions["dna_contact_sites"].append(i + 1)
        
        # Interface prediction: hydrophobic patches
        window_size = 7
        for i in range(len(sequence) - window_size):
            window = sequence[i:i+window_size]
            hydrophobic_count = sum(1 for aa in window if aa in "AILMVFWY")
            if hydrophobic_count >= 4:  # Hydrophobic patch
                predictions["interface_likelihood"].extend(range(i+1, i+window_size+1))
        
        # Flexible loops: regions between domains
        domains = sorted(features.get("domains", []), key=lambda x: x["start"])
        for i in range(len(domains) - 1):
            gap_start = domains[i]["end"] + 1
            gap_end = domains[i+1]["start"] - 1
            if gap_end - gap_start > 5:  # Significant gap
                predictions["flexible_loops"].extend(range(gap_start, gap_end + 1))
        
        return predictions
    
    def annotate_protein(self, gene_name: str, uniprot_id: Optional[str] = None) -> Dict:
        """Complete protein annotation pipeline"""
        if not uniprot_id:
            # Try to find UniProt ID from gene name
            uniprot_id = self._find_uniprot_id(gene_name)
            if not uniprot_id:
                return {"error": f"Could not find UniProt ID for {gene_name}"}
        
        print(f"Annotating {gene_name} ({uniprot_id})...")
        
        # Get UniProt features
        features = self.get_uniprot_features(uniprot_id)
        if "error" in features:
            return features
        
        sequence = features.get("sequence", "")
        if not sequence:
            return {"error": "No sequence found"}
        
        # Detect sequence motifs
        motifs = self.detect_sequence_motifs(sequence)
        features.update(motifs)
        
        # Predict functional sites
        predictions = self.predict_functional_sites(sequence, features)
        features.update(predictions)
        
        # Clean up and format for Nova's context system
        nova_format = self._format_for_nova(gene_name, features)
        
        return nova_format
    
    def _find_uniprot_id(self, gene_name: str) -> Optional[str]:
        """Find UniProt ID from gene name"""
        url = f"{self.uniprot_base}/uniprotkb/search"
        params = {
            "query": f"gene:{gene_name} AND organism_id:9606",
            "format": "json",
            "size": 1
        }
        
        try:
            response = requests.get(url, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            results = data.get("results", [])
            if results:
                return results[0].get("primaryAccession")
                
        except Exception as e:
            print(f"UniProt search failed for {gene_name}: {e}")
        
        return None
    
    def _format_for_nova(self, gene_name: str, features: Dict) -> Dict:
        """Format annotations for Nova's context system"""
        nova_format = {
            "uniprot_id": features.get("uniprot_id"),
            "sequence_length": features.get("length", 0),
            "function": features.get("function", ""),
            "sequence": features.get("sequence", ""),
            "domains": features.get("domains", [])
        }
        
        # Map to Nova's expected fields
        if features.get("dna_contact_sites"):
            nova_format["dna_contact_sites"] = features["dna_contact_sites"]
        
        if features.get("active_sites"):
            nova_format["known_active_or_binding_sites"] = features["active_sites"]
        
        if features.get("transmembrane"):
            tm = features["transmembrane"]
            if tm:
                nova_format["transmembrane_domain"] = [tm[0]["start"], tm[0]["end"]]
        
        if features.get("collagen_repeats"):
            cr = features["collagen_repeats"]
            if cr:
                nova_format["collagen_repeats"] = cr[0]
        
        if features.get("disulfide_bonds"):
            nova_format["disulfide_pairs"] = features["disulfide_bonds"]
        
        if features.get("flexible_loops"):
            nova_format["flexible_loops"] = features["flexible_loops"][:10]  # Limit size
        
        if features.get("interface_likelihood"):
            # Convert to ranges for efficiency
            positions = sorted(set(features["interface_likelihood"]))
            if positions:
                nova_format["interface_regions"] = [positions[0], positions[-1]]
        
        return nova_format


def main():
    """Test the universal annotator"""
    annotator = UniversalProteinAnnotator()
    
    # Test with a few proteins
    test_proteins = [
        ("TP53", "P04637"),
        ("BRCA1", "P38398"),
        ("CFTR", "P13569")
    ]
    
    results = {}
    for gene, uniprot_id in test_proteins:
        print(f"\n=== Annotating {gene} ===")
        result = annotator.annotate_protein(gene, uniprot_id)
        results[gene] = result
        
        if "error" not in result:
            print(f"✅ Success: {len(result)} features extracted")
        else:
            print(f"❌ Failed: {result['error']}")
    
    # Save results
    with open("universal_annotations.json", "w") as f:
        json.dump({"proteins": results}, f, indent=2)
    
    print(f"\n🚀 Universal annotations saved to universal_annotations.json")


if __name__ == "__main__":
    main()
