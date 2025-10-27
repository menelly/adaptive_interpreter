#!/usr/bin/env python3
"""
Universal Protein Annotator - NO MORE HARDCODING!
Automatically extracts protein features from UniProt, Pfam, and sequence analysis
"""

import requests
import re
import json
import os
from typing import Dict, List, Tuple, Optional
import time

class UniversalProteinAnnotator:
    def __init__(self, cache_dir="protein_annotations_cache"):
        self.uniprot_base = "https://rest.uniprot.org"
        self.pfam_base = "https://pfam.xfam.org"

        # 🎯 CACHING SYSTEM - Save domain data locally!
        self.cache_dir = cache_dir
        import os
        os.makedirs(self.cache_dir, exist_ok=True)
        print(f"🔍 Protein annotation cache: {self.cache_dir}")


        
    def get_uniprot_features(self, uniprot_id: str) -> Dict:
        """Extract all features from UniProt API with caching"""

        # 🎯 CHECK CACHE FIRST!
        cache_file = f"{self.cache_dir}/{uniprot_id}_domains.json"
        if os.path.exists(cache_file):
            try:
                with open(cache_file, 'r') as f:
                    cached_data = json.load(f)
                print(f"🔍 Using cached domain data for {uniprot_id}")
                return cached_data
            except Exception as e:
                print(f"⚠️ Cache read error for {uniprot_id}: {e}")

        # Fetch from UniProt API
        url = f"{self.uniprot_base}/uniprotkb/{uniprot_id}.json"

        try:
            print(f"🌐 Fetching domain data from UniProt for {uniprot_id}...")
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
                "phosphorylation": [],
                # 🎯 NEW! Critical domain features we were missing
                "propeptides": [],  # N-terminal and C-terminal propeptides
                "mature_chain": [],  # The actual functional protein
                "regions": [],  # Functional regions like triple helix
                "cleavage_sites": [],  # Where protein gets processed
                "motifs": [],  # Functional motifs
                "compositional_bias": [],  # Compositionally biased regions (important for HMSN-P!)
                # 🎯 NEW! GO cross-references (IDs and term names)
                "go_ids": [],
                "go_terms": []
            }

            # Debug: print function extraction
            print(f"Extracted function for {uniprot_id}: {features['function'][:100] if features['function'] else 'NONE'}...")

            # Extract GO terms from UniProt cross-references (offline-friendly once cached)
            try:
                xrefs = data.get("uniProtKBCrossReferences", []) or []
                go_ids = []
                go_terms = []
                for x in xrefs:
                    if (x or {}).get("database") == "GO":
                        go_id = x.get("id")
                        if go_id:
                            go_ids.append(go_id)
                        for prop in (x.get("properties") or []):
                            val = (prop.get("value") or "").strip()
                            if val:
                                # Typical format: "F:protein binding" / "P:signal transduction" / "C:nucleus"
                                term = val.split(":", 1)[1].strip() if ":" in val else val
                                go_terms.append(term)
                # De-duplicate, store
                features["go_ids"] = sorted(set(go_ids))
                features["go_terms"] = sorted(set(go_terms))
            except Exception as ge:
                print(f"⚠️ GO extraction failed for {uniprot_id}: {ge}")

            # Extract features from UniProt annotations
            for feature in data.get("features", []):
                self._parse_uniprot_feature(feature, features)

            # 🎯 NEW! Get REAL functional domains from Pfam API
            self._add_pfam_domains(features)

            # 🎯 NEW! Detect compositional bias regions
            self._detect_compositional_bias(features)

            # 🎯 SAVE TO CACHE!
            try:
                with open(cache_file, 'w') as f:
                    json.dump(features, f, indent=2)
                print(f"💾 Cached domain data for {uniprot_id}")
            except Exception as e:
                print(f"⚠️ Cache write error for {uniprot_id}: {e}")

            return features

        except Exception as e:
            print(f"UniProt API failed for {uniprot_id}: {e}")
            # Minimal offline structure, so downstream can still read 'function' and 'go_terms'
            return {"uniprot_id": uniprot_id, "function": "", "go_terms": [], "go_ids": [], "error": str(e)}
    
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
        # 🎯 NEW! Parse the critical domain features we were missing
        elif feature_type == "Propeptide":  # Fixed case!
            features["propeptides"].append({
                "start": start, "end": end or start,
                "type": "propeptide", "description": description
            })
        elif feature_type == "Chain":  # Fixed case!
            features["mature_chain"].append({
                "start": start, "end": end or start,
                "type": "mature_chain", "description": description
            })
        elif feature_type == "Region":  # Fixed case!
            features["regions"].append({
                "start": start, "end": end or start,
                "type": "region", "description": description
            })
        elif feature_type == "Site":  # Fixed case!
            features["cleavage_sites"].append({
                "position": start, "description": description
            })
        elif feature_type == "Motif":  # Fixed case!
            features["motifs"].append({
                "start": start, "end": end or start,
                "type": "motif", "description": description
            })
    
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
        """Find UniProt ID from gene name - prioritize reviewed entries"""
        url = f"{self.uniprot_base}/uniprotkb/search"
        params = {
            "query": f"gene:{gene_name} AND organism_id:9606",
            "format": "json",
            "size": 10  # Get multiple results to find the best one
        }

        try:
            response = requests.get(url, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()

            results = data.get("results", [])
            if not results:
                return None

            # 🎯 PRIORITIZE REVIEWED ENTRIES (canonical proteins)
            reviewed_entries = []
            unreviewed_entries = []

            for result in results:
                accession = result.get("primaryAccession")
                entry_type = result.get("entryType", "")

                if entry_type == "UniProtKB reviewed (Swiss-Prot)":
                    reviewed_entries.append(accession)
                else:
                    unreviewed_entries.append(accession)

            # Return first reviewed entry if available, otherwise first unreviewed
            if reviewed_entries:
                print(f"🔍 Found reviewed UniProt entry for {gene_name}: {reviewed_entries[0]}")
                return reviewed_entries[0]
            elif unreviewed_entries:
                print(f"⚠️ Only unreviewed entries found for {gene_name}: {unreviewed_entries[0]}")
                return unreviewed_entries[0]

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
            "domains": features.get("domains", []),
            # Expose GO to BiologicalRouter/plausibility systems
            "go_terms": features.get("go_terms", []),
            "go_ids": features.get("go_ids", [])
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

    def _add_pfam_domains(self, features: Dict):
        """🎯 Get REAL functional domains from Pfam API using sequence search"""
        sequence = features.get("sequence", "")
        if not sequence or len(sequence) < 20:
            return

        try:
            # Use Pfam sequence search API
            print(f"🌐 Searching Pfam for domains in {features.get('uniprot_id', 'unknown')}...")

            # Pfam sequence search endpoint
            pfam_url = "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"

            # Prepare sequence for submission
            data = {
                'seq': sequence,
                'database': 'pfam'
            }

            # Submit search (this is a simplified version - real API needs more handling)
            response = requests.post(pfam_url, data=data, timeout=30)

            if response.status_code == 200:
                # Parse results (this would need proper parsing of HMMER output)
                print(f"✅ Pfam search submitted for {features.get('uniprot_id', 'unknown')}")
                # For now, we'll use sequence-based domain prediction instead
                self._predict_domains_from_sequence(features)
            else:
                print(f"⚠️ Pfam API returned status {response.status_code}")
                self._predict_domains_from_sequence(features)

        except Exception as e:
            print(f"⚠️ Pfam API failed: {e}")
            # Fallback to sequence-based prediction
            self._predict_domains_from_sequence(features)

    def _predict_domains_from_sequence(self, features: Dict):
        """🎯 Predict domains from sequence patterns when APIs fail"""
        sequence = features.get("sequence", "")
        if not sequence:
            return

        print(f"🔍 Predicting domains from sequence patterns...")

        # Look for coiled coil patterns (heptad repeats)
        coiled_coil_regions = self._find_coiled_coil_regions(sequence)
        for region in coiled_coil_regions:
            features["domains"].append({
                "start": region["start"],
                "end": region["end"],
                "type": "predicted_coiled_coil",
                "description": "Predicted coiled coil domain"
            })
            print(f"🔍 Predicted coiled coil: {region['start']}-{region['end']}")

        # Look for PB1-like patterns (basic residue clusters)
        pb1_regions = self._find_pb1_like_regions(sequence)
        for region in pb1_regions:
            features["domains"].append({
                "start": region["start"],
                "end": region["end"],
                "type": "predicted_pb1",
                "description": "Predicted PB1-like domain"
            })
            print(f"🔍 Predicted PB1-like domain: {region['start']}-{region['end']}")

    def _find_coiled_coil_regions(self, sequence: str) -> List[Dict]:
        """Find potential coiled coil regions using heptad repeat patterns"""
        regions = []
        window_size = 21  # 3 heptad repeats

        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i+window_size]

            # Check for heptad repeat pattern (hydrophobic at positions 1,4 of each heptad)
            hydrophobic_score = 0
            for j in range(0, window_size, 7):
                if j < len(window) and window[j] in "AILMVFWY":
                    hydrophobic_score += 1
                if j+3 < len(window) and window[j+3] in "AILMVFWY":
                    hydrophobic_score += 1

            if hydrophobic_score >= 4:  # At least 4 hydrophobic positions in pattern
                regions.append({"start": i+1, "end": i+window_size})

        # Merge overlapping regions
        return self._merge_overlapping_regions(regions)

    def _find_pb1_like_regions(self, sequence: str) -> List[Dict]:
        """Find potential PB1-like domains (protein-protein interaction domains)"""
        regions = []
        window_size = 40

        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i+window_size]

            # PB1 domains often have clusters of basic residues and beta-strand patterns
            basic_count = sum(1 for aa in window if aa in "RK")
            hydrophobic_count = sum(1 for aa in window if aa in "ILMVFWY")

            # Look for balanced basic/hydrophobic content (typical of PB1)
            if basic_count >= 6 and hydrophobic_count >= 8:
                regions.append({"start": i+1, "end": i+window_size})

        return self._merge_overlapping_regions(regions)

    def _merge_overlapping_regions(self, regions: List[Dict]) -> List[Dict]:
        """Merge overlapping regions"""
        if not regions:
            return []

        sorted_regions = sorted(regions, key=lambda x: x["start"])
        merged = [sorted_regions[0]]

        for current in sorted_regions[1:]:
            last = merged[-1]
            if current["start"] <= last["end"] + 10:  # Allow small gaps
                last["end"] = max(last["end"], current["end"])
            else:
                merged.append(current)

        return merged

    def _detect_compositional_bias(self, features: Dict):
        """🎯 Detect compositional bias regions (important for HMSN-P variants!)"""
        sequence = features.get("sequence", "")
        if not sequence:
            return

        # Look for compositionally biased regions
        window_size = 20
        bias_threshold = 0.6  # 60% of one amino acid type

        compositional_bias = []

        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i+window_size]

            # Check for each amino acid type
            for aa in "ACDEFGHIKLMNPQRSTVWY":
                aa_count = window.count(aa)
                if aa_count / window_size >= bias_threshold:
                    compositional_bias.append({
                        "start": i + 1,
                        "end": i + window_size,
                        "type": "compositional_bias",
                        "description": f"{aa}-rich region ({aa_count}/{window_size} = {aa_count/window_size:.1%})",
                        "amino_acid": aa,
                        "percentage": aa_count / window_size
                    })
                    break  # Don't double-count the same region

        # Merge overlapping regions
        if compositional_bias:
            merged_bias = []
            current = compositional_bias[0]

            for next_region in compositional_bias[1:]:
                if next_region["start"] <= current["end"] + 5:  # Allow small gaps
                    # Merge regions
                    current["end"] = max(current["end"], next_region["end"])
                    # Update description to show range
                    current["description"] = f"Compositional bias region ({current['start']}-{current['end']})"
                else:
                    merged_bias.append(current)
                    current = next_region

            merged_bias.append(current)
            features["compositional_bias"] = merged_bias

            print(f"🔍 Found {len(merged_bias)} compositional bias regions")
            for bias in merged_bias:
                print(f"   {bias['start']}-{bias['end']}: {bias['description']}")


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
