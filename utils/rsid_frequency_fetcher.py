#!/usr/bin/env python3
"""
🔑 rsID FREQUENCY FETCHER - Nova's Brilliant Solution!
Skip gnomAD coordinate hell - use rsIDs as universal passports!

Built by Ace (2025) based on Nova's genius insight
"rsIDs are the universal passport for variants"
"""

import requests
import time
import json
import csv
import re
import os
from typing import Dict, List, Optional, Any
from pathlib import Path

class RSIDFrequencyFetcher:
    """Fetch variant frequencies using rsIDs via Ensembl REST API"""
    
    def __init__(self, cache_file: str = "rsid_frequency_cache.json", delay: float = 0.2):
        self.cache_file = cache_file
        self.delay = delay  # Be kind to Ensembl servers
        self.cache = self.load_cache()
        self.session = requests.Session()
        self.session.headers.update({"Content-Type": "application/json"})
        
        print("🔑 rsID Frequency Fetcher initialized!")
        print(f"📊 Loaded {len(self.cache)} cached frequencies")
    
    def load_cache(self) -> Dict[str, Any]:
        """Load cached frequency data"""
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"⚠️ Error loading cache: {e}")
        return {}
    
    def save_cache(self):
        """Save frequency cache"""
        try:
            with open(self.cache_file, 'w') as f:
                json.dump(self.cache, f, indent=2)
            print(f"💾 Saved {len(self.cache)} frequencies to cache")
        except Exception as e:
            print(f"⚠️ Error saving cache: {e}")
    
    def fetch_frequency(self, rsid: str) -> Optional[Dict[str, Any]]:
        """
        Fetch frequency data for a single rsID
        
        Args:
            rsid: rsID (e.g., "rs2043293738")
            
        Returns:
            Dict with frequency data or None if failed
        """
        
        # Check cache first
        if rsid in self.cache:
            return self.cache[rsid]
        
        # Clean rsID
        if not rsid.startswith('rs'):
            rsid = f"rs{rsid}"
        
        url = f"https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"
        
        try:
            print(f"🔍 Fetching {rsid}...")
            response = self.session.get(url)
            
            if response.status_code != 200:
                print(f"⚠️ Failed for {rsid}: HTTP {response.status_code}")
                result = {"rsid": rsid, "max_af": 0.0, "error": f"HTTP {response.status_code}"}
                self.cache[rsid] = result
                return result
            
            data = response.json()
            
            # Extract population frequencies
            pop_genotypes = data.get("population_genotypes", [])
            
            # Find maximum allele frequency across all populations
            max_af = 0.0
            cohort_freqs = {}
            
            for pop in pop_genotypes:
                freq = pop.get("frequency", 0.0)
                population = pop.get("population", "unknown")
                
                if freq > max_af:
                    max_af = freq
                
                cohort_freqs[population] = freq
            
            # Extract gnomAD frequencies specifically if available
            gnomad_af = 0.0
            for pop in pop_genotypes:
                if "gnomad" in pop.get("population", "").lower():
                    gnomad_af = max(gnomad_af, pop.get("frequency", 0.0))
            
            result = {
                "rsid": rsid,
                "max_af": max_af,
                "gnomad_af": gnomad_af,
                "cohort_freqs": cohort_freqs,
                "total_populations": len(pop_genotypes)
            }
            
            # Cache result
            self.cache[rsid] = result
            
            print(f"✅ {rsid}: max_af={max_af:.6f}, gnomad_af={gnomad_af:.6f}")
            
            # Be kind to Ensembl
            time.sleep(self.delay)
            
            return result
            
        except Exception as e:
            print(f"❌ Error fetching {rsid}: {e}")
            result = {"rsid": rsid, "max_af": 0.0, "error": str(e)}
            self.cache[rsid] = result
            return result
    
    def extract_rsids_from_csv(self, csv_file: str, rsid_column: str = None) -> List[str]:
        """
        Extract rsIDs from CSV file
        
        Args:
            csv_file: Path to CSV file
            rsid_column: Column name containing rsIDs (auto-detect if None)
            
        Returns:
            List of unique rsIDs
        """
        
        rsids = set()
        
        try:
            with open(csv_file, 'r') as f:
                # Try to detect delimiter
                sample = f.read(1024)
                f.seek(0)
                
                delimiter = '\t' if '\t' in sample else ','
                reader = csv.DictReader(f, delimiter=delimiter)
                
                # Auto-detect rsID column if not specified
                if rsid_column is None:
                    headers = reader.fieldnames
                    rsid_candidates = [h for h in headers if 'rsid' in h.lower() or 'rs' in h.lower()]
                    if rsid_candidates:
                        rsid_column = rsid_candidates[0]
                        print(f"🔍 Auto-detected rsID column: {rsid_column}")
                    else:
                        print("⚠️ No rsID column found, searching all columns...")
                
                for row in reader:
                    if rsid_column and rsid_column in row:
                        rsid = row[rsid_column].strip()
                        if rsid and rsid.startswith('rs'):
                            rsids.add(rsid)
                    else:
                        # Search all columns for rsIDs
                        for value in row.values():
                            if value and isinstance(value, str):
                                # Look for rsID patterns
                                rs_matches = re.findall(r'rs\d+', value)
                                rsids.update(rs_matches)
        
        except Exception as e:
            print(f"❌ Error reading CSV: {e}")
            return []
        
        rsid_list = sorted(list(rsids))
        print(f"🔍 Found {len(rsid_list)} unique rsIDs")
        return rsid_list
    
    def filter_variants_by_frequency(self, csv_file: str, output_file: str, 
                                   patho_af_threshold: float = 0.01,
                                   benign_af_threshold: float = 0.001):
        """
        Filter variants based on frequency thresholds
        
        Args:
            csv_file: Input CSV file
            output_file: Output filtered CSV file
            patho_af_threshold: Max AF for pathogenic variants (default: 1%)
            benign_af_threshold: Min AF for benign variants (default: 0.1%)
        """
        
        print(f"🔍 Filtering variants from {csv_file}")
        print(f"📊 Thresholds: Patho AF < {patho_af_threshold}, Benign AF > {benign_af_threshold}")
        
        # Extract rsIDs
        rsids = self.extract_rsids_from_csv(csv_file)
        
        if not rsids:
            print("❌ No rsIDs found in CSV file")
            return
        
        # Fetch frequencies for all rsIDs
        print(f"🚀 Fetching frequencies for {len(rsids)} rsIDs...")
        
        for rsid in rsids:
            self.fetch_frequency(rsid)
        
        # Save cache after fetching
        self.save_cache()
        
        # Filter CSV based on frequencies
        self.apply_frequency_filter(csv_file, output_file, patho_af_threshold, benign_af_threshold)
    
    def apply_frequency_filter(self, input_file: str, output_file: str,
                             patho_af_threshold: float, benign_af_threshold: float):
        """Apply frequency filtering to CSV file"""
        
        filtered_rows = []
        total_rows = 0
        
        try:
            with open(input_file, 'r') as f:
                # Detect delimiter
                sample = f.read(1024)
                f.seek(0)
                delimiter = '\t' if '\t' in sample else ','
                
                reader = csv.DictReader(f, delimiter=delimiter)
                headers = reader.fieldnames
                
                for row in reader:
                    total_rows += 1
                    
                    # Find rsID in row
                    rsid = None
                    for value in row.values():
                        if value and isinstance(value, str):
                            rs_match = re.search(r'rs\d+', value)
                            if rs_match:
                                rsid = rs_match.group()
                                break
                    
                    if not rsid or rsid not in self.cache:
                        # Keep variants without rsIDs (might be very rare)
                        filtered_rows.append(row)
                        continue
                    
                    freq_data = self.cache[rsid]
                    max_af = freq_data.get("max_af", 0.0)
                    
                    # Determine if variant should be kept based on classification
                    classification = ""
                    for value in row.values():
                        if value and isinstance(value, str):
                            if any(x in value.upper() for x in ['PATHOGENIC', 'BENIGN']):
                                classification = value.upper()
                                break
                    
                    keep_variant = True
                    
                    if 'PATHOGENIC' in classification and 'BENIGN' not in classification:
                        # Pathogenic variants should be rare
                        if max_af > patho_af_threshold:
                            print(f"⚠️ Suspicious pathogenic {rsid}: AF={max_af:.6f} > {patho_af_threshold}")
                            keep_variant = False
                    
                    elif 'BENIGN' in classification and 'PATHOGENIC' not in classification:
                        # Benign variants should have some frequency
                        if max_af < benign_af_threshold:
                            print(f"⚠️ Suspicious benign {rsid}: AF={max_af:.6f} < {benign_af_threshold}")
                            keep_variant = False
                    
                    if keep_variant:
                        # Add frequency data to row
                        row['rsid'] = rsid
                        row['max_af'] = max_af
                        row['gnomad_af'] = freq_data.get("gnomad_af", 0.0)
                        filtered_rows.append(row)
            
            # Write filtered results
            if filtered_rows:
                with open(output_file, 'w', newline='') as f:
                    # Add frequency columns to headers
                    new_headers = list(headers) + ['rsid', 'max_af', 'gnomad_af']
                    writer = csv.DictWriter(f, fieldnames=new_headers, delimiter=delimiter)
                    writer.writeheader()
                    writer.writerows(filtered_rows)
                
                print(f"✅ Filtered {total_rows} → {len(filtered_rows)} variants")
                print(f"📊 Results saved to {output_file}")
            else:
                print("❌ No variants passed filtering")
        
        except Exception as e:
            print(f"❌ Error filtering variants: {e}")

def main():
    """Test the rsID frequency fetcher"""
    fetcher = RSIDFrequencyFetcher()
    
    # Test with a few rsIDs
    test_rsids = ["rs2043293738", "rs2044669685", "rs2141266590"]
    
    print("🧪 Testing rsID frequency fetching...")
    for rsid in test_rsids:
        result = fetcher.fetch_frequency(rsid)
        if result:
            print(f"  {rsid}: {result}")
    
    fetcher.save_cache()

if __name__ == "__main__":
    main()
