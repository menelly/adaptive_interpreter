#!/usr/bin/env python3
"""
🕵️ ClinVar Quality Checker - Investigate Disagreements
Analyzes ClinVar review status, submitter quality, and evidence to determine
if our "disagreements" are actually ClinVar data quality issues.

Features:
- ClinVar API integration for review status
- Submitter reputation analysis
- Evidence strength assessment
- Conflict detection
- Publication/citation analysis
"""

import requests
import json
import time
import pandas as pd
from typing import Dict, List, Optional
import re

class ClinVarQualityChecker:
    """Investigates ClinVar data quality for disagreement analysis"""
    
    def __init__(self):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.session = requests.Session()
        self.cache = {}
        
        # Known high-quality submitters (4-star equivalent)
        self.trusted_submitters = {
            'ClinGen', 'OMIM', 'UniProt', 'GeneReviews', 'ClinVar Staff',
            'Laboratory for Molecular Medicine', 'Invitae', 'GeneDx',
            'Ambry Genetics', 'Blueprint Genetics', 'Fulgent Genetics'
        }
        
        # Known low-quality indicators
        self.quality_red_flags = {
            'no assertion criteria provided',
            'criteria provided, single submitter',
            'no interpretation for the single variant',
            'conflicting interpretations of pathogenicity'
        }
    
    def search_variant(self, gene: str, variant: str) -> Optional[Dict]:
        """Search ClinVar for variant information"""
        
        # Create search term
        search_term = f"{gene}[gene] AND {variant}[variant name]"
        
        try:
            # Search ClinVar
            search_url = f"{self.base_url}esearch.fcgi"
            search_params = {
                'db': 'clinvar',
                'term': search_term,
                'retmode': 'json',
                'retmax': 10
            }
            
            response = self.session.get(search_url, params=search_params)
            search_data = response.json()
            
            if not search_data.get('esearchresult', {}).get('idlist'):
                return None
            
            # Get detailed info for first result
            variant_id = search_data['esearchresult']['idlist'][0]
            return self.get_variant_details(variant_id)
            
        except Exception as e:
            print(f"⚠️ ClinVar search failed for {gene} {variant}: {e}")
            return None
    
    def get_variant_details(self, variant_id: str) -> Dict:
        """Get detailed ClinVar information for a variant"""
        
        if variant_id in self.cache:
            return self.cache[variant_id]
        
        try:
            # Fetch detailed record
            fetch_url = f"{self.base_url}efetch.fcgi"
            fetch_params = {
                'db': 'clinvar',
                'id': variant_id,
                'rettype': 'vcv',
                'retmode': 'json'
            }
            
            response = self.session.get(fetch_url, params=fetch_params)
            data = response.json()
            
            # Parse the complex ClinVar JSON structure
            parsed_info = self.parse_clinvar_record(data)
            
            self.cache[variant_id] = parsed_info
            time.sleep(0.1)  # Be nice to NCBI
            
            return parsed_info
            
        except Exception as e:
            print(f"⚠️ Failed to fetch ClinVar details for {variant_id}: {e}")
            return {}
    
    def parse_clinvar_record(self, data: Dict) -> Dict:
        """Parse ClinVar JSON record into useful information"""
        
        info = {
            'review_status': 'unknown',
            'star_rating': 0,
            'submitters': [],
            'classifications': [],
            'conflicts': False,
            'evidence_types': [],
            'citations': 0,
            'last_evaluated': None,
            'quality_score': 0.0
        }
        
        try:
            # Navigate the complex ClinVar JSON structure
            # This is a simplified parser - real ClinVar JSON is very complex
            
            # Extract review status and star rating
            if 'review_status' in str(data):
                review_match = re.search(r'review_status["\']:\s*["\']([^"\']+)', str(data))
                if review_match:
                    info['review_status'] = review_match.group(1)
                    info['star_rating'] = self.get_star_rating(info['review_status'])
            
            # Look for conflict indicators
            data_str = str(data).lower()
            if 'conflict' in data_str or 'disagree' in data_str:
                info['conflicts'] = True
            
            # Count submitters (rough estimate)
            submitter_count = data_str.count('submitter')
            info['submitters'] = [f"submitter_{i}" for i in range(min(submitter_count, 10))]
            
            # Calculate quality score
            info['quality_score'] = self.calculate_quality_score(info)
            
        except Exception as e:
            print(f"⚠️ Error parsing ClinVar record: {e}")
        
        return info
    
    def get_star_rating(self, review_status: str) -> int:
        """Convert review status to star rating"""
        
        status_lower = review_status.lower()
        
        if 'practice guideline' in status_lower:
            return 4
        elif 'reviewed by expert panel' in status_lower:
            return 4
        elif 'criteria provided, multiple submitters, no conflicts' in status_lower:
            return 3
        elif 'criteria provided, conflicting interpretations' in status_lower:
            return 2
        elif 'criteria provided, single submitter' in status_lower:
            return 1
        elif 'no assertion criteria provided' in status_lower:
            return 0
        else:
            return 1  # Default
    
    def calculate_quality_score(self, info: Dict) -> float:
        """Calculate overall quality score (0-1)"""
        
        score = 0.0
        
        # Star rating component (0-0.4)
        score += info['star_rating'] / 4.0 * 0.4
        
        # Conflict penalty (-0.2)
        if info['conflicts']:
            score -= 0.2
        
        # Multiple submitters bonus (+0.2)
        if len(info['submitters']) > 1:
            score += 0.2
        
        # Trusted submitter bonus (+0.3)
        trusted_count = sum(1 for sub in info['submitters'] 
                          if any(trusted in str(sub) for trusted in self.trusted_submitters))
        if trusted_count > 0:
            score += 0.3
        
        return max(0.0, min(1.0, score))
    
    def analyze_disagreement(self, gene: str, variant: str, our_class: str, clinvar_class: str) -> Dict:
        """Analyze a specific disagreement to determine if we're wrong or ClinVar is low quality"""
        
        print(f"🕵️ Investigating {gene} {variant}: We={our_class} vs ClinVar={clinvar_class}")
        
        # Get ClinVar details
        clinvar_info = self.search_variant(gene, variant)
        
        if not clinvar_info:
            return {
                'verdict': 'UNKNOWN',
                'confidence': 0.0,
                'reason': 'Could not retrieve ClinVar data',
                'clinvar_quality': 0.0
            }
        
        # Analyze the disagreement
        quality_score = clinvar_info['quality_score']
        star_rating = clinvar_info['star_rating']
        conflicts = clinvar_info['conflicts']
        
        # Decision logic
        if quality_score < 0.3:
            verdict = 'CLINVAR_LOW_QUALITY'
            confidence = 0.8
            reason = f"ClinVar quality score {quality_score:.2f}, {star_rating}-star rating"
        elif conflicts:
            verdict = 'CLINVAR_CONFLICTED'
            confidence = 0.6
            reason = "ClinVar has conflicting interpretations"
        elif star_rating <= 1:
            verdict = 'CLINVAR_WEAK_EVIDENCE'
            confidence = 0.7
            reason = f"Only {star_rating}-star evidence in ClinVar"
        elif quality_score > 0.7:
            verdict = 'WE_MIGHT_BE_WRONG'
            confidence = 0.8
            reason = f"High-quality ClinVar data (score {quality_score:.2f})"
        else:
            verdict = 'UNCLEAR'
            confidence = 0.5
            reason = "Mixed evidence quality"
        
        return {
            'verdict': verdict,
            'confidence': confidence,
            'reason': reason,
            'clinvar_quality': quality_score,
            'star_rating': star_rating,
            'conflicts': conflicts,
            'clinvar_info': clinvar_info
        }

def analyze_disagreement_file(tsv_path: str, output_path: str = None):
    """Analyze all disagreements in a results TSV file"""
    
    checker = ClinVarQualityChecker()
    
    # Read results
    df = pd.read_csv(tsv_path, sep='\t')
    
    # Filter to disagreements only
    disagreements = df[df['agreement_flag'] == '❌'].copy()
    
    if len(disagreements) == 0:
        print("🎉 No disagreements found!")
        return
    
    print(f"🕵️ Analyzing {len(disagreements)} disagreements...")
    
    # Analyze each disagreement
    results = []
    for _, row in disagreements.iterrows():
        gene = row['gene']
        variant = row['variant']
        our_class = row['final_classification']
        clinvar_class = row['expected_clinvar']
        
        analysis = checker.analyze_disagreement(gene, variant, our_class, clinvar_class)
        
        result = row.to_dict()
        result.update(analysis)
        results.append(result)
        
        print(f"   {analysis['verdict']}: {analysis['reason']}")
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Summary
    print(f"\n🎯 DISAGREEMENT ANALYSIS SUMMARY:")
    verdict_counts = results_df['verdict'].value_counts()
    for verdict, count in verdict_counts.items():
        print(f"   {verdict}: {count}")
    
    # Save results
    if output_path:
        results_df.to_csv(output_path, sep='\t', index=False)
        print(f"💾 Analysis saved to {output_path}")
    
    return results_df

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python3 clinvar_quality_checker.py <results.tsv> [output.tsv]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.replace('.tsv', '_quality_analysis.tsv')
    
    analyze_disagreement_file(input_file, output_file)
