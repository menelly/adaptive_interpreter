#!/usr/bin/env python3
"""
🔥 HOTSPOT DATABASE - Known Pathogenic and Activating Regions

Manages hotspot region data for genes with known pathogenic clusters.
Extracted from cascade_analyzer.py for modularity.

Built by Ace (2025-10-02) during the Great Refactoring
"""

from typing import Dict, List, Optional


class HotspotDatabase:
    """🔥 HOTSPOT INTELLIGENCE: Known pathogenic and activating regions"""

    def __init__(self):
        # Known hotspot regions from literature and our analysis
        self.hotspots = {
            # Collagen genes - glycine-rich regions are poison hotspots
            'COL1A1': [
                {'start': 359, 'end': 362, 'type': 'dominant_cluster', 'mechanism': 'collagen_poison', 'confidence': 0.9},
                {'start': 988, 'end': 991, 'type': 'dominant_cluster', 'mechanism': 'collagen_poison', 'confidence': 0.8},
            ],
            'COL1A2': [
                {'start': 359, 'end': 362, 'type': 'dominant_cluster', 'mechanism': 'collagen_poison', 'confidence': 0.9},
            ],
            # Ion channels - pore and gate regions
            'SCN1A': [
                {'start': 1616, 'end': 1620, 'type': 'activating_hotspot', 'mechanism': 'channel_disruption', 'confidence': 0.8},
                {'start': 373, 'end': 375, 'type': 'dominant_cluster', 'mechanism': 'channel_poison', 'confidence': 0.7},
            ],
            'KCNQ2': [
                {'start': 265, 'end': 267, 'type': 'activating_hotspot', 'mechanism': 'channel_disruption', 'confidence': 0.9},
            ],
            # Tumor suppressors - DNA binding domains
            'TP53': [
                {'start': 270, 'end': 280, 'type': 'dominant_cluster', 'mechanism': 'dna_binding_loss', 'confidence': 0.95},
            ],
            # Muscle genes - critical structural regions
            'FKRP': [
                {'start': 338, 'end': 340, 'type': 'dominant_cluster', 'mechanism': 'structural_disruption', 'confidence': 0.8},
            ],
        }

    def get_hotspots(self, gene: str) -> List[Dict]:
        """Get all hotspots for a gene"""
        return self.hotspots.get(gene, [])

    def check_variant_in_hotspot(self, gene: str, position: int) -> Optional[Dict]:
        """Check if variant position overlaps any hotspot"""
        hotspots = self.get_hotspots(gene)
        for hotspot in hotspots:
            if hotspot['start'] <= position <= hotspot['end']:
                return hotspot
        return None

