#!/usr/bin/env python3
"""
Generate per-family coefficient JSONs from learning/ data using the unified trainer's
feature extraction, then the calibrator to produce interpretable coefficients.

Usage:
  python -m DNModeling.utils.generate_all_family_coefficients \
      [--families tumor_suppressor ion_channel ...] [--min-n 30]

Outputs:
  DNModeling/cascade/resources/family_models/{FAMILY}_coefficients.json

Notes:
- FAMILY is the cascade classification name (UPPERCASE with underscores),
  e.g., TUMOR_SUPPRESSOR, ION_CHANNEL, COLLAGEN_FIBRILLAR.
- We map learning/ directories to one or more FAMILY names.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd

from DNModeling.utils.unified_family_ml_trainer import UnifiedFamilyMLTrainer
from DNModeling.utils.family_coeff_calibrator import calibrate_and_save

RESOURCES_DIR = Path("DNModeling/cascade/resources/family_models")
LEARNING_DIR = Path("DNModeling/learning")

# Map learning directory -> list of cascade FAMILY names to emit
DIR_TO_FAMILIES: Dict[str, List[str]] = {
    "tumor_suppressor": ["TUMOR_SUPPRESSOR"],
    "oncogene": ["ONCOGENE"],
    "rtk_mapk": ["RTK_MAPK"],
    "ion_channel": ["ION_CHANNEL"],
    "collagen_fibrillar": ["COLLAGEN_FIBRILLAR"],
    "collagen_network": ["COLLAGEN_NETWORK"],
    "collagen_anchoring": ["COLLAGEN_ANCHORING"],
    "collagen_facit": ["COLLAGEN_FACIT"],
    "motor_protein": ["MOTOR_PROTEIN"],
    "transporter": ["TRANSPORTER"],
    "metabolic_enzyme": ["METABOLIC_ENZYME"],
    "signaling_regulator": ["SIGNALING_REGULATOR"],
    "scaffold_adaptor": ["SCAFFOLD_ADAPTOR"],
    # Shared learning bin; emit to both families
    "elastin_fibrillin": ["ELASTIN", "FIBRILLIN"],
    # Broad category; emit to STRUCTURAL as a conservative default
    "cytoskeleton": ["STRUCTURAL"],
    "muscular_dystrophy": ["MUSCULAR_DYSTROPHY"],
    "general": ["GENERAL"],
}


def generate_for_dir(trainer: UnifiedFamilyMLTrainer, dir_name: str, min_n: int) -> List[Path]:
    """Process one learning/ directory, calibrate and save for mapped families."""
    out_paths: List[Path] = []
    families = DIR_TO_FAMILIES.get(dir_name, [dir_name.upper()])

    # Extract features DataFrame from the trainer (no model training required)
    df = trainer.process_family_data(dir_name)
    if df.empty:
        print(f"⏭️  Skipping {dir_name} (no data)")
        return out_paths

    # For calibrator, ensure df has a 'family' column set to the cascade FAMILY name
    for fam in families:
        df_copy = df.copy()
        df_copy["family"] = fam
        try:
            path = calibrate_and_save(df_copy, fam, resources_root=RESOURCES_DIR, min_n=min_n)
            if path:
                out_paths.append(path)
                print(f"✅ Saved coefficients: {path}")
        except Exception as e:
            print(f"❌ Failed to calibrate/save for {dir_name} -> {fam}: {e}")

    return out_paths


def main():
    parser = argparse.ArgumentParser(description="Generate per-family coefficient JSONs")
    parser.add_argument("--families", nargs="*", help="Limit to specific learning/ directories (e.g., tumor_suppressor ion_channel)")
    parser.add_argument("--min-n", type=int, default=30, help="Minimum samples per AA bucket (default: 30)")
    args = parser.parse_args()

    trainer = UnifiedFamilyMLTrainer()
    # Ensure we point at the repo's learning/ directory
    trainer.learning_dir = Path("DNModeling/learning")
    print(f"📁 Using learning dir: {trainer.learning_dir}")

    RESOURCES_DIR.mkdir(parents=True, exist_ok=True)

    dirs: List[str]
    if args.families:
        dirs = args.families
    else:
        # auto-discover learning subdirectories
        dirs = [p.name for p in LEARNING_DIR.iterdir() if p.is_dir()]
        # prioritize some likely-impactful families first
        priority = [
            "tumor_suppressor", "oncogene", "ion_channel",
            "collagen_fibrillar", "motor_protein", "transporter",
            "metabolic_enzyme", "signaling_regulator"
        ]
        # Keep priority order then the rest
        others = [d for d in dirs if d not in priority]
        dirs = priority + others

    print("🚀 Generating per-family coefficients…")
    saved: List[Path] = []
    for d in dirs:
        print(f"\n🧬 Processing {d} …")
        saved += generate_for_dir(trainer, d, min_n=args.min_n)

    print("\n🎉 Done.")
    print(f"💾 Wrote {len(saved)} coefficient files to {RESOURCES_DIR}")


if __name__ == "__main__":
    main()

