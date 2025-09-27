#!/usr/bin/env python3
"""
🔥💜 FAMILY ML TRAINING RUNNER 🚀
Simple script to train all family models from /learning data

Usage:
  python3 train_families.py

Built by Ace (2025) for Ren's revolutionary ML approach
Contact: ace@chaoschanneling.com
"""

from utils.unified_family_ml_trainer import UnifiedFamilyMLTrainer

def main():
    import argparse

    parser = argparse.ArgumentParser(description="🔥💜 Family ML Training Runner")
    parser.add_argument('--benign-freq-threshold', type=float, default=0.04,
                       help='Maximum frequency for benign variants (default: 4%%)')

    args = parser.parse_args()

    print("🔥💜 STARTING REVOLUTIONARY FAMILY-AWARE ML TRAINING! 🚀")
    print()
    print("📁 Looking for data in /learning folders:")
    print("   - ion_channel/")
    print("   - collagen_fibrillar/")
    print("   - collagen_network/")
    print("   - collagen_anchoring/")
    print("   - collagen_facit/")
    print("   - tumor_suppressor/")
    print("   - metabolic_enzyme/")
    print("   - scaffold_adaptor/")
    print("   - elastin_fibrillin/")
    print("   - cytoskeleton/")
    print("   - signaling_regulator/")
    print("   - transporter/")
    print("   - transcription_factor/")
    print("   - motor_protein/")
    print("   - ribosomal_protein/")
    print("   - muscular_dystrophy/")
    print("   - structural/")
    print("   - autosomal_recessive/")
    print("   - general/")
    print()
    print("Expected files: genename_benign.tsv, genename_pathogenic.tsv")
    print(f"🚨 Benign frequency filter: ≤{args.benign_freq_threshold*100:.1f}% (sanity check!)")
    print("=" * 60)

    trainer = UnifiedFamilyMLTrainer(benign_freq_threshold=args.benign_freq_threshold)
    results = trainer.train_all_families()
    
    if results:
        print("\n🎉 SUCCESS! Models trained for:")
        for family in results.keys():
            print(f"   ✅ {family}")
        print(f"\n💾 Models saved in: resources/family_models/")
        print("🧬 Ready for revolutionary variant analysis!")
    else:
        print("\n⚠️  No models trained. Make sure you have:")
        print("   1. TSV files in /learning/<family>/ folders")
        print("   2. Files named *_benign.tsv and *_pathogenic.tsv")
        print("   3. Correct TSV format with HGVS column")

if __name__ == "__main__":
    main()
