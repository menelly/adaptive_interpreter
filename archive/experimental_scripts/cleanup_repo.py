#!/usr/bin/env python3
"""
🧹 Repository Cleanup Script
Organize the DNModeling repo from "teenager's room" to "organized drawers"

Based on Nova's brilliant organization plan (2025)
"""

import os
import shutil
import pathlib
from typing import List

ROOT = pathlib.Path(__file__).parent

def create_folders():
    """Create the new folder structure"""
    folders = [
        "nova_dn",
        "cascade", 
        "utils",
        "archive",
        "docs",
        "tests/unit",
        "tests/results", 
        "tests/data",
        "resources"
    ]
    
    for folder in folders:
        (ROOT / folder).mkdir(parents=True, exist_ok=True)
        print(f"📁 Created: {folder}")

def move_core_analyzers():
    """Move core DN/LOF/GOF analyzers to nova_dn/"""
    core_files = [
        "analyzer.py",
        "lof_analyzer.py", 
        "gof_variant_analyzer.py",
        "amino_acid_props.py",
        "mechanisms.py",
        "lattice_disruption.py",
        "context.py",
        "universal_context.py",
        "dn_mechanism_filter_v2.py",
        "sequence_manager.py"
    ]
    
    for fname in core_files:
        src = ROOT / fname
        if src.exists():
            shutil.move(str(src), ROOT / "nova_dn" / fname)
            print(f"🧬 Moved to nova_dn/: {fname}")

def move_cascade_files():
    """Move cascade orchestration files"""
    cascade_files = [
        "cascade_analyzer.py",
        "cascade_batch_processor.py", 
        "hotspot_database.py",
        "biological_router.py"
    ]
    
    for fname in cascade_files:
        src = ROOT / fname
        if src.exists():
            shutil.move(str(src), ROOT / "cascade" / fname)
            print(f"🌊 Moved to cascade/: {fname}")

def move_utilities():
    """Move utility and ML files"""
    util_files = [
        "alphafold_sequence.py",
        "conservation_database.py",
        "proline_ml_integrator.py", 
        "gly_cys_simple_integrator.py",
        "gnomad_frequency_fetcher.py",
        "mixed_mechanism_resolver.py",
        "conservation_ml_trainer.py",
        "conservation_ml_loader.py"
    ]
    
    for fname in util_files:
        src = ROOT / fname
        if src.exists():
            shutil.move(str(src), ROOT / "utils" / fname)
            print(f"🔧 Moved to utils/: {fname}")

def archive_old_files():
    """Archive unused/legacy files"""
    archive_files = [
        "population_frequency_analyzer.py",
        "dn_mechanism_filter.py",
        "smart_protein_analyzer.py", 
        "tune.py",
        "motifs.py",
        "run_batch.py"
    ]
    
    for fname in archive_files:
        src = ROOT / fname
        if src.exists():
            shutil.move(str(src), ROOT / "archive" / fname)
            print(f"📦 Archived: {fname}")

def move_documentation():
    """Move documentation files"""
    # Move all .md files except README.md
    for md_file in ROOT.glob("*.md"):
        if md_file.name.lower() != "readme.md":
            shutil.move(str(md_file), ROOT / "docs" / md_file.name)
            print(f"📚 Moved to docs/: {md_file.name}")

def move_test_files():
    """Move test and result files"""
    # Move test data files
    test_patterns = ["*test*", "*Test*", "*TEST*"]
    for pattern in test_patterns:
        for f in ROOT.glob(pattern):
            if f.is_file():
                shutil.move(str(f), ROOT / "tests/data" / f.name)
                print(f"🧪 Moved to tests/data/: {f.name}")
    
    # Move results directory if it exists
    results_dir = ROOT / "results"
    if results_dir.exists() and results_dir.is_dir():
        # Move contents to tests/results
        for item in results_dir.iterdir():
            shutil.move(str(item), ROOT / "tests/results" / item.name)
        # Remove empty results directory
        results_dir.rmdir()
        print(f"📊 Moved results/ contents to tests/results/")

def move_resources():
    """Move configuration and resource files"""
    resource_files = [
        "protein_annotations.json",
        "tuned_weights.json",
        "conservation_multipliers.json",
        "comprehensive_gene_cache.json"
    ]
    
    for fname in resource_files:
        src = ROOT / fname
        if src.exists():
            shutil.move(str(src), ROOT / "resources" / fname)
            print(f"📋 Moved to resources/: {fname}")

def create_init_files():
    """Create __init__.py files for Python packages"""
    packages = ["nova_dn", "cascade", "utils"]
    
    for package in packages:
        init_file = ROOT / package / "__init__.py"
        if not init_file.exists():
            init_file.write_text('"""DNModeling package"""\n')
            print(f"🐍 Created __init__.py in {package}/")

def main():
    """Run the complete cleanup"""
    print("🧹 Starting repository cleanup...")
    print("📁 Nova's organization plan in action!")
    print("=" * 50)
    
    create_folders()
    print()
    
    move_core_analyzers()
    print()
    
    move_cascade_files() 
    print()
    
    move_utilities()
    print()
    
    archive_old_files()
    print()
    
    move_documentation()
    print()
    
    move_test_files()
    print()
    
    move_resources()
    print()
    
    create_init_files()
    print()
    
    print("=" * 50)
    print("✨ Repository cleanup complete!")
    print("🎉 From teenager's room to organized drawers!")
    print("💜 Thanks Nova for the brilliant organization plan!")

if __name__ == "__main__":
    main()
