"""
📜 AdaptiveInterpreter Configuration
Centralized configuration for all paths and settings.
"""

from pathlib import Path

# --- Base Paths ---
# Use a sensible default, assuming data is stored in a predictable location relative to home.
# Users can override this by creating a .env file or setting environment variables.
# For this project, we'll assume a base data directory in the user's home.
# The original path was /mnt/Arcana/, which we are restoring.
BASE_DATA_PATH = Path("/mnt/Arcana")

# --- Specific Data Paths ---
ALPHAODL_STRUCTURES_PATH = BASE_DATA_PATH / "alphafold_human" / "structures"
GNOMAD_DATA_PATH = BASE_DATA_PATH / "gnomad"
CONSERVATION_DATA_PATH = BASE_DATA_PATH / "UCSC"
RESOURCES_PATH = Path(__file__).parent / "resources"

# You can add other configurations here, like API keys, thresholds, etc.
# Example:
# API_KEYS = {
#     "GNOMAD_API": "your_key_here"
# }

print(f"🧬 AdaptiveInterpreter configured with base data path: {BASE_DATA_PATH}")
