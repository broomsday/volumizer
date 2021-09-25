"""
Define project paths.
"""

from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT_DIR / "data"
RCSB_CLUSTER_FILE = DATA_DIR / "rcsb_cluster" / "bc-90.out"
