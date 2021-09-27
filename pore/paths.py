"""
Define project paths.
"""

from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT_DIR / "data"
PDB_DIR = DATA_DIR / "pdbs"
PROCESSED_PDB_DIR = DATA_DIR / "processed_pdbs"
ALIGNED_PORE_DIR = DATA_DIR / "aligned_pores"
RCSB_CLUSTER_FILE = DATA_DIR / "rcsb_cluster" / "bc-90.out"
