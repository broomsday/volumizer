"""
Define project paths.
"""

from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]

C_CODE_DIR = ROOT_DIR / "source"

DATA_DIR = ROOT_DIR / "data"
DOWNLOADED_PDB_DIR = DATA_DIR / "downloaded_pdbs"
PREPARED_PDB_DIR = DATA_DIR / "prepared_pdbs"
ANNOTATED_PDB_DIR = DATA_DIR / "annotated_pdbs"
ANNOTATED_DF_DIR = DATA_DIR / "annotated_dfs"
PDB_FILTERING_METRIC_DIR = DATA_DIR / "pdb_filter_metrics"

RCSB_CLUSTER_FILE = DATA_DIR / "rcsb_cluster" / "bc-90.out"
RCSB_CCD_FILE = DATA_DIR / "rcsb_components" / "components.cif"
PROTEIN_COMPONENTS_FILE = DATA_DIR / "rcsb_components" / "protein_components.txt"
