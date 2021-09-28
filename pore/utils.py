"""
Various utility functions
"""


from pore.paths import PDB_DIR, PROCESSED_PDB_DIR


def is_pdb_downloaded(pdb_id: str) -> bool:
    """
    Does the downloaded PDB archive exist.
    """
    return (PDB_DIR / f"{pdb_id}.pdb1.gz").is_file()


def is_pdb_processed(pdb_id: str) -> bool:
    """
    Has the PDB been fully processed.
    """
    return (PROCESSED_PDB_DIR / f"{pdb_id}.pdb").is_file()
