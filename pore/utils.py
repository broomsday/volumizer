"""
Various utility functions
"""


import gzip

from pore.paths import PDB_DIR, CLEANED_PDB_DIR


def is_pdb_downloaded(pdb_id: str) -> bool:
    """
    Does the downloaded PDB archive exist.
    """
    return (PDB_DIR / f"{pdb_id}.pdb1.gz").is_file()


def is_pdb_cleaned(pdb_id: str) -> bool:
    """
    Has the PDB been fully processed.
    """
    return (CLEANED_PDB_DIR / f"{pdb_id}.pdb").is_file()


def decompress_pdb(pdb_id: str) -> None:
    """
    Decompress the gzipped PDB file in the PDB_DIR and save
    the decompressed file into the PROCESSED_PDB_DIR
    """
    with gzip.open(PDB_DIR / f"{pdb_id}.pdb1.gz", mode="rb") as fi:
        pdb_data = fi.read()
    with open(CLEANED_PDB_DIR / f"{pdb_id}.pdb", mode="wb") as fo:
        fo.write(pdb_data)
