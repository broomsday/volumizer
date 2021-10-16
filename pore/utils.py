"""
Various utility functions
"""


from pathlib import Path

import gzip

from pore.paths import PDB_DIR, CLEANED_PDB_DIR, PROCESSED_PDB_DIR


def get_raw_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return PDB_DIR / f"{pdb_id}.pdb1.gz"


def get_clean_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return CLEANED_PDB_DIR / f"{pdb_id}.pdb"


def get_processed_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return PROCESSED_PDB_DIR / f"{pdb_id}.pdb"


def is_pdb_downloaded(pdb_id: str) -> bool:
    """
    Does the downloaded PDB archive exist.
    """
    return get_raw_pdb_path(pdb_id).is_file()


def is_pdb_cleaned(pdb_id: str) -> bool:
    """
    Has the PDB been fully processed.
    """
    return get_clean_pdb_path(pdb_id).is_file()


def is_pdb_processed(pdb_id: str) -> bool:
    """
    Has the PDB been fully processed.
    """
    return get_processed_pdb_path(pdb_id).is_file()


def decompress_pdb(pdb_id: str) -> None:
    """
    Decompress the gzipped PDB file in the PDB_DIR and save
    the decompressed file into the PROCESSED_PDB_DIR
    """
    with gzip.open(PDB_DIR / f"{pdb_id}.pdb1.gz", mode="rb") as fi:
        pdb_data = fi.read()
    with open(CLEANED_PDB_DIR / f"{pdb_id}.pdb", mode="wb") as fo:
        fo.write(pdb_data)


