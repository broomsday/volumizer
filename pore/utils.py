"""
Various utility functions
"""


from pathlib import Path

import gzip

from pore.paths import DOWNLOADED_PDB_DIR, PREPARED_PDB_DIR, ANNOTATED_PDB_DIR, PROTEIN_COMPONENTS_FILE
from pore import rcsb


def get_downloaded_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb"


def get_prepared_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return PREPARED_PDB_DIR / f"{pdb_id}.pdb"


def get_annotated_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return ANNOTATED_PDB_DIR / f"{pdb_id}.pdb"


def is_pdb_downloaded(pdb_id: str) -> bool:
    """
    Does the downloaded PDB archive exist.
    """
    return get_downloaded_pdb_path(pdb_id).is_file()


def is_pdb_prepared(pdb_id: str) -> bool:
    """
    Has the PDB been fully processed.
    """
    return get_prepared_pdb_path(pdb_id).is_file()


def is_pdb_annotated(pdb_id: str) -> bool:
    """
    Has the PDB been fully processed.
    """
    return get_annotated_pdb_path(pdb_id).is_file()


def decompress_pdb(pdb_id: str) -> None:
    """
    Decompress the gzipped PDB file in the PDB_DIR and save
    the decompressed file into the PROCESSED_PDB_DIR
    """
    zipped_path = DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb1.gz"
    unzipped_path = DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb"

    with gzip.open(zipped_path, mode="rb") as fi:
        pdb_data = fi.read()
    with open(unzipped_path, mode="wb") as fo:
        fo.write(pdb_data)

    zipped_path.unlink()


def setup_dirs():
    """
    Make sure the base data directories exist.
    """
    DOWNLOADED_PDB_DIR.mkdir(parents=True, exist_ok=True)
    PREPARED_PDB_DIR.mkdir(parents=True, exist_ok=True)
    ANNOTATED_PDB_DIR.mkdir(parents=True, exist_ok=True)


def load_protein_components() -> set[str]:
    """
    Load the processed protein components file.
    """
    with open(PROTEIN_COMPONENTS_FILE, mode="r", encoding="utf-8") as fi:
        return set([component.rstrip("\n") for component in fi.readlines()])


def ensure_protein_components_file() -> None:
    """
    Make sure the processed components file exists.
    If not already present, generate it.
    """
    if not PROTEIN_COMPONENTS_FILE.is_file():
        components = rcsb.get_components()
        protein_components = rcsb.build_protein_component_set(components)
        with open(PROTEIN_COMPONENTS_FILE, mode="w", encoding="utf-8") as fo:
            fo.write("\n".join([component for component in protein_components]))
