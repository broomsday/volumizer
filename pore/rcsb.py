"""
Functions for using the RCSB.
"""

from pathlib import Path

from pore.paths import DATA_DIR
from pore.constants import PDB_ID_LENGTH


def get_cluster_file() -> Path:
    """
    """
    cluster_file = DATA_DIR / "rcsb_cluster" / "bc-90.out"
    assert cluster_file.is_file()

    return cluster_file


def parse_cluster_file(lines: list[str]) -> set[str]:
    """
    Take the lines from an RCSB cluster file and return a list of all the PDB IDs in the file.
    """
    return {pdb[:PDB_ID_LENGTH] for pdb in lines}


def build_pdb_list() -> set[str]:
    """
    """

    with open(get_cluster_file(), mode="r", encoding="utf-8") as cluster_file:
        pdbs = parse_cluster_file(cluster_file.readlines())

    return pdbs