"""
Functions for using the RCSB.
"""

from pathlib import Path

from urllib import request

from pore.paths import RCSB_CLUSTER_FILE
from pore.constants import PDB_ID_LENGTH, RCSB_CLUSTER_URL


def cluster_file_exists(cluster_file: Path) -> bool:
    """
    Check if the RCSB cluster file exists locally.
    """
    return cluster_file.is_file()


def download_cluster_file() -> None:
    """
    Download the RCSB cluster file.
    """
    request.urlretrieve(RCSB_CLUSTER_URL, RCSB_CLUSTER_FILE)


def get_rcsb_cluster_file() -> None:
    """
    Download the RCSB cluster file if it doesn't exist locally.

    Return the path to the RCSB cluster file.
    """
    download_attempts = 0
    while (not cluster_file_exists(RCSB_CLUSTER_FILE)) and (download_attempts < 10):
        download_cluster_file()
        download_attempts += 1

    assert cluster_file_exists(RCSB_CLUSTER_FILE)
        

def parse_cluster_file(lines: list[str]) -> set[str]:
    """
    Take the lines from an RCSB cluster file and return a list of all the PDB IDs in the file.
    """
    return {pdb[:PDB_ID_LENGTH] for pdb in lines}


def build_pdb_set(cluster_file: Path) -> set[str]:
    """
    Get a set of all the PDB IDs we want to download and process.
    """
    assert cluster_file_exists(cluster_file)
    with open(cluster_file, mode="r", encoding="utf-8") as fi:
        pdbs = parse_cluster_file(fi.readlines())

    return pdbs