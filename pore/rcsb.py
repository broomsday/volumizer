"""
Functions for using the RCSB.
"""

from pathlib import Path
from urllib import request
from urllib.error import HTTPError, URLError
from time import sleep

from pymongo import database
from tqdm import tqdm

from pore.paths import RCSB_CLUSTER_FILE, PDB_DIR
from pore.constants import PDB_ID_LENGTH, RCSB_CLUSTER_URL, RCSB_BIOUNIT_URL, RCSB_STRUCTURE_URL
from pore import mongo
from pore import utils


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
    RCSB_CLUSTER_FILE.parent.mkdir(exist_ok=True, parents=True)
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


def download_biological_assembly(pdb_id: str, retries: int=10) -> bool:
    """
    Check to see if the biological assembly is available at the RCSB.
    If so, download it.

    If not, download the standard PDB.

    In both cases the file downloaded is compressed.
    """
    biounit_url = RCSB_BIOUNIT_URL + f"{pdb_id[1:3].lower()}/{pdb_id.lower()}.pdb1.gz"
    structure_url = RCSB_STRUCTURE_URL + f"{pdb_id[1:3].lower()}/pdb{pdb_id.lower()}.ent.gz"

    for _ in range(retries):
        try:
            request.urlretrieve(biounit_url, PDB_DIR / f"{pdb_id}.pdb1.gz")
            return True
        except HTTPError:
            try:
                request.urlretrieve(structure_url, PDB_DIR / f"{pdb_id}.pdb1.gz")
                return True
            except HTTPError:
                return False
        except URLError:
            sleep(1)

    return False


def download_biological_assemblies(db: database.Database):
    """
    If we have not already downloaded a PDB, do so now.
    """
    PDB_DIR.mkdir(exist_ok=True, parents=True)
    for pdb in tqdm(list(db.pdbs.find()), "Downloading PDBs"):
        if not pdb["downloaded"]:
            mongo.update_downloaded(db, pdb, download_biological_assembly(pdb["pdb_id"]))
