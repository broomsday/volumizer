"""
Functions for using the RCSB.
"""

from pathlib import Path
from urllib import request
from urllib.error import HTTPError, URLError
from time import sleep
import re

from pymongo import database
from tqdm import tqdm

from pore.paths import RCSB_CLUSTER_FILE, PDB_DIR, RCSB_CCD_FILE
from pore.constants import PDB_ID_LENGTH, RCSB_CLUSTER_URL, RCSB_BIOUNIT_URL, RCSB_STRUCTURE_URL, RCSB_CCD_URL
from pore.types import ComponentData
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


def get_cluster_file() -> None:
    """
    Download the RCSB cluster file if it doesn't exist locally.
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


def download_biological_assemblies(db: database.Database) -> None:
    """
    If we have not already downloaded a PDB, do so now.
    """
    PDB_DIR.mkdir(exist_ok=True, parents=True)
    for pdb in tqdm(list(db.pdbs.find()), "Downloading PDBs"):
        if not pdb["downloaded"]:
            if utils.is_pdb_downloaded(pdb["pdb_id"]):
                mongo.update_pdb_downloaded(db, pdb, True)
            else:
                mongo.update_pdb_downloaded(db, pdb, download_biological_assembly(pdb["pdb_id"]))


def component_file_exists(component_file: Path) -> bool:
    """
    Check if the Chemical Component Dictionary file exists locally.
    """
    return component_file.is_file()


def download_component_file() -> None:
    """
    Download the Chemical Component Dictionary file.
    """
    request.urlretrieve(RCSB_CCD_URL, RCSB_CCD_FILE)


def get_component_file() -> None:
    """
    Download the Chemical Component Dictionary file if it doesn't exist locally.
    """
    RCSB_CCD_FILE.parent.mkdir(exist_ok=True, parents=True)
    download_attempts = 0
    while (not component_file_exists(RCSB_CCD_FILE)) and (download_attempts < 10):
        download_component_file()
        download_attempts += 1

    assert component_file_exists(RCSB_CCD_FILE)


def format_component_id_line(component_id_line: str) -> str:
    """
    Pull out just the 3-letter ID from a CCD component ID line.
    """
    return component_id_line.strip().split()[-1]


def format_component_type_line(component_type_line: str) -> str:
    """
    Pull out just the type without quotes from a CCD component type line.
    """
    return " ".join(component_type_line.strip().split()[1:]).strip("\"")


def parse_component_file(lines: list[str]) -> list[ComponentData]:
    """
    Read lines from a CCD .cif file and return a list of all component ids and a lsit of all component types.
    """
    id_search = re.compile("_chem_comp.id*")
    id_lines = list(filter(id_search.match, lines))

    type_search = re.compile("_chem_comp.type*")
    type_lines = list(filter(type_search.match, lines))

    component_ids = [format_component_id_line(id_line) for id_line in id_lines]
    component_types = [format_component_type_line(type_line) for type_line in type_lines]

    return [ComponentData(component_id=component_id, component_type=component_type) for component_id, component_type in zip(component_ids, component_types)]


def is_component_protein(component: ComponentData) -> bool:
    """
    Return whether the supplied component is a protein residue.

    This is based on having L-LINKING or D-LINKING properties
    """
    return ("PEPTIDE LINKING" in component.component_type) or ("petide linking" in component.component_type)


def process_all_components(db: database.Database, components: list[ComponentData]) -> None:
    """
    Assign to the component database if the type of each component is peptide linking or not.
    """
    processed_component_ids = [component["component_id"] for component in db.components.find({"processed": True})]
    if len(processed_component_ids) != len(components):
        for component in tqdm(components, "Processing Components"):
            if not component.component_id in processed_component_ids:
                mongo.update_component_is_protein(db, component, is_component_protein(component))


def build_component_set(component_file: Path) -> list[ComponentData]:
    """
    Get a set of all the component IDs we want to download and process.
    """
    assert component_file_exists(component_file)
    with open(component_file, mode="r", encoding="utf-8") as fi:
        component_data = parse_component_file(fi.readlines())

    return component_data
