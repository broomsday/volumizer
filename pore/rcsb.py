"""
Functions for using the RCSB.
"""

from pathlib import Path
from urllib import request
from urllib.error import HTTPError, URLError
from time import sleep
import re

from pore.paths import RCSB_CLUSTER_FILE, DOWNLOADED_PDB_DIR, RCSB_CCD_FILE
from pore.constants import (
    PDB_ID_LENGTH,
    RCSB_CLUSTER_URL,
    RCSB_BIOUNIT_URL,
    RCSB_BUNDLE_URL,
    RCSB_STRUCTURE_URL,
    RCSB_CCD_URL,
)
from pore.types import ComponentData


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


def download_biological_assembly(pdb_id: str, retries: int = 10) -> bool:
    """
    Check to see if the biological assembly is available at the RCSB.
    If so, download it.

    If not, check to see if the PDB-like bundle is available.
    If not, just down the standard PDB.
    """
    biounit_url = RCSB_BIOUNIT_URL + f"{pdb_id[1:3].lower()}/{pdb_id.lower()}.pdb1.gz"
    bundle_url = RCSB_BUNDLE_URL + f"{pdb_id[1:3].lower()}/{pdb_id.lower()}/{pdb_id.lower()}-pdb-bundle.tar.gz"
    structure_url = RCSB_STRUCTURE_URL + f"{pdb_id[1:3].lower()}/pdb{pdb_id.lower()}.ent.gz"

    for _ in range(retries):
        try:
            request.urlretrieve(biounit_url, DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb1.gz")
            return True
        except HTTPError:
            try:
                request.urlretrieve(bundle_url, DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb1.tar.gz")
                return True
            except HTTPError:
                try:
                    request.urlretrieve(structure_url, DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb1.gz")
                    return True
                except HTTPError:
                    return False
        except URLError:
            sleep(1)

    return False


def component_file_exists() -> bool:
    """
    Check if the Chemical Component Dictionary file exists locally.
    """
    return RCSB_CCD_FILE.is_file()


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
    while (not component_file_exists()) and (download_attempts < 10):
        download_component_file()
        download_attempts += 1

    assert component_file_exists()


def format_component_id_line(component_id_line: str) -> str:
    """
    Pull out just the 3-letter ID from a CCD component ID line.
    """
    return component_id_line.strip().split()[-1]


def format_component_type_line(component_type_line: str) -> str:
    """
    Pull out just the type without quotes from a CCD component type line.
    """
    return " ".join(component_type_line.strip().split()[1:]).strip('"')


def parse_component_file(lines: list[str]) -> list[ComponentData]:
    """
    Read lines from a CCD .cif file and return a list of all component ids and types.
    """
    id_search = re.compile("_chem_comp.id*")
    id_lines = list(filter(id_search.match, lines))

    type_search = re.compile("_chem_comp.type*")
    type_lines = list(filter(type_search.match, lines))

    component_ids = [format_component_id_line(id_line) for id_line in id_lines]
    component_types = [format_component_type_line(type_line) for type_line in type_lines]

    return [
        ComponentData(component_id=component_id, component_type=component_type)
        for component_id, component_type in zip(component_ids, component_types)
    ]


def is_component_protein(component: ComponentData) -> bool:
    """
    Return whether the supplied component is a protein residue.

    This is based on having L-LINKING or D-LINKING properties
    """
    return ("PEPTIDE LINKING" in component.component_type) or ("petide linking" in component.component_type)


def build_protein_component_set(components: list[ComponentData]) -> set[str]:
    """
    Build a set of the component id (name) of all components that are protein
    """
    return {component.component_id for component in components if is_component_protein(component)}


def get_components() -> list[ComponentData]:
    """
    Get all components and their types.
    """
    if not component_file_exists():
        get_component_file()
    with open(RCSB_CCD_FILE, mode="r", encoding="utf-8") as fi:
        component_data = parse_component_file(fi.readlines())

    return component_data
