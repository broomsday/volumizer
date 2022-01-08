"""
Various utility functions
"""


from pathlib import Path

import gzip
import numpy as np

from pore.paths import DOWNLOADED_PDB_DIR, PREPARED_PDB_DIR, ANNOTATED_PDB_DIR, PROTEIN_COMPONENTS_FILE
from pore import rcsb
from pore import constants
from pore.types import Annotation, VoxelGroup


VOXEL_SIZE = constants.VOXEL_SIZE


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


def set_resolution(resolution: float) -> None:
    """
    Set the value of the VOXEL_SIZE global constant.
    """
    global VOXEL_SIZE
    VOXEL_SIZE = resolution


def get_volume_summary(voxel_group_dict: dict[int, VoxelGroup], summary_type: str = "total") -> float:
    """
    Compute a summary value for the volume
    """
    volumes = [voxel_group.volume for voxel_group in voxel_group_dict.values() if voxel_group.volume is not None]
    if volumes:
        if summary_type == "total":
            return sum(volumes)
        elif summary_type == "max":
            return max(volumes)
        elif summary_type == "mean":
            return np.mean(volumes)
        else:
            print("Unsupported volume summary type")

    return 0.0


def print_annotation(annotation: Annotation) -> None:
    """
    Print the annotation to the terminal in a readable manner.
    """
    print("")
    print(f"Number of pores: {annotation.num_pores}")
    print(f"Number of cavities: {annotation.num_cavities}")
    print(f"Number of pockets: {annotation.num_pockets}")
    print("")
    print(f"Largest pore volume: {annotation.largest_pore_volume}")
    print(f"Largest cavity volume: {annotation.largest_cavity_volume}")
    print(f"Largest pocket volume: {annotation.largest_pocket_volume}")
    print("")
    print(f"Total pore volume: {annotation.total_pore_volume}")
    print(f"Total cavity volume: {annotation.total_cavity_volume}")
    print(f"Total pocket volume: {annotation.total_pocket_volume}")


def save_annotated_pdb(pdb_name: str, annotated_lines: list[str]) -> None:
    """
    Save a PDB formatted coordinate file of the voxels.
    Individual atoms/voxels are labelled according to type
    """
    # get the original PDB lines so that our annotation can be appended
    with open(PREPARED_PDB_DIR / f"{pdb_name}.pdb", mode="r", encoding="utf-8") as input_pdb_file:
        input_pdb_lines = input_pdb_file.readlines()
    input_pdb_lines = [line.rstrip("\n") for line in input_pdb_lines]

    # remove the terminal END line
    for line in input_pdb_lines[::-1]:
        if ("END" in line) or (line.strip() == ""):
            input_pdb_lines.pop()
    
    # add the annotated lines and save
    with open(ANNOTATED_PDB_DIR / f"{pdb_name}.pdb", mode="w", encoding="utf-8") as annotated_pdb_file:
        annotated_pdb_file.write("\n".join([*input_pdb_lines, "END", *annotated_lines]))


def sort_voxelgroups_by_volume(voxelgroups: dict[int, VoxelGroup]) -> dict[int, VoxelGroup]:
    """
    Take a dictionary with indices as the keys and VoxelGroups as the values.
    Reassign the keys such that the VoxelGroups are sorted by volume.
    """
    return {i: voxelgroup for i, voxelgroup in enumerate(sorted(voxelgroups.values(), key=lambda group: group.volume, reverse=True))}
