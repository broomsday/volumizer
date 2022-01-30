"""
Various utility functions
"""


from pathlib import Path
from typing import Optional

import gzip
import tarfile
import numpy as np
import pandas as pd

from pore import rcsb, constants, paths
from pore.types import Annotation, VoxelGroup


VOXEL_SIZE = constants.VOXEL_SIZE
VOXEL_VOLUME = VOXEL_SIZE ** 3
KEEP_MODELS = True  # by default keep all models in case they are other biological assembly units
KEEP_NON_PROTEIN = False  # by default only keep protein residues
KEEP_HYDROGENS = False  # by default remove hydrogens


def get_downloaded_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return paths.DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb"


def get_prepared_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return paths.PREPARED_PDB_DIR / f"{pdb_id}.pdb"


def get_annotated_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return paths.ANNOTATED_PDB_DIR / f"{pdb_id}.pdb"


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
    Decompress the gzipped PDB file in the DOWNLOADED_PDB_DIR.
    Save as a text file.
    Delete the original compressed object.
    """
    tar_gz_path = paths.DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb1.tar.gz"
    gz_path = paths.DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb1.gz"
    unzipped_path = paths.DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb"

    if gz_path.is_file():
        with gzip.open(gz_path, mode="rb") as fi:
            pdb_data = fi.read()
        with open(unzipped_path, mode="wb") as fo:
            fo.write(pdb_data)
        gz_path.unlink()
    elif tar_gz_path.is_file():
        with tarfile.open(tar_gz_path, mode="r") as fi:
            name = [name for name in fi.getnames() if "bundle" in name][0]
            member = fi.getmember(name)
            pdb_data = fi.extractfile(member).readlines()
        with open(unzipped_path, mode="wb") as fo:
            fo.writelines(pdb_data)
        tar_gz_path.unlink()


def setup_dirs():
    """
    Make sure the base data directories exist.
    """
    paths.DOWNLOADED_PDB_DIR.mkdir(parents=True, exist_ok=True)
    paths.PREPARED_PDB_DIR.mkdir(parents=True, exist_ok=True)
    paths.ANNOTATED_PDB_DIR.mkdir(parents=True, exist_ok=True)
    paths.ANNOTATED_DF_DIR.mkdir(parents=True, exist_ok=True)


def load_protein_components() -> set[str]:
    """
    Load the processed protein components file.
    """
    with open(paths.PROTEIN_COMPONENTS_FILE, mode="r", encoding="utf-8") as fi:
        return set([component.rstrip("\n") for component in fi.readlines()])


def ensure_protein_components_file() -> None:
    """
    Make sure the processed components file exists.
    If not already present, generate it.
    """
    if not paths.PROTEIN_COMPONENTS_FILE.is_file():
        components = rcsb.get_components()
        protein_components = rcsb.build_protein_component_set(components)
        with open(paths.PROTEIN_COMPONENTS_FILE, mode="w", encoding="utf-8") as fo:
            fo.write("\n".join([component for component in protein_components]))


def set_resolution(resolution: float) -> None:
    """
    Set the value of the VOXEL_SIZE and VOXEL_VOLUME global constants.
    """
    global VOXEL_SIZE
    global VOXEL_VOLUME

    VOXEL_SIZE = resolution
    VOXEL_VOLUME = VOXEL_SIZE ** 3


def set_non_protein(non_protein: bool) -> None:
    """
    Set whether to include non-protein residues during the calculations.
    """
    global KEEP_NON_PROTEIN
    KEEP_NON_PROTEIN = non_protein


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


def print_annotation(annotation: Annotation, annotation_df: pd.DataFrame) -> None:
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
    print("")

    print(annotation_df)


def sort_voxelgroups_by_volume(voxelgroups: dict[int, VoxelGroup]) -> dict[int, VoxelGroup]:
    """
    Take a dictionary with indices as the keys and VoxelGroups as the values.
    Reassign the keys such that the VoxelGroups are sorted by volume.
    """
    return {
        i: voxelgroup
        for i, voxelgroup in enumerate(sorted(voxelgroups.values(), key=lambda group: group.volume, reverse=True))
    }


def filter_voxelgroups_by_volume(
    voxelgroups: dict[int, VoxelGroup], min_volume: Optional[float] = None, min_voxels: Optional[int] = None
) -> dict[int, VoxelGroup]:
    """
    Remove voxel groups based on min volumem and min number of voxels cutoff
    """
    if min_volume is not None:
        voxelgroups = {
            i: voxelgroup for i, voxelgroup in enumerate(voxelgroups.values()) if voxelgroup.volume >= min_volume
        }
    if min_voxels is not None:
        voxelgroups = {
            i: voxelgroup for i, voxelgroup in enumerate(voxelgroups.values()) if voxelgroup.num_voxels >= min_voxels
        }

    return voxelgroups


def make_annotation_dataframe(annotation: Annotation) -> pd.DataFrame:
    """
    Return a dataframe of the annotation values for each volume type.
    """
    pores_df = pd.DataFrame.from_dict(
        [
            {
                "id": i,
                "type": "pore",
                "volume": annotation.pore_volumes[i],
                "x": annotation.pore_dimensions[i][0],
                "y": annotation.pore_dimensions[i][1],
                "z": annotation.pore_dimensions[i][2],
            }
            for i in annotation.pore_volumes.keys()
        ]
    )

    cavities_df = pd.DataFrame.from_dict(
        [
            {
                "id": i,
                "type": "cavity",
                "volume": annotation.cavity_volumes[i],
                "x": annotation.cavity_dimensions[i][0],
                "y": annotation.cavity_dimensions[i][1],
                "z": annotation.cavity_dimensions[i][2],
            }
            for i in annotation.cavity_volumes.keys()
        ]
    )

    pockets_df = pd.DataFrame.from_dict(
        [
            {
                "id": i,
                "type": "pocket",
                "volume": annotation.pocket_volumes[i],
                "x": annotation.pocket_dimensions[i][0],
                "y": annotation.pocket_dimensions[i][1],
                "z": annotation.pocket_dimensions[i][2],
            }
            for i in annotation.pocket_volumes.keys()
        ]
    )

    return pd.concat([pores_df, cavities_df, pockets_df], ignore_index=True)


def save_annotation_dataframe(name: str, annotation_df: pd.DataFrame):
    """
    Save the annotation dataframe.
    """
    annotation_df.to_json(paths.ANNOTATED_DF_DIR / f"{name}.json")
