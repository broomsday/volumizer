"""
Various utility functions
"""


from pathlib import Path
from typing import Optional

import gzip
import tarfile
import numpy as np
import pandas as pd
import json

from volumizer import rcsb, constants, paths
from volumizer.types import Annotation, VoxelGroup
from volumizer.paths import C_CODE_DIR


VOXEL_SIZE = constants.VOXEL_SIZE
VOXEL_VOLUME = VOXEL_SIZE ** 3
KEEP_MODELS = True  # by default keep all models in case they are other biological assembly units
KEEP_NON_PROTEIN = False  # by default only keep protein residues
KEEP_HYDROGENS = False  # by default remove hydrogens


def get_downloaded_pdb_path(pdb_id: str) -> Optional[Path]:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    path = paths.DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb"
    if path.is_file():
        return path

    return None


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
    downloaded_path = get_downloaded_pdb_path(pdb_id)
    if downloaded_path is None:
        return False

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
    paths.PDB_FILTERING_METRIC_DIR.mkdir(parents=True, exist_ok=True)


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
    hubs_df = pd.DataFrame.from_dict(
        [
            {
                "id": i,
                "type": "hub",
                "volume": annotation.hub_volumes[i],
                "x": annotation.hub_dimensions[i][0],
                "y": annotation.hub_dimensions[i][1],
                "z": annotation.hub_dimensions[i][2],
            }
            for i in annotation.hub_volumes.keys()
        ]
    )

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

    return pd.concat([hubs_df, pores_df, cavities_df, pockets_df], ignore_index=True)


def save_annotation_dataframe(name: str, annotation_df: pd.DataFrame):
    """
    Save the annotation dataframe.
    """
    annotation_df.to_json(paths.ANNOTATED_DF_DIR / f"{name}.{str(VOXEL_SIZE)}.json")


def have_annotation(file_stem: str) -> bool:
    """
    If we have already completed the annotation of this file, return True.
    False otherwise.
    """
    pdb_path = paths.ANNOTATED_PDB_DIR / f"{file_stem}.{VOXEL_SIZE}.pdb"
    df_path = paths.ANNOTATED_DF_DIR / f"{file_stem}.{VOXEL_SIZE}.json"
    if pdb_path.is_file() and df_path.is_file():
        return True

    return False


def load_annotation_df(file_stem: str) -> pd.DataFrame:
    """
    Return the annotation dataframe associated with this file-stem.
    """
    return pd.read_json(paths.ANNOTATED_DF_DIR / f"{file_stem}.{VOXEL_SIZE}.json")


def using_performant() -> bool:
    """
    Return True if the performant .so libraries exist.
    """
    if (C_CODE_DIR / "voxel.so").is_file() and (C_CODE_DIR / "fib_sphere.so").is_file():
        return True

    return False


def have_pdb_size_metrics_on_file(pdb_id: str) -> bool:
    """
    Check to see if the PDB size metrics have been recorded previously.
    """
    return (paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_size.json").is_file()


def save_pdb_size_metrics(pdb_id: str, metrics: dict[str, int]) -> None:
    """
    Save the number of atoms, residues, and chains in a PDB assembly.
    """
    with open(paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_size.json", mode="w", encoding="utf-8") as out_file:
        json.dump(metrics, out_file)


def load_pdb_size_metrics(pdb_id: str) -> Optional[dict[str, int]]:
    """
    Load the number of atoms, residues, and chains for a PDB
    """
    pdb_size_metric_path = paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_size.json"
    if not pdb_size_metric_path.is_file():
        return None

    with open(pdb_size_metric_path, mode="r", encoding="utf-8") as in_file:
        return json.load(in_file)


def have_stoichiometry_on_file(pdb_id: str) -> bool:
    """
    Check to see if the stoichiometry of the PDB assembly has been recorded previously.
    """
    return (paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_stoichiometry.json").is_file()


def save_stoichiometry(pdb_id: str, metrics: dict[int, int]) -> None:
    """
    Save the stoichiometry of the PDB assembly.
    """
    with open(paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_stoichiometry.json", mode="w", encoding="utf-8") as out_file:
        json.dump(metrics, out_file)


def load_stoichiometry(pdb_id: str) -> Optional[dict[int, int]]:
    """
    Load the stoichiometry of the PDB assembly.
    """
    pdb_size_metric_path = paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_stoichiometry.json"
    if not pdb_size_metric_path.is_file():
        return None

    with open(pdb_size_metric_path, mode="r", encoding="utf-8") as in_file:
        return json.load(in_file)


def have_secondary_structure_on_file(pdb_id: str) -> bool:
    """
    Check to see if the secondary structure has been recorded previously.
    """
    return (paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_secondary_structure.json").is_file()


def save_secondary_structure(pdb_id: str, metrics: dict[int, int]) -> None:
    """
    Save the secondary structure fractions.
    """
    with open(
        paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_secondary_structure.json", mode="w", encoding="utf-8"
    ) as out_file:
        json.dump(metrics, out_file)


def load_secondary_structure(pdb_id: str) -> Optional[dict[int, int]]:
    """
    Load the secondary structure fractions.
    """
    pdb_size_metric_path = paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_secondary_structure.json"
    if not pdb_size_metric_path.is_file():
        return None

    with open(pdb_size_metric_path, mode="r", encoding="utf-8") as in_file:
        return json.load(in_file)
