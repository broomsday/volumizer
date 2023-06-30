"""
Various utility functions
"""


from typing import Optional

import numpy as np
import pandas as pd

from volumizer import constants, protein_components
from volumizer.types import Annotation, VoxelGroup
from volumizer.paths import C_CODE_DIR


VOXEL_SIZE = constants.VOXEL_SIZE
VOXEL_VOLUME = VOXEL_SIZE ** 3
KEEP_MODELS = True  # by default keep all models in case they are other biological assembly units
KEEP_NON_PROTEIN = False  # by default only keep protein residues
KEEP_HYDROGENS = False  # by default remove hydrogens
UTILS_PROTEIN_COMPONENTS = protein_components.PROTEIN_COMPONENTS


def get_protein_components() -> set[str]:
    """
    Return the protein components definition currently in use.
    """
    return UTILS_PROTEIN_COMPONENTS


def reset_protein_components():
    """
    Reset the protein components definition to the default.
    """
    global UTILS_PROTEIN_COMPONENTS

    UTILS_PROTEIN_COMPONENTS = protein_components.PROTEIN_COMPONENTS

def add_protein_components(additional_components: set[str]) -> None:
    """
    Extend the protein components definition currently in use with additional components
    """
    global UTILS_PROTEIN_COMPONENTS
    
    UTILS_PROTEIN_COMPONENTS = UTILS_PROTEIN_COMPONENTS.union(additional_components)


def remove_protein_components(disallowed_components: set[str]) -> None:
    """
    Reduce the protein components definition currently in use with blacklisted components.
    """
    global UTILS_PROTEIN_COMPONENTS
    
    UTILS_PROTEIN_COMPONENTS = UTILS_PROTEIN_COMPONENTS - disallowed_components


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


def using_performant() -> bool:
    """
    Return True if the performant .so libraries exist.
    """
    if (C_CODE_DIR / "voxel.so").is_file() and (C_CODE_DIR / "fib_sphere.so").is_file():
        return True

    return False
