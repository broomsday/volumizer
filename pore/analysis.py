"""
Module holding functions for post poration analysis of volume objects.
"""


from pathlib import Path

import pandas as pd

from pore.paths import ANNOTATED_DF_DIR
from pore.utils import VOXEL_SIZE


def get_annotations_by_id(pdb_ids: list[str]) -> tuple[list[Path], int]:
    """
    Given a list of PDB IDs, return the corresponding list of annotated
    dataframe file pathes.
    """
    annotation_paths = [
        (ANNOTATED_DF_DIR / f"{pdb_id}.{VOXEL_SIZE}.json")
        if (ANNOTATED_DF_DIR / f"{pdb_id}.{VOXEL_SIZE}.json").is_file()
        else None
        for pdb_id in pdb_ids
    ]
    missing_annotations = annotation_paths.count(None)
    annotation_paths = [annotation_path for annotation_path in annotation_paths if annotation_path is not None]

    return annotation_paths, missing_annotations


def get_pdb_annotations(annotation_paths: list[Path]) -> dict[str, pd.DataFrame]:
    """
    Create a dictionary of annotation dataframes keyed by the PDB ID and resolution.
    """
    return {annotation_path.stem: pd.read_json(annotation_path) for annotation_path in annotation_paths}


def annotation_has_types(annotation: pd.DataFrame, pores: bool, pockets: bool, cavities: bool) -> bool:
    """
    If the given annotation contains any of pores/pockets/cavities that are True, return True.
    """
    if annotation.empty:
        return False

    accepted_types = []
    if pores:
        accepted_types.append("pore")
    if pockets:
        accepted_types.append("pocket")
    if cavities:
        accepted_types.append("cavity")

    for accepted_type in accepted_types:
        if accepted_type in list(annotation["type"]):
            return True

    return False


def select_annotations_by_type(
    pdb_annotations: dict[str, pd.DataFrame], pores: bool, pockets: bool, cavities: bool
) -> dict[str, pd.DataFrame]:
    """
    Subset the input dictionary into only those where the value annotation has at least one type
    matching at least one of pores/pockets/cavities that are True
    """
    return {
        pdb: annotation
        for pdb, annotation in pdb_annotations.items()
        if annotation_has_types(annotation, pores, pockets, cavities)
    }
