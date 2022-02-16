"""
Module holding functions for post poration analysis of volume objects.
"""


from pathlib import Path

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