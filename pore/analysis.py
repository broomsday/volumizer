"""
Module holding functions for post poration analysis of volume objects.
"""


from pathlib import Path
from typing import Union

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


def compile_accepted_types(metrics: dict[str, Union[bool, float]]) -> set[str]:
    """
    convert metrics to a list of accepted volume types.
    """
    accepted_types = set()
    if metrics["pores"]:
        accepted_types.add("pore")
    if metrics["pockets"]:
        accepted_types.add("pocket")
    if metrics["cavities"]:
        accepted_types.add("cavity")

    return accepted_types


def is_metric_in_range(
    metrics: Union[pd.Series, dict], metric_name: str, metric_cutoffs: dict[str, Union[bool, float]]
) -> bool:
    """
    If the given metric is within the range of min->max inclusive, return True.
    """
    metric_min = metric_cutoffs[f"min_{metric_name}"]
    metric_max = metric_cutoffs[f"max_{metric_name}"]

    if metric_max is None:
        return metrics[metric_name] >= metric_min

    return (metrics[metric_name] >= metric_min) and (metrics[metric_name] <= metric_max)


def annotation_satisfies_metrics(
    annotation: pd.DataFrame, metrics: dict[str, Union[bool, float]], accepted_types: set[str]
) -> bool:
    """
    Ensure at least one volume object in annotation satisfies all metrics.
    """
    for _, row in annotation.iterrows():
        if (
            str(row["type"]) in accepted_types
            and is_metric_in_range(row, "volume", metrics)
            and is_metric_in_range(row, "x", metrics)
            and is_metric_in_range(row, "y", metrics)
            and is_metric_in_range(row, "z", metrics)
        ):
            return True

    return False


def pdb_satisfies_metrics(pdb_metrics: dict[str, int], metric_cutoffs: dict[str, int]) -> bool:
    """
    If all metrics in `pdb_metrics` fall within ranges in `metric_cutoffs` return True.
    False otherwise.
    """
    if (
        is_metric_in_range(pdb_metrics, "atoms", metric_cutoffs)
        and is_metric_in_range(pdb_metrics, "residues", metric_cutoffs)
        and is_metric_in_range(pdb_metrics, "chains", metric_cutoffs)
    ):
        return True

    return False


def select_annotations_by_metrics(
    pdb_annotations: dict[str, pd.DataFrame], metrics: dict[str, Union[bool, float]]
) -> dict[str, pd.DataFrame]:
    """
    Subset the input dictionary into only those where the value annotation has metrics matching those in `metrics`.
    """
    accepted_types = compile_accepted_types(metrics)
    return {
        pdb: annotation
        for pdb, annotation in pdb_annotations.items()
        if annotation_satisfies_metrics(annotation, metrics, accepted_types)
    }
