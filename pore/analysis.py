"""
Module holding functions for post poration analysis of volume objects.
"""


from pathlib import Path
from typing import Union

import pandas as pd
from tqdm import tqdm

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
        for pdb_id in tqdm(pdb_ids, desc="Getting dataframe locations")
    ]
    missing_annotations = annotation_paths.count(None)
    annotation_paths = [annotation_path for annotation_path in annotation_paths if annotation_path is not None]

    return annotation_paths, missing_annotations


def get_pdb_annotations(annotation_paths: list[Path]) -> dict[str, pd.DataFrame]:
    """
    Create a dictionary of annotation dataframes keyed by the PDB ID and resolution.
    """
    return {annotation_path.stem: pd.read_json(annotation_path) for annotation_path in tqdm(annotation_paths, desc="Compiling dataframes")}


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
        for pdb, annotation in tqdm(pdb_annotations.items(), desc="Filtering PDBs")
        if annotation_satisfies_metrics(annotation, metrics, accepted_types)
    }


def is_stoichiometry_factorable(stoichiometry: dict[int, int]) -> bool:
    """
    If all chains counts are factors of the largest chain count, return True.
    False otherwise
    """
    chain_counts = sorted(list(stoichiometry.values()))
    initial_chain_counts = len(chain_counts)
    if 1 in chain_counts:
        return False

    factorable_counts = [chain_counts.pop()]
    while len(chain_counts) > 0:
        query_count = chain_counts.pop()
        for factorable_count in factorable_counts:
            if factorable_count % query_count == 0:
                factorable_counts.append(query_count)
                break

    if initial_chain_counts == len(factorable_counts):
        return True
    return False


def pdb_satisfies_stoichiometry(stoichiometry: dict[int, int], metric_cutoffs: dict[str, Union[int, bool]]) -> bool:
    """
    If stoichiometry values fall within ranges in `metric_cutoffs` return True.
    False otherwise.
    """
    unique_chains = len(stoichiometry)
    chain_lengths = list(stoichiometry.values())

    # check repeat chains
    if min(chain_lengths) < metric_cutoffs["min_chain_repeats"]:
        return False
    elif (metric_cutoffs["max_chain_repeats"] is not None) and (
        max(chain_lengths) > metric_cutoffs["max_chain_repeats"]
    ):
        return False

    # check unique chains
    if unique_chains < metric_cutoffs["min_unique_chains"]:
        return False
    elif (metric_cutoffs["max_unique_chains"] is not None) and (unique_chains > metric_cutoffs["max_unique_chains"]):
        return False

    # check factorable stoichiometry
    if metric_cutoffs["stoichiometry_factorable"]:
        if not is_stoichiometry_factorable(stoichiometry):
            return False

    return True


def pdb_satisfies_secondary_structure(secondary_structure: dict[str, float], metric_cutoffs: dict[str, float]) -> bool:
    """
    If secondary structure fractions fall within ranges in `metric_cutoffs` return True.
    False otherwise.
    """
    if (
        is_metric_in_range(secondary_structure, "helix", metric_cutoffs)
        and is_metric_in_range(secondary_structure, "strand", metric_cutoffs)
        and is_metric_in_range(secondary_structure, "coil", metric_cutoffs)
    ):
        return True

    return False
