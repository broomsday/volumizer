"""
Entry-level functions to find pores and cavities.
"""


from pathlib import Path
from time import perf_counter
from typing import Callable, Optional, TypeVar

import biotite.structure as bts
import pandas as pd

from volumizer import align, fib_sphere, pdb, utils, voxel
from volumizer.types import Annotation


_T = TypeVar("_T")


def _timed_call(
    stage_timings: dict[str, float] | None,
    stage_name: str,
    fn: Callable[..., _T],
    *args,
    **kwargs,
) -> _T:
    """
    Execute a callable and optionally accumulate elapsed time under `stage_name`.
    """
    if stage_timings is None:
        return fn(*args, **kwargs)

    start = perf_counter()
    result = fn(*args, **kwargs)
    stage_timings[stage_name] = stage_timings.get(stage_name, 0.0) + (
        perf_counter() - start
    )
    return result


def annotate_structure_volumes(
    structure: bts.AtomArray,
    min_voxels: Optional[int] = 2,
    min_volume: Optional[float] = None,
    stage_timings: dict[str, float] | None = None,
) -> tuple[pd.DataFrame, bts.AtomArray]:
    """
    Perform analysis of a prepared structure.
    """
    coords = _timed_call(
        stage_timings,
        "get_structure_coords",
        pdb.get_structure_coords,
        structure,
    )
    coords = _timed_call(
        stage_timings,
        "add_extra_points",
        fib_sphere.add_extra_points,
        coords,
        utils.VOXEL_SIZE,
    )

    cloud = _timed_call(
        stage_timings,
        "coords_to_point_cloud",
        voxel.coords_to_point_cloud,
        coords,
    )
    cloud, voxel_grid_id = _timed_call(
        stage_timings,
        "add_voxel_grid",
        voxel.add_voxel_grid,
        cloud,
    )
    voxel_grid = _timed_call(
        stage_timings,
        "get_voxel_grid",
        voxel.get_voxel_grid,
        cloud,
        voxel_grid_id,
    )

    protein_solvent_voxels = _timed_call(
        stage_timings,
        "get_protein_solvent_voxel_array",
        voxel.get_protein_solvent_voxel_array,
        voxel_grid,
    )
    protein_voxels, solvent_voxels = _timed_call(
        stage_timings,
        "get_protein_and_solvent_voxels",
        voxel.get_protein_and_solvent_voxels,
        protein_solvent_voxels,
        voxel_grid.x_y_z,
    )

    exposed_voxels, buried_voxels = _timed_call(
        stage_timings,
        "get_exposed_and_buried_voxels",
        voxel.get_exposed_and_buried_voxels,
        solvent_voxels,
        protein_voxels,
        voxel_grid.x_y_z,
    )

    # Subset exposed voxels into those directly neighboring a buried voxel.
    first_shell_exposed_voxels = _timed_call(
        stage_timings,
        "get_first_shell_exposed_voxels",
        voxel.get_first_shell_exposed_voxels,
        exposed_voxels,
        buried_voxels,
        voxel_grid,
    )

    if stage_timings is None:
        hubs, pores, pockets, cavities, occluded = (
            voxel.get_pores_pockets_cavities_occluded(
                buried_voxels,
                first_shell_exposed_voxels,
                voxel_grid,
            )
        )
    else:
        hubs, pores, pockets, cavities, occluded = (
            voxel.get_pores_pockets_cavities_occluded(
                buried_voxels,
                first_shell_exposed_voxels,
                voxel_grid,
                stage_timings=stage_timings,
            )
        )


    def _filter_and_sort_groups(group_dict):
        return utils.sort_voxelgroups_by_volume(
            utils.filter_voxelgroups_by_volume(
                group_dict,
                min_voxels=min_voxels,
                min_volume=min_volume,
            )
        )

    hubs = _timed_call(stage_timings, "filter_sort_hubs", _filter_and_sort_groups, hubs)
    pores = _timed_call(stage_timings, "filter_sort_pores", _filter_and_sort_groups, pores)
    pockets = _timed_call(
        stage_timings,
        "filter_sort_pockets",
        _filter_and_sort_groups,
        pockets,
    )
    cavities = _timed_call(
        stage_timings,
        "filter_sort_cavities",
        _filter_and_sort_groups,
        cavities,
    )
    occluded = _timed_call(
        stage_timings,
        "filter_sort_occluded",
        _filter_and_sort_groups,
        occluded,
    )

    def _build_annotation() -> Annotation:
        return Annotation(
            total_hub_volume=utils.get_volume_summary(hubs, "total"),
            total_pore_volume=utils.get_volume_summary(pores, "total"),
            total_cavity_volume=utils.get_volume_summary(cavities, "total"),
            total_pocket_volume=utils.get_volume_summary(pockets, "total"),
            largest_hub_volume=utils.get_volume_summary(hubs, "max"),
            largest_pore_volume=utils.get_volume_summary(pores, "max"),
            largest_cavity_volume=utils.get_volume_summary(cavities, "max"),
            largest_pocket_volume=utils.get_volume_summary(pockets, "max"),
            num_hubs=len(hubs),
            num_pores=len(pores),
            num_cavities=len(cavities),
            num_pockets=len(pockets),
            hub_volumes={i: hub.volume for i, hub in hubs.items()},
            pore_volumes={i: pore.volume for i, pore in pores.items()},
            cavity_volumes={i: cavity.volume for i, cavity in cavities.items()},
            pocket_volumes={i: pocket.volume for i, pocket in pockets.items()},
            hub_dimensions={i: hub.axial_lengths for i, hub in hubs.items()},
            pore_dimensions={i: pore.axial_lengths for i, pore in pores.items()},
            cavity_dimensions={i: cavity.axial_lengths for i, cavity in cavities.items()},
            pocket_dimensions={i: pocket.axial_lengths for i, pocket in pockets.items()},
        )

    annotation = _timed_call(stage_timings, "build_annotation", _build_annotation)
    annotation_df = _timed_call(
        stage_timings,
        "make_annotation_dataframe",
        utils.make_annotation_dataframe,
        annotation,
    )

    annotation_structure = _timed_call(
        stage_timings,
        "volumes_to_structure",
        pdb.volumes_to_structure,
        voxel_grid,
        hubs,
        pores,
        pockets,
        cavities,
        occluded,
    )

    return annotation_df, annotation_structure


def prepare_pdb_structure(structure: bts.AtomArray) -> bts.AtomArray:
    """
    Prepares a biotite AtomArray object for analysis by:

    1. Cleaning extraneous atoms
    2. Aligning along the principal axis

    Returns a biotite AtomArray
    """
    cleaned_structure = pdb.clean_structure(structure)
    cleaned_aligned_structure = align.align_structure(cleaned_structure)

    return cleaned_aligned_structure


def volumize_structure(
    structure: bts.AtomArray,
    stage_timings: dict[str, float] | None = None,
) -> tuple[pd.DataFrame, bts.AtomArray, bts.AtomArray]:
    """
    Perform the volumizer pipeline on a biotite structure and return the:
    annotation_df, pdb_structure, and annotated_structure
    """
    prepared_structure = _timed_call(
        stage_timings,
        "prepare_structure",
        prepare_pdb_structure,
        structure,
    )

    annotation_df, annotation_structure = annotate_structure_volumes(
        prepared_structure,
        stage_timings=stage_timings,
    )

    return annotation_df, prepared_structure, annotation_structure


def volumize_pdb(
    pdb_file: Path,
    stage_timings: dict[str, float] | None = None,
) -> tuple[pd.DataFrame, bts.AtomArray, bts.AtomArray]:
    """
    Convenience function to volumize a PDB file.
    """
    pdb_structure = _timed_call(
        stage_timings,
        "load_structure",
        pdb.load_structure,
        pdb_file,
    )
    return volumize_structure(pdb_structure, stage_timings=stage_timings)


def volumize_pdb_and_save(
    in_pdb_file: Path,
    out_pdb_file: Path,
    out_annotation_df: Path,
) -> None:
    """
    Convenience function to perform volumization and save the output.
    """
    annotation_df, prepared_pdb_structure, annotation_structure = volumize_pdb(
        in_pdb_file
    )
    out_pdb_lines = pdb.make_volumized_pdb_lines(
        [prepared_pdb_structure, annotation_structure]
    )

    annotation_df.to_json(out_annotation_df)
    pdb.save_pdb_lines(out_pdb_lines, out_pdb_file)
