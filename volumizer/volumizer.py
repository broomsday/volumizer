"""
Entry-level functions to find pores and cavities.
"""


from pathlib import Path
from typing import Optional

import biotite.structure as bts
import pandas as pd

from volumizer import utils, pdb, voxel, fib_sphere, align
from volumizer.types import Annotation


def annotate_structure_volumes(
    structure: bts.AtomArray,
    min_voxels: Optional[int] = 2,
    min_volume: Optional[float] = None,
) -> tuple[Annotation, list[str]]:
    """
    Perform analysis of a prepared structure.
    """
    coords = pdb.get_structure_coords(structure)
    coords = fib_sphere.add_extra_points(coords, utils.VOXEL_SIZE)

    cloud = voxel.coords_to_point_cloud(coords)
    cloud, voxel_grid_id = voxel.add_voxel_grid(cloud)
    voxel_grid = voxel.get_voxel_grid(cloud, voxel_grid_id)

    protein_solvent_voxels = voxel.get_protein_solvent_voxel_array(voxel_grid)
    protein_voxels, solvent_voxels = voxel.get_protein_and_solvent_voxels(
        protein_solvent_voxels, voxel_grid.x_y_z
    )

    exposed_voxels, buried_voxels = voxel.get_exposed_and_buried_voxels(
        solvent_voxels, protein_voxels, voxel_grid.x_y_z
    )

    # get sub-selection of exposed voxels that are NEXT to a buried voxel
    first_shell_exposed_voxels = voxel.get_first_shell_exposed_voxels(
        exposed_voxels, buried_voxels, voxel_grid
    )

    (
        hubs,
        pores,
        pockets,
        cavities,
        occluded,
    ) = voxel.get_pores_pockets_cavities_occluded(
        buried_voxels, first_shell_exposed_voxels, voxel_grid
    )
    hubs = utils.sort_voxelgroups_by_volume(
        utils.filter_voxelgroups_by_volume(
            hubs, min_voxels=min_voxels, min_volume=min_volume
        )
    )
    pores = utils.sort_voxelgroups_by_volume(
        utils.filter_voxelgroups_by_volume(
            pores, min_voxels=min_voxels, min_volume=min_volume
        )
    )
    pockets = utils.sort_voxelgroups_by_volume(
        utils.filter_voxelgroups_by_volume(
            pockets, min_voxels=min_voxels, min_volume=min_volume
        )
    )
    cavities = utils.sort_voxelgroups_by_volume(
        utils.filter_voxelgroups_by_volume(
            cavities, min_voxels=min_voxels, min_volume=min_volume
        )
    )
    occluded = utils.sort_voxelgroups_by_volume(
        utils.filter_voxelgroups_by_volume(
            occluded, min_voxels=min_voxels, min_volume=min_volume
        )
    )

    annotation = Annotation(
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
    annotation_df = utils.make_annotation_dataframe(annotation)

    annotation_structure = pdb.volumes_to_structure(
        voxel_grid, hubs, pores, pockets, cavities, occluded
    )

    return (annotation_df, annotation_structure)


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
) -> tuple[pd.DataFrame, bts.AtomArray, bts.AtomArray]:
    """
    Perform the volumizer pipeline on a biotite structure and return the:
    annotation_df, pdb_structure, and annotated_structure
    """
    prepared_structure = prepare_pdb_structure(structure)

    annotation_df, annotation_structure = annotate_structure_volumes(prepared_structure)

    return annotation_df, prepared_structure, annotation_structure


def volumize_pdb(pdb_file: Path) -> tuple[pd.DataFrame, bts.AtomArray, bts.AtomArray]:
    """
    Convenience function to volumize a PDB file.
    """
    pdb_structure = pdb.load_structure(pdb_file)
    return volumize_structure(pdb_structure)


def volumize_pdb_and_save(
    in_pdb_file: Path, out_pdb_file: Path, out_annotation_df: Path
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
