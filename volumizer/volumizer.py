"""
Entry-level functions to find pores and cavities.
"""


from pathlib import Path
from typing import Optional

import biotite.structure as bts

from volumizer import rcsb, utils, pdb, voxel, fib_sphere, align
from volumizer.types import Annotation
from volumizer.paths import PREPARED_PDB_DIR


def annotate_pdb_structure(
    structure: bts.AtomArray, min_voxels: Optional[int] = 2, min_volume: Optional[float] = None
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
    protein_voxels, solvent_voxels = voxel.get_protein_and_solvent_voxels(protein_solvent_voxels, voxel_grid.x_y_z)

    exposed_voxels, buried_voxels = voxel.get_exposed_and_buried_voxels(
        solvent_voxels, protein_voxels, voxel_grid.x_y_z
    )

    # get sub-selection of exposed voxels that are NEXT to a buried voxel
    first_shell_exposed_voxels = voxel.get_first_shell_exposed_voxels(exposed_voxels, buried_voxels, voxel_grid)

    hubs, pores, pockets, cavities, occluded = voxel.get_pores_pockets_cavities_occluded(
        buried_voxels, first_shell_exposed_voxels, voxel_grid
    )
    hubs = utils.sort_voxelgroups_by_volume(
        utils.filter_voxelgroups_by_volume(hubs, min_voxels=min_voxels, min_volume=min_volume)
    )
    pores = utils.sort_voxelgroups_by_volume(
        utils.filter_voxelgroups_by_volume(pores, min_voxels=min_voxels, min_volume=min_volume)
    )
    pockets = utils.sort_voxelgroups_by_volume(
        utils.filter_voxelgroups_by_volume(pockets, min_voxels=min_voxels, min_volume=min_volume)
    )
    cavities = utils.sort_voxelgroups_by_volume(
        utils.filter_voxelgroups_by_volume(cavities, min_voxels=min_voxels, min_volume=min_volume)
    )
    occluded = utils.sort_voxelgroups_by_volume(
        utils.filter_voxelgroups_by_volume(occluded, min_voxels=min_voxels, min_volume=min_volume)
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

    pdb_lines = pdb.points_to_pdb(voxel_grid, hubs, pores, pockets, cavities, occluded)
    pdb_lines.extend(pdb.generate_resolution_remarks())
    pdb_lines.extend(pdb.generate_volume_remarks(annotation))

    return (annotation, pdb_lines)


def prepare_pdb_structure(structure: bts.AtomArray) -> tuple[bts.AtomArray, str]:
    """
    Prepares a biotite AtomArray object for analysis by:

    1. Cleaning extraneous atoms
    2. Aligning along the principal axis

    Returns a biotite AtomArray
    """
    protein_components = utils.load_protein_components()    # TODO: make this a single datafile or in constants.

    structure = pdb.clean_structure(structure, protein_components)

    structure, rotation, translation = align.align_structure(structure)
    remarks = pdb.generate_rotation_translation_remarks(rotation, translation)

    return structure, remarks


def prepare_pdb_file(pdb_file: Path) -> tuple[bts.AtomArray, str]:
    """
    Prepares a PDB file at the given path for analysis by:

    1. Cleaning extraneous atoms
    2. Adding hydrogens
    3. Aligning along the principal axis

    Returns a biopython Structure object.
    """
    return prepare_pdb_structure(pdb.load_pdb(pdb_file))


def download_pdb_file(pdb_id: str) -> Path:
    """
    Download the biological assembly from the RCSB.
    Unzip and save the PDB.
    """
    if not utils.is_pdb_downloaded(pdb_id):
        if rcsb.download_biological_assembly(pdb_id):
            utils.decompress_pdb(pdb_id)

    return utils.get_downloaded_pdb_path(pdb_id)


def process_pdb_file(pdb_file: Path) -> tuple[Annotation, list[str]]:
    """
    Perform end-to-end pipeline on a single PDB file.
    """
    # TODO: add optional parameter for output PDB path
    prepared_structure, remarks = prepare_pdb_file(pdb_file)
    pdb.save_pdb(prepared_structure, PREPARED_PDB_DIR / pdb_file.name, remarks=remarks)

    return annotate_pdb_structure(prepared_structure)
