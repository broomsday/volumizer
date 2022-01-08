"""
Entry-level functions to find pores and cavities.
"""


from pathlib import Path

from Bio.PDB import Structure

from pore import rcsb, utils, pdb, voxel, fib_sphere
from pore.types import Annotation
from pore.paths import PREPARED_PDB_DIR


def annotate_pdb_structure(structure: Structure) -> tuple[Annotation, list[str]]:
    """
    Perform analysis of a prepared structure.
    """
    coords = pdb.get_structure_coords(structure)
    # TODO here need to get the elements for the coords, or in the above function even
    coords = fib_sphere.add_extra_points(coords, utils.VOXEL_SIZE)

    cloud = voxel.coords_to_point_cloud(coords)
    cloud, voxel_grid_id = voxel.add_voxel_grid(cloud)
    voxel_grid = voxel.get_voxel_grid(cloud, voxel_grid_id)

    protein_solvent_voxels = voxel.get_protein_solvent_voxel_array(voxel_grid)
    protein_voxels, solvent_voxels = voxel.get_protein_and_solvent_voxels(protein_solvent_voxels, voxel_grid.x_y_z)

    exposed_voxels, buried_voxels = voxel.get_exposed_and_buried_voxels(
        solvent_voxels, protein_voxels, voxel_grid.x_y_z
    )

    pores, pockets, cavities, occluded = voxel.get_pores_pockets_cavities_occluded(buried_voxels, exposed_voxels, voxel_grid.x_y_z)
    pores = utils.sort_voxelgroups_by_volume(pores)
    pockets = utils.sort_voxelgroups_by_volume(pockets)
    cavities = utils.sort_voxelgroups_by_volume(cavities)
    occluded = utils.sort_voxelgroups_by_volume(occluded)

    pdb_lines = pdb.points_to_pdb(voxel_grid, exposed_voxels, buried_voxels, pores, pockets, cavities, occluded)

    annotation = Annotation(
        total_pore_volume=utils.get_volume_summary(pores, "total"),
        total_cavity_volume=utils.get_volume_summary(cavities, "total"),
        total_pocket_volume=utils.get_volume_summary(pockets, "total"),
        largest_pore_volume=utils.get_volume_summary(pores, "max"),
        largest_cavity_volume=utils.get_volume_summary(cavities, "max"),
        largest_pocket_volume=utils.get_volume_summary(pockets, "max"),
        num_pores=len(pores),
        num_cavities=len(cavities),
        num_pockets=len(pockets),
        pore_volumes={i: pore.volume for i, pore in pores.items()},
        cavity_volumes={i: cavity.volume for i, cavity in cavities.items()},
        pocket_volumes={i: pocket.volume for i, pocket in pockets.items()},
    )

    return (annotation, pdb_lines)


def prepare_pdb_structure(structure: Structure) -> Structure:
    """
    Prepares a biopython Structure object for analysis by:

    1. Cleaning extraneous atoms
    2. Adding hydrogens
    3. Aligning along the principal axis

    Returns a biopython Structure object.
    """
    protein_components = utils.load_protein_components()

    structure = pdb.clean_structure(structure, protein_components)

    # TODO align

    return structure


def prepare_pdb_file(pdb_file: Path) -> Structure:
    """
    Prepares a PDB file at the given path for analysis by:

    1. Cleaning extraneous atoms
    2. Adding hydrogens
    3. Aligning along the principal axis

    Returns a biopython Structure object.
    """
    return prepare_pdb_structure(pdb.load_pdb(pdb_file))


def process_pdb_file(pdb_file: Path) -> tuple[Annotation, list[str]]:
    """
    Perform end-to-end pipeline on a single PDB file.
    """
    prepared_structure = prepare_pdb_file(pdb_file)
    pdb.save_pdb(prepared_structure, PREPARED_PDB_DIR / pdb_file.name)

    return annotate_pdb_structure(prepared_structure)


def download_pdb_file(pdb_id: str) -> Path:
    """
    Download the biological assembly from the RCSB.
    Unzip and save the PDB.
    """
    if not utils.is_pdb_downloaded(pdb_id):
        rcsb.download_biological_assembly(pdb_id)
        utils.decompress_pdb(pdb_id)

    return utils.get_downloaded_pdb_path(pdb_id)