"""
Entry-level functions to find pores and cavities.
"""


from pathlib import Path

from Bio.PDB import Structure

from pore import rcsb, utils, pdb, voxel, fib_sphere
from pore.constants import VOXEL_SIZE
from pore.types import Annotation
from pore.paths import PREPARED_PDB_DIR, ANNOTATED_PDB_DIR


def annotate_pdb_structure(structure: Structure) -> Annotation:
    """
    Perform analysis of a prepared structure.
    """
    coords = pdb.get_structure_coords(structure)
    # TODO here need to get the elements for the coords, or in the above function even
    coords = fib_sphere.add_extra_points(coords, VOXEL_SIZE)

    cloud = voxel.coords_to_point_cloud(coords)
    cloud, voxel_grid_id = voxel.add_voxel_grid(cloud)
    voxel_grid = voxel.get_voxel_grid(cloud, voxel_grid_id)

    protein_solvent_voxels = voxel.get_protein_solvent_voxel_array(voxel_grid)
    protein_voxels, solvent_voxels = voxel.get_protein_and_solvent_voxels(protein_solvent_voxels, voxel_grid.x_y_z)

    exposed_voxels, buried_voxels = voxel.get_exposed_and_buried_voxels(
        solvent_voxels, protein_voxels, voxel_grid.x_y_z
    )
    pores, cavities = voxel.get_pores_and_cavities(buried_voxels, exposed_voxels, voxel_grid.x_y_z)
    # TODO what about void voxels? e.g. pores that DO NOT span the box
        # these should actually be called 'cavities' and the current should be called 'pockets'
        # buried solvent that fails a minimium volume cutoff should just be reclassed as part of the bulk solvent

    print("\n")
    print(voxel_grid_id)
    print("Protein Voxels:", protein_voxels.num_voxels)
    print("Total Solvent Voxels:", solvent_voxels.num_voxels)
    print("Exposed Solvent Voxels:", exposed_voxels.num_voxels)
    print("Buried Solvent Voxels:", buried_voxels.num_voxels)
    print("Pores:", len(pores))
    print("Cavities:", len(cavities))

    pdb_lines = pdb.points_to_pdb(voxel_grid, exposed_voxels, buried_voxels, pores, cavities)
    with open(ANNOTATED_PDB_DIR / f"{structure.id}.pdb", mode="w", encoding="utf-8") as pdb_file:
        pdb_file.write("\n".join(pdb_lines))

    # TODO need to return an actual annotation, generated from the voxels

    quit()


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


def process_pdb_file(pdb_file: Path) -> Annotation:
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