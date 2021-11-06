"""
Entry-level functions to find pores and cavities.
"""


from pathlib import Path

from Bio.PDB import Structure

from pore import rcsb, utils, pdb, voxel
from pore.types import Annotation
from pore.paths import PREPARED_PDB_DIR, ANNOTATED_PDB_DIR


def annotate_pdb_structure(structure: Structure) -> Annotation:
    """
    Perform analysis of a prepared structure.
    """
    coords = pdb.get_structure_coords(structure)
    cloud = voxel.coords_to_point_cloud(coords)
    cloud, voxel_grid_id = voxel.add_voxel_grid(cloud)
    voxel_grid = voxel.get_voxel_grid(cloud, voxel_grid_id)
    protein_solvent_voxels = voxel.get_protein_solvent_voxels(voxel_grid)
    protein_voxels, solvent_voxels = voxel.get_protein_and_solvent_voxels(protein_solvent_voxels)

    exposed_solvent_voxels, buried_solvent_voxels = voxel.get_exposed_and_buried_solvent_voxels(
        solvent_voxels, protein_voxels, voxel_grid.x_y_z
    )
    pore_voxels, cavity_voxels = voxel.get_pore_and_cavity_voxels(buried_solvent_voxels, exposed_solvent_voxels)

    print("\n")
    print(voxel_grid_id)
    print("Protein Voxels:", protein_voxels[0].size)
    print("Total Solvent Voxels:", solvent_voxels[0].size)
    print("Exposed Solvent Voxels:", exposed_solvent_voxels[0].size)
    print("Buried Solvent Voxels:", buried_solvent_voxels[0].size)
    print("Pores:", len(pore_voxels))
    print("Cavities:", len(cavity_voxels))

    pdb_lines = pdb.points_to_pdb(voxel_grid, exposed_solvent_voxels, buried_solvent_voxels, pore_voxels, cavity_voxels)
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

    # TODO add hydrogens

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