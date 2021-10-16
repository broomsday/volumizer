"""
Functions for parsing, cleaning, and modifying PDBs.
"""

import warnings

from pymongo import database
from tqdm import tqdm
from Bio.PDB import PDBParser, Select
from Bio.PDB.PDBIO import PDBIO
from pyntcloud import PyntCloud
import pandas as pd
import numpy as np

from pore import utils, mongo
from pore.paths import CLEANED_PDB_DIR, PROCESSED_PDB_DIR
from pore.constants import VOXEL_ATOM_NAMES, VOXEL_SIZE


PDB_IN = PDBParser()
PDB_OUT = PDBIO()


class ProteinSelect(Select):
    def __init__(self, components):
        Select.__init__(self)
        self.components = components

    def accept_model(self, model):
        if model.serial_num == 1:
            return True
        return False

    def accept_residue(self, residue):
        if residue.resname in self.components:
            return True
        return False

    def accept_atom(self, atom):
        if atom.element != "H":
            return True
        if (not atom.is_disordered()) or (atom.get_altloc() == "A"):
            atom.set_altloc(" ")
            return True
        return False


def clean_one_pdb(pdb_id: str, protein_components: set[str]) -> bool:
    """
    Decompress and cleanup a single PDB.

    Return False if the final processed PDB could not be created.
    """
    utils.decompress_pdb(pdb_id)
    structure = PDB_IN.get_structure(pdb_id, utils.get_clean_pdb_path(pdb_id))

    PDB_OUT.set_structure(structure)
    PDB_OUT.save(str(utils.get_clean_pdb_path(pdb_id)), select=ProteinSelect(protein_components))

    return utils.is_pdb_cleaned(pdb_id)


def clean_all_pdbs(db: database.Database) -> None:
    """
    Decompress and cleanup all downloaded PDBs, and align to their major axis.

    Enter their state in the database.
    """
    CLEANED_PDB_DIR.mkdir(exist_ok=True, parents=True)

    protein_components = mongo.get_protein_components(db)

    # biopython tends to have many warnings about PDB construction which can be ignored
    warnings.filterwarnings("ignore")
    for pdb in tqdm(list(db.pdbs.find({"downloaded": True, "cleaned": False})), "Cleaning PDBs"):
        if utils.is_pdb_cleaned(pdb["pdb_id"]):
            mongo.update_pdb_cleaned(db, pdb, True)
        else:
            mongo.update_pdb_cleaned(db, pdb, clean_one_pdb(pdb["pdb_id"], protein_components))
    warnings.filterwarnings("default")


def get_pdb_coords(pdb_id: str) -> pd.DataFrame:
    """
    Load the PDB file using Biopython.

    Return the coordinates of backbone heavy atoms and CB atoms as a dataframe.
    """
    structure = PDB_IN.get_structure(pdb_id, utils.get_clean_pdb_path(pdb_id))
    coordinates = [atom.coord for atom in structure.get_atoms() if atom.name in VOXEL_ATOM_NAMES]
    return pd.DataFrame(coordinates, columns=["x", "y", "z"])


def coords_to_point_cloud(coords: pd.DataFrame) -> PyntCloud:
    """
    Produce a point-cloud from atomic coordinates.
    """
    return PyntCloud(coords)


def add_voxel_grid(cloud: PyntCloud) -> tuple[PyntCloud, str]:
    """
    Generate a voxel grid surrounding the point-cloud.
    """
    voxel_grid_id = cloud.add_structure(
        "voxelgrid", size_x=VOXEL_SIZE, size_y=VOXEL_SIZE, size_z=VOXEL_SIZE, regular_bounding_box=False
    )
    return cloud, voxel_grid_id


def create_protein_solvent_voxel_grid(cloud: PyntCloud, voxel_grid_id: str) -> np.ndarray:
    """
    Generate an array representing a binary voxel grid where zero represents no protein
    atoms in that voxel (e.g. solvent) and a non-zero represents a protein atom in that voxel.
    """
    voxel_grid = cloud.structures[voxel_grid_id]
    return voxel_grid.get_feature_vector(mode="binary")


def invert_binary_array(array: np.ndarray) -> np.ndarray:
    """
    Make all ones into zeros, and zeros into ones.

    This is only sensible if the array contains only zeros and ones
    """
    return (array * -1) + 1


def get_protein_and_solvent_voxels(
    binary_voxel_grid: np.ndarray,
) -> tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
    """
    From the voxel grid, return the indices of the protein containing voxels,
    and those not containing protein (e.g. solvent).
    """
    protein_voxels = binary_voxel_grid.nonzero()
    solvent_voxels = invert_binary_array(binary_voxel_grid).nonzero()

    return protein_voxels, solvent_voxels


def process_one_pdb(pdb_id: str) -> bool:
    """
    Generate a voxelized representation of the protein.

    Label voxels as:
    a) Protein
    b) Buried Solvent (potentially a pore)
    c) Exposed Solvent

    Save the voxelized represenation.
    """

    coords = get_pdb_coords(pdb_id)
    cloud = coords_to_point_cloud(coords)
    cloud, voxel_grid_id = add_voxel_grid(cloud)
    protein_solvent_voxel_grid = create_protein_solvent_voxel_grid(cloud, voxel_grid_id)
    protein_voxels, solvent_voxels = get_protein_and_solvent_voxels(protein_solvent_voxel_grid)

    # TODO sub-divide the solvent voxels into exposed vs. buried

    # TODO save the voxelized representation (may need colors to identify different label types)

    return False


def process_all_pdbs(db: database.Database) -> None:
    """
    Generate a voxelized representation of all proteins in the database.

    Label voxels as:
    a) Protein
    b) Buried Solvent (potentially a pore)
    c) Exposed Solvent

    Save the voxelized representations and labels for each protein.
    """
    PROCESSED_PDB_DIR.mkdir(exist_ok=True, parents=True)

    for pdb in tqdm(db.pdbs.find({"downloaded": True, "cleaned": True, "processed": False})):
        if utils.is_pdb_processed(pdb["pdb_id"]):
            mongo.update_pdb_processed(db, pdb, True)
        else:
            mongo.update_pdb_processed(db, pdb, process_one_pdb(pdb["pdb_id"]))
