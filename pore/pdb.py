"""
Functions for parsing, cleaning, and modifying PDBs.
"""

import warnings

from pymongo import database
from tqdm import tqdm
from Bio.PDB import PDBParser, Select
from Bio.PDB.PDBIO import PDBIO
from pyntcloud import PyntCloud
from pyntcloud.structures.voxelgrid import VoxelGrid
import pandas as pd
import numpy as np
from tqdm import tqdm

from pore import utils, mongo
from pore.paths import CLEANED_PDB_DIR, DATA_DIR, PROCESSED_PDB_DIR
from pore.constants import VOXEL_ATOM_NAMES, VOXEL_SIZE, OCCLUDED_DIMENSION_LIMIT


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


def get_voxel_grid(cloud: PyntCloud, voxel_grid_id: str) -> VoxelGrid:
    """
    Generate an array representing a binary voxel grid where zero represents no protein
    atoms in that voxel (e.g. solvent) and a non-zero represents a protein atom in that voxel.
    """
    return cloud.structures[voxel_grid_id]


def get_protein_solvent_voxels(voxel_grid: VoxelGrid) -> np.ndarray:
    """
    Generate an array representing a binary voxel grid where zero represents no protein
    atoms in that voxel (e.g. solvent) and a non-zero represents a protein atom in that voxel.
    """
    return voxel_grid.get_feature_vector(mode="binary")


def get_protein_and_solvent_voxels(
    binary_voxel_grid: np.ndarray,
) -> tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
    """
    From the voxel grid, return the indices of the protein containing voxels,
    and those not containing protein (e.g. solvent).
    """
    protein_voxels = np.nonzero(binary_voxel_grid)
    solvent_voxels = np.nonzero(binary_voxel_grid == 0)

    return protein_voxels, solvent_voxels


def get_num_occluded_dimensions(
    query_voxel: tuple[int, int, int],
    possible_occluding_x: list[int],
    possible_occluding_y: list[int],
    possible_occluding_z: list[int],
) -> int:
    """
    Determine how many ordinal axes are occluded.
    """
    num_occluded_dimensions = 0

    if len(possible_occluding_z) != 0:
        if min(possible_occluding_z) < query_voxel[0]:
            num_occluded_dimensions += 1
        elif max(possible_occluding_z) > query_voxel[0]:
            num_occluded_dimensions += 1
    if len(possible_occluding_y) != 0:
        if min(possible_occluding_y) < query_voxel[0]:
            num_occluded_dimensions += 1
        elif max(possible_occluding_y) > query_voxel[0]:
            num_occluded_dimensions += 1
    if len(possible_occluding_x) != 0:
        if min(possible_occluding_x) < query_voxel[0]:
            num_occluded_dimensions += 1
        elif max(possible_occluding_x) > query_voxel[0]:
            num_occluded_dimensions += 1

    return num_occluded_dimensions


def get_exposed_and_buried_solvent_voxels(
    solvent_voxels: tuple[np.ndarray, ...],
    protein_voxels: tuple[np.ndarray, ...],
) -> tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
    """
    Use simple geometric heuristics to determine if a given solvent voxel is buried or exposed.
    """
    buried_solvent_voxels = ([], [], [])
    exposed_solvent_voxels = ([], [], [])

    previous_voxel = (-1, -1, -1)
    for i in tqdm(range(solvent_voxels[0].size), desc="Searching voxels"):
        # TODO numpy.where() is somewhat slow, might be a way to improve
        query_voxel = (solvent_voxels[0][i], solvent_voxels[1][i], solvent_voxels[2][i])
        if query_voxel[0] != previous_voxel[0]:
            match_x_protein_voxel_indices = set(np.where(protein_voxels[0] == query_voxel[0])[0])
        if query_voxel[1] != previous_voxel[1]:
            match_y_protein_voxel_indices = set(np.where(protein_voxels[1] == query_voxel[1])[0])
        if query_voxel[2] != previous_voxel[2]:
            match_z_protein_voxel_indices = set(np.where(protein_voxels[2] == query_voxel[2])[0])
        previous_voxel = query_voxel

        possible_occluding_z_indices = match_x_protein_voxel_indices.intersection(match_y_protein_voxel_indices)
        possible_occluding_y_indices = match_x_protein_voxel_indices.intersection(match_z_protein_voxel_indices)
        possible_occluding_x_indices = match_y_protein_voxel_indices.intersection(match_z_protein_voxel_indices)

        possible_occluding_z = [protein_voxels[2][i] for i in possible_occluding_z_indices]
        possible_occluding_y = [protein_voxels[1][i] for i in possible_occluding_y_indices]
        possible_occluding_x = [protein_voxels[0][i] for i in possible_occluding_x_indices]

        num_occluded_dimensions = get_num_occluded_dimensions(
            query_voxel, possible_occluding_x, possible_occluding_y, possible_occluding_z
        )

        if num_occluded_dimensions >= OCCLUDED_DIMENSION_LIMIT:
            buried_solvent_voxels[0].append(query_voxel[0])
            buried_solvent_voxels[1].append(query_voxel[1])
            buried_solvent_voxels[2].append(query_voxel[2])
        else:
            exposed_solvent_voxels[0].append(query_voxel[0])
            exposed_solvent_voxels[1].append(query_voxel[1])
            exposed_solvent_voxels[2].append(query_voxel[2])

    return (
        (
            np.array(exposed_solvent_voxels[0]),
            np.array(exposed_solvent_voxels[1]),
            np.array(exposed_solvent_voxels[2]),
        ),
        (
            np.array(buried_solvent_voxels[0]),
            np.array(buried_solvent_voxels[1]),
            np.array(buried_solvent_voxels[2]),
        ),
    )


def make_atom_line(
    point: np.ndarray,
    exposed_solvent_voxel_indices: set[int],
    buried_solvent_voxel_indices: set[int],
    index: int,
) -> str:
    """
    Make each exposed solvent point a hydrogen, buried solvent ponit an oxygen, and protein point a carbon
    """
    if index in exposed_solvent_voxel_indices:
        return f"ATOM  {index:>5d}  H   EXP A   1    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00  0.00           H"
    elif index in buried_solvent_voxel_indices:
        return f"ATOM  {index:>5d}  O   BUR B   1    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00  0.00           O"
    else:
        return f"ATOM  {index:>5d}  C   PTN C   1    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00  0.00           C"


def compute_voxel_indices(voxels: tuple[np.ndarray, ...], grid_dimensions: np.ndarray) -> set[int]:
    """
    Given a 3D array of voxels, compute the 1D set of their indices.
    """
    return {
        (voxels[0][i] * grid_dimensions[1] * grid_dimensions[2] + voxels[1][i] * grid_dimensions[2] + voxels[2][i])
        for i in range(len(voxels[0]))
    }


def points_to_pdb(
    voxel_grid: VoxelGrid,
    exposed_solvent_voxels: tuple[np.ndarray, ...],
    buried_solvent_voxels: tuple[np.ndarray, ...],
) -> None:
    """
    Write out points as though it was a PDB file.

    Carbon -> protein
    Nitrogen -> exposed solvent
    Oxygen -> buried solvent
    """
    exposed_solvent_voxel_indices = compute_voxel_indices(exposed_solvent_voxels, voxel_grid.x_y_z)
    buried_solvent_voxel_indices = compute_voxel_indices(buried_solvent_voxels, voxel_grid.x_y_z)

    output_lines = [
        make_atom_line(point, exposed_solvent_voxel_indices, buried_solvent_voxel_indices, i)
        for i, point in enumerate(voxel_grid.voxel_centers)
    ]

    with open(DATA_DIR / "tmp.pdb", mode="w", encoding="utf-8") as pdb_file:
        pdb_file.write("\n".join(output_lines))


def process_one_pdb(pdb_id: str) -> bool:
    """
    Generate a voxelized representation of the protein.

    Label voxels as:
    a) Protein
    b) Buried Solvent (potentially a pore)
    c) Exposed Solvent

    Save the voxelized represenation.
    """
    pdb_id = "6MRT"  # TODO: temporary hard-coding of structure that should have obvious pore

    coords = get_pdb_coords(pdb_id)
    cloud = coords_to_point_cloud(coords)
    cloud, voxel_grid_id = add_voxel_grid(cloud)
    voxel_grid = get_voxel_grid(cloud, voxel_grid_id)
    protein_solvent_voxels = get_protein_solvent_voxels(voxel_grid)  # TODO simplify this
    protein_voxels, solvent_voxels = get_protein_and_solvent_voxels(protein_solvent_voxels)

    exposed_solvent_voxels, buried_solvent_voxels = get_exposed_and_buried_solvent_voxels(
        solvent_voxels, protein_voxels
    )

    print("\n")
    print(pdb_id)
    print(voxel_grid_id)
    print("Protein Voxels:", protein_voxels[0].size)
    print("Total Solvent Voxels:", solvent_voxels[0].size)
    print("Exposed Solvent Voxels:", exposed_solvent_voxels[0].size)
    print("Buried Solvent Voxels:", buried_solvent_voxels[0].size)

    points_to_pdb(voxel_grid, exposed_solvent_voxels, buried_solvent_voxels)

    quit()
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
