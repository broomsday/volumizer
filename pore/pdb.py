"""
Functions for parsing, cleaning, and modifying PDBs.
"""

from pathlib import Path
import tempfile

from Bio.PDB import PDBParser, Select, Structure
from Bio.PDB.PDBIO import PDBIO
from pyntcloud.structures.voxelgrid import VoxelGrid
import pandas as pd
import numpy as np

from pore.constants import VOXEL_ATOM_NAMES
from pore.types import VoxelGroup


PDB_IN = PDBParser()
PDB_OUT = PDBIO()


class ProteinSelect(Select):
    def __init__(self, components):
        Select.__init__(self)
        self.components = components

    def accept_model(self, model):
        if (model.serial_num == 0) or (model.serial_num == 1):
            return True
        return False

    def accept_residue(self, residue):
        if residue.resname in self.components:
            return True
        return False

    def accept_atom(self, atom):
        if (not atom.is_disordered()) or (atom.get_altloc() == "A"):
            atom.set_altloc(" ")

        if atom.element != "H":
                return True
        return True


def save_pdb(structure: Structure, pdb_file: Path) -> None:
    """
    Save a biopython PDB structure to a PDB file.
    """
    PDB_OUT.set_structure(structure)
    PDB_OUT.save(str(pdb_file))


def load_pdb(pdb_file: Path) -> Structure:
    """
    Load a PDB file as a biopython PDB structure.
    """
    return PDB_IN.get_structure(pdb_file.stem, pdb_file)


def clean_structure(structure: Structure, protein_components: set[str]) -> Structure:
    """
    Clean a PDB by saving it to a temporary file with the select class and then reloading.
    """
    PDB_OUT.set_structure(structure)
    with tempfile.SpooledTemporaryFile(mode="w+") as pdb_file:
        PDB_OUT.save(pdb_file, select=ProteinSelect(protein_components))
        pdb_file.seek(0)
        structure = PDB_IN.get_structure(structure.id, pdb_file)

    return structure


def get_structure_coords(structure: Structure) -> pd.DataFrame:
    """
    Load the PDB file using Biopython.

    Return the coordinates of backbone heavy atoms and CB atoms as a dataframe.
    """
    # TODO used of VOXEL_ATOM_NAMES should be optional
    coordinates = [atom.coord for atom in structure.get_atoms() if atom.name in VOXEL_ATOM_NAMES]
    return pd.DataFrame(coordinates, columns=["x", "y", "z"])


def make_atom_line(
    point: np.ndarray,
    exposed_solvent_voxel_indices: set[int],
    buried_solvent_voxel_indices: set[int],
    pore_voxel_index_map: dict[int, int],
    pocket_voxel_index_map: dict[int, int],
    cavity_voxel_index_map: dict[int, int],
    occluded_voxel_index_map: dict[int, int],
    index: int,
) -> str:
    """
    Make each exposed solvent point a hydrogen, buried solvent ponit an oxygen, and protein point a carbon
    """
    if index in exposed_solvent_voxel_indices:
        return f"ATOM  {index:>5d}  H   EXP B   1    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00  0.00           H"
    elif index in buried_solvent_voxel_indices:
        if index in list(pore_voxel_index_map.keys()):
            return f"ATOM  {index:>5d}  O   POR E{pore_voxel_index_map[index]:>4d}    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00  0.00           O"
        elif index in list(pocket_voxel_index_map.keys()):
            return f"ATOM  {index:>5d}  N   POK D{pocket_voxel_index_map[index]:>4d}    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00  0.00           N"
        elif index in list(cavity_voxel_index_map.keys()):
            return f"ATOM  {index:>5d}  S   CAV F{cavity_voxel_index_map[index]:>4d}    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00  0.00           S"
        elif index in list(occluded_voxel_index_map.keys()):
            return f"ATOM  {index:>5d}  H   OCC C   1    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00  0.00           H"
    else:
        return f"ATOM  {index:>5d}  C   PTN A   1    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00  0.00           C"


def points_to_pdb(
    voxel_grid: VoxelGrid,
    exposed_voxels: VoxelGroup,
    buried_voxels: VoxelGroup,
    pores: dict[int, VoxelGroup],
    pockets: dict[int, VoxelGroup],
    cavities: dict[int, VoxelGroup],
    occluded: dict[int, VoxelGroup],
) -> list[str]:
    """
    Write out points as though it was a PDB file.

    Carbon -> protein
    Nitrogen -> exposed solvent
    Oxygen -> buried solvent
    """
    pore_voxel_indices = {i: voxels.indices for i, voxels in pores.items()}
    pocket_voxel_indices = {i: voxels.indices for i, voxels in pockets.items()}
    cavity_voxel_indices = {i: voxels.indices for i, voxels in cavities.items()}
    occluded_voxel_indices = {i: voxels.indices for i, voxels in occluded.items()}

    pore_voxel_index_map = {}
    for id, indices in pore_voxel_indices.items():
        for index in indices:
            pore_voxel_index_map[index] = id

    pocket_voxel_index_map = {}
    for id, indices in pocket_voxel_indices.items():
        for index in indices:
            pocket_voxel_index_map[index] = id

    cavity_voxel_index_map = {}
    for id, indices in cavity_voxel_indices.items():
        for index in indices:
            cavity_voxel_index_map[index] = id

    occluded_voxel_index_map = {}
    for id, indices in occluded_voxel_indices.items():
        for index in indices:
            occluded_voxel_index_map[index] = id

    return [
        make_atom_line(
            point,
            exposed_voxels.indices,
            buried_voxels.indices,
            pore_voxel_index_map,
            pocket_voxel_index_map,
            cavity_voxel_index_map,
            occluded_voxel_index_map,
            i,
        )
        for i, point in enumerate(voxel_grid.voxel_centers)
    ]

