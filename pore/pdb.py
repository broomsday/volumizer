"""
Functions for parsing, cleaning, and modifying PDBs.
"""

from pathlib import Path
import tempfile
import itertools

from Bio.PDB import PDBParser, Select, Structure
from Bio.PDB.PDBIO import PDBIO
from pyntcloud.structures.voxelgrid import VoxelGrid
import pandas as pd
import numpy as np

from pore.constants import VOXEL_ATOM_NAMES, VOXEL_TYPE_CHAIN_MAP, VOXEL_TYPE_ATOM_MAP
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
    voxel_type: str,
    resnum: int,
    voxel_index: int,
    point: np.ndarray,
) -> str:
    """
    Make a single atom line for a single voxel
    """
    return f"ATOM  {voxel_index:>5d}  {VOXEL_TYPE_ATOM_MAP[voxel_type]}   {voxel_type} {VOXEL_TYPE_CHAIN_MAP[voxel_type]}{resnum:>4d}    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00  0.00           {VOXEL_TYPE_ATOM_MAP[voxel_type]}"


def make_atom_lines(
    voxel_type: str,
    resnum: int,
    voxel_indices: set[int],
    voxel_grid_centers: np.ndarray,
) -> list[str]:
    """
    Make all atom lines for a given set of voxel indices
    """
    return [
        make_atom_line(voxel_type, resnum, voxel_index, voxel_grid_centers[voxel_index])
        for voxel_index in voxel_indices
    ]


def points_to_pdb(
    voxel_grid: VoxelGrid,
    pores: dict[int, VoxelGroup],
    pockets: dict[int, VoxelGroup],
    cavities: dict[int, VoxelGroup],
    occluded: dict[int, VoxelGroup],
) -> list[str]:
    """
    Write out voxels of our volumes as though they were atoms in a PDB file.
    """
    pore_voxel_indices = {i: voxels.indices for i, voxels in pores.items()}
    pocket_voxel_indices = {i: voxels.indices for i, voxels in pockets.items()}
    cavity_voxel_indices = {i: voxels.indices for i, voxels in cavities.items()}
    occluded_voxel_indices = {i: voxels.indices for i, voxels in occluded.items()}

    return list(
        itertools.chain(
            *[
                make_atom_lines("POK", voxel_group_index, voxel_indices, voxel_grid.voxel_centers)
                for voxel_group_index, voxel_indices in pocket_voxel_indices.items()
            ],
            *[
                make_atom_lines("POR", voxel_group_index, voxel_indices, voxel_grid.voxel_centers)
                for voxel_group_index, voxel_indices in pore_voxel_indices.items()
            ],
            *[
                make_atom_lines("CAV", voxel_group_index, voxel_indices, voxel_grid.voxel_centers)
                for voxel_group_index, voxel_indices in cavity_voxel_indices.items()
            ],
            *[
                make_atom_lines("OCC", voxel_group_index, voxel_indices, voxel_grid.voxel_centers)
                for voxel_group_index, voxel_indices in occluded_voxel_indices.items()
            ],
        )
    )
