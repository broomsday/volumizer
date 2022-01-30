"""
Functions for parsing, cleaning, and modifying PDBs.
"""

from pathlib import Path
import tempfile
import itertools
from typing import Optional

from Bio.PDB import PDBParser, Select, Structure
from Bio.PDB.PDBIO import PDBIO
from pyntcloud.structures.voxelgrid import VoxelGrid
import pandas as pd
import numpy as np

from pore.constants import VOXEL_ATOM_NAMES, VOXEL_TYPE_CHAIN_MAP, VOXEL_TYPE_ATOM_MAP
from pore.paths import PREPARED_PDB_DIR, ANNOTATED_PDB_DIR
from pore.types import VoxelGroup, Annotation
from pore import utils


PDB_IN = PDBParser()
PDB_OUT = PDBIO()


class ProteinSelect(Select):
    def __init__(self, components):
        Select.__init__(self)
        self.components = components

    def accept_model(self, model):
        if (model.serial_num == 0) or (model.serial_num == 1):
            return True
        return utils.KEEP_MODELS

    def accept_residue(self, residue):
        if residue.resname in self.components:
            return True
        return utils.KEEP_NON_PROTEIN

    def accept_atom(self, atom):
        if (not atom.is_disordered()) or (atom.get_altloc() == "A"):
            atom.set_altloc(" ")

        if atom.element != "H":
            return True
        return utils.KEEP_HYDROGENS


def save_pdb(structure: Structure, pdb_file: Path, remarks: Optional[str] = None) -> None:
    """
    Save a biopython PDB structure to a PDB file.
    """
    PDB_OUT.set_structure(structure)
    PDB_OUT.save(str(pdb_file))

    # if we are adding remarks, read in the just written file, add the remarks and overwrite
    if remarks is not None:
        with open(pdb_file, mode="r", encoding="utf-8") as in_file:
            lines = in_file.readlines()

        with open(pdb_file, mode="w", encoding="utf-8") as out_file:
            out_file.write(remarks + "".join(lines))


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


def get_structure_coords(structure: Structure, canonical_protein_atoms_only: bool = False) -> pd.DataFrame:
    """
    Load the PDB file using Biopython.

    Return the coordinates of backbone heavy atoms and CB atoms as a dataframe.
    """
    if canonical_protein_atoms_only:
        coordinates = [atom.coord for atom in structure.get_atoms() if atom.name in VOXEL_ATOM_NAMES]
        elements = [atom.element for atom in structure.get_atoms() if atom.name in VOXEL_ATOM_NAMES]
    else:
        coordinates = [atom.coord for atom in structure.get_atoms()]
        elements = [atom.element for atom in structure.get_atoms()]

    coordinates = pd.DataFrame(coordinates, columns=["x", "y", "z"])
    coordinates["element"] = elements

    return coordinates


def make_atom_line(
    voxel_type: str,
    resnum: int,
    voxel_index: int,
    point: np.ndarray,
    beta: float,
) -> str:
    """
    Make a single atom line for a single voxel
    """
    return f"ATOM  {voxel_index:>5d}  {VOXEL_TYPE_ATOM_MAP[voxel_type]}   {voxel_type} {VOXEL_TYPE_CHAIN_MAP[voxel_type]}{resnum:>4d}    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00{beta:>6.2f}           {VOXEL_TYPE_ATOM_MAP[voxel_type]}"


def make_atom_lines(
    voxel_type: str,
    resnum: int,
    voxel_indices: set[int],
    surface_indices: set[int],
    voxel_grid_centers: np.ndarray,
) -> list[str]:
    """
    Make all atom lines for a given set of voxel indices
    """
    if voxel_type == "POR":
        print(voxel_indices)
        print(surface_indices)

    return [
        make_atom_line(voxel_type, resnum, voxel_index, voxel_grid_centers[voxel_index], 50.0)
        if voxel_index in surface_indices
        else make_atom_line(voxel_type, resnum, voxel_index, voxel_grid_centers[voxel_index], 0.0)
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
    return list(
        itertools.chain(
            *[
                make_atom_lines("POK", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in pockets.items()
            ],
            *[
                make_atom_lines("POR", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in pores.items()
            ],
            *[
                make_atom_lines("CAV", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in cavities.items()
            ],
            *[
                make_atom_lines("OCC", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in occluded.items()
            ],
        )
    )


def save_annotated_pdb(pdb_name: str, annotated_lines: list[str]) -> None:
    """
    Save a PDB formatted coordinate file of the voxels.
    Individual atoms/voxels are labelled according to type
    """
    # get the original PDB lines so that our annotation can be appended
    with open(PREPARED_PDB_DIR / f"{pdb_name}.pdb", mode="r", encoding="utf-8") as input_pdb_file:
        pdb_lines = input_pdb_file.readlines()
    pdb_lines = [line.rstrip("\n") for line in pdb_lines]

    # remove the terminal END line (NOTE: also removes terminal ENDMDL line)
    for line in pdb_lines[::-1]:
        if line.strip() == "":
            pdb_lines.pop()
        if ("ENDMDL" not in line) and ("END" in line):
            pdb_lines.pop()
            break

    # add the annotated lines and save
    with open(ANNOTATED_PDB_DIR / f"{pdb_name}.pdb", mode="w", encoding="utf-8") as annotated_pdb_file:
        annotated_pdb_file.write("\n".join([*pdb_lines, "END", *annotated_lines]))


def generate_rotation_translation_remarks(rotation: np.ndarray, translation: np.ndarray) -> str:
    """
    Given a rotation matrix and translation vector, generate a string allowing them to be
    written as REMARKS to a PDB file.
    """
    return (
        f"REMARK TRANSLATION {str(translation)}\n"
        f"REMARK ROTATION [{str(rotation[0])}, {str(rotation[1])}, {str(rotation[2])}]\n"
    )


def generate_volume_remarks(annotation: Annotation) -> list[str]:
    """
    Return REMARKS for a PDB file to hold all volume information from an annotation object.
    """
    pore_remarks = [
        f"REMARK PORE VOLUME {annotation.pore_volumes[i]} DIMENSIONS {annotation.pore_dimensions[i][0], annotation.pore_dimensions[i][1], annotation.pore_dimensions[i][2]}"
        for i in annotation.pore_volumes
    ]

    cavity_remarks = [
        f"REMARK CAVITY VOLUME {annotation.cavity_volumes[i]} DIMENSIONS {annotation.cavity_dimensions[i][0], annotation.cavity_dimensions[i][1], annotation.cavity_dimensions[i][2]}"
        for i in annotation.cavity_volumes
    ]

    pocket_remarks = [
        f"REMARK POCKET VOLUME {annotation.pocket_volumes[i]} DIMENSIONS {annotation.pocket_dimensions[i][0], annotation.pocket_dimensions[i][1], annotation.pocket_dimensions[i][2]}"
        for i in annotation.pocket_volumes
    ]

    return pore_remarks + cavity_remarks + pocket_remarks
