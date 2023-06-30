"""
Functions for parsing, cleaning, and modifying PDBs.
"""

from pathlib import Path
import itertools

import biotite.structure as bts
from biotite.structure.io import load_structure # in use by modules importing pdb.py
from pyntcloud.structures.voxelgrid import VoxelGrid
import pandas as pd
import numpy as np

from volumizer.constants import (
    VOXEL_TYPE_CHAIN_MAP,
    VOXEL_TYPE_ATOM_MAP,
)
from volumizer.types import VoxelGroup, Annotation
from volumizer import utils


def save_pdb_lines(pdb_lines: list[str], output: Path) -> None:
    """
    Save individual pdb lines to a file.
    """
    with open(output, mode="w", encoding="utf-8") as out_file:
        out_file.writelines("%s\n" % line for line in pdb_lines)


def clean_structure(structure: bts.AtomArray) -> bts.AtomArray:
    """
    Clean the AtomArray of a PDB based on selected preferences.
    """
    if isinstance(structure, bts.AtomArrayStack):
        if not utils.KEEP_MODELS:
            structure = structure[0]
        elif structure.stack_depth() > 1:
            new_structure = structure[0]
            for model in structure[1:]:
                new_structure += model
            structure = new_structure

    if not utils.KEEP_NON_PROTEIN:
        structure = structure[np.isin(structure.res_name, list(utils.get_protein_components()))]

    if not utils.KEEP_HYDROGENS:
        structure = structure[~np.isin(structure.element, ["H"])]

    return structure


def get_structure_coords(structure: bts.AtomArray) -> pd.DataFrame:
    """
    Return the coordinates of an atom array as a dataframe.

    # TODO: shouldn't we just keep the coords as an np array?  why are we converting to a less performant format?
    # TODO: in fact, we could just keep the whole structure and thereby keep the elements along with it
    """
    coordinates = structure.coord
    elements = structure.element

    coordinates = pd.DataFrame(coordinates, columns=["x", "y", "z"])
    coordinates["element"] = elements

    return coordinates


def structure_to_lines(structure: bts.AtomArray) -> list[str]:
    """
    Convert an AtomArray into individual PDB lines
    """
    try:
        return [make_atom_line(atom.atom_name, atom.res_name, atom.res_id, atom.chain_id, idx, atom.coord, atom.element, beta=atom.b_factor) for idx, atom in enumerate(structure)]
    except AttributeError:
        return [make_atom_line(atom.atom_name, atom.res_name, atom.res_id, atom.chain_id, idx, atom.coord, atom.element) for idx, atom in enumerate(structure)]


def format_atom_name(atom_name: str) -> str:
    """
    Format an atom-name to occupy 4 characters the weird way PDBs do it.
    """
    if len(atom_name) == 1:
        return f" {atom_name}  "
    elif len(atom_name) == 2:
        return f" {atom_name} "
    elif len(atom_name) == 3:
        return f" {atom_name}"
    return atom_name[:4]


def make_atom_line(
    atom_name: str,
    res_name: str,
    res_num: int,
    chain_id: str,
    atom_index: int,
    coord: np.ndarray,
    element: str,
    beta: float = 0.0,
) -> str:
    """
    Make a single atom line for a single voxel
    """
    return f"ATOM  {atom_index:>5d} {format_atom_name(atom_name)} {res_name} {chain_id}{res_num:>4d}    {coord[0]:>8.3f}{coord[1]:>8.3f}{coord[2]:>8.3f}  1.00{beta:>6.2f}           {element}"
 

def make_voxel_lines(
    voxel_type: str,
    resnum: int,
    voxel_indices: set[int],
    surface_indices: set[int],
    voxel_grid_centers: np.ndarray,
) -> list[str]:
    """
    Make all atom lines for a given set of voxel indices
    """
    return [
        make_atom_line(VOXEL_TYPE_ATOM_MAP[voxel_type], voxel_type, resnum, VOXEL_TYPE_CHAIN_MAP[voxel_type], voxel_index, voxel_grid_centers[voxel_index], VOXEL_TYPE_ATOM_MAP[voxel_type], 50.0)
        if voxel_index in surface_indices
        else make_atom_line(VOXEL_TYPE_ATOM_MAP[voxel_type], voxel_type, resnum, VOXEL_TYPE_CHAIN_MAP[voxel_type], voxel_index, voxel_grid_centers[voxel_index], VOXEL_TYPE_ATOM_MAP[voxel_type], 0.0)
        for voxel_index in voxel_indices
    ]


def points_to_pdb(
    voxel_grid: VoxelGrid,
    hubs: dict[int, VoxelGroup],
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
                make_voxel_lines("HUB", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in hubs.items()
            ],
            *[
                make_voxel_lines("POR", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in pores.items()
            ],
            *[
                make_voxel_lines("POK", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in pockets.items()
            ],
            *[
                make_voxel_lines("CAV", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in cavities.items()
            ],
            *[
                make_voxel_lines("OCC", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in occluded.items()
            ],
        )
    )


def generate_rotation_translation_remarks(rotation: np.ndarray, translation: np.ndarray) -> list[str]:
    """
    Given a rotation matrix and translation vector, generate a string allowing them to be
    written as REMARKS to a PDB file.
    """
    return [
        f"REMARK TRANSLATION {str(translation)}",
        f"REMARK ROTATION [{str(rotation[0])}, {str(rotation[1])}, {str(rotation[2])}]",
    ]


def generate_volume_remarks(annotation: Annotation) -> list[str]:
    """
    Return REMARKS for a PDB file to hold all volume information from an annotation object.
    """
    hub_remarks = [
        f"REMARK HUB VOLUME {annotation.hub_volumes[i]} DIMENSIONS {annotation.hub_dimensions[i][0], annotation.hub_dimensions[i][1], annotation.hub_dimensions[i][2]}"
        for i in annotation.hub_volumes
    ]

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

    return hub_remarks + pore_remarks + cavity_remarks + pocket_remarks


def generate_resolution_remarks() -> list[str]:
    """
    Add a REMARK line with the resolution used fo rthis annotation
    """
    return [f"REMARK RESOLUTION {str(utils.VOXEL_SIZE)}"]


def merge_pdb_lines(lines: list[list[str]]) -> list[str]:
    """
    Merge together the various PDB-line elements into the final output file.
    """
    final_lines = lines[0]
    for lines in lines[1:]:
        final_lines.extend(lines)

    return final_lines
