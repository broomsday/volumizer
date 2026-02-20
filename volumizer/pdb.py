"""
Functions for parsing, cleaning, and modifying PDBs.
"""

from pathlib import Path

import biotite.structure as bts
from biotite.structure.io import load_structure as biotite_load_structure
from biotite.structure.io import mmtf, pdbx, pdb
from biotite import InvalidFileError
from pyntcloud.structures.voxelgrid import VoxelGrid
import pandas as pd
import numpy as np

from volumizer.constants import (
    VOXEL_TYPE_CHAIN_MAP,
    VOXEL_TYPE_ATOM_MAP,
    VOXEL_TYPE_ELEMENT_MAP,
)
from volumizer.types import VoxelGroup
from volumizer import utils


def load_structure(file_path: Path) -> bts.AtomArray:
    """
    Load various structure formats.
    """
    try:
        if file_path.suffix == ".pdb":
            file = pdb.PDBFile.read(file_path)
            assembly = pdb.get_assembly(file, model=1)
        elif file_path.suffix == ".cif":
            file = pdbx.PDBxFile.read(file_path)
            assembly = pdbx.get_assembly(file, model=1)
        elif file_path.suffix == ".mmtf":
            file = mmtf.MMTFFile.read(file_path)
            assembly = mmtf.get_assembly(file, model=1)
        else:
            return biotite_load_structure(file_path)
        return assembly
    except (InvalidFileError, NotImplementedError):
        return biotite_load_structure(file_path)


def save_structure(structure: bts.AtomArray, output: Path | str) -> None:
    """
    Save a structure to PDB or CIF based on file suffix.
    """
    output_path = Path(output)
    suffix = output_path.suffix.lower()

    if suffix == ".pdb":
        pdb_file = pdb.PDBFile()
        pdb_file.set_structure(structure)
        pdb_file.write(output_path)
        return

    if suffix in {".cif", ".mmcif"}:
        pdbx_file = pdbx.PDBxFile()
        pdbx.set_structure(pdbx_file, structure, data_block="structure")
        pdbx_file.write(output_path)
        return

    raise ValueError(f"Unsupported structure output format: {suffix}")


def save_pdb_lines(pdb_lines: list[str], output: Path | str) -> None:
    """
    Save individual pdb lines to a file.
    """
    if isinstance(output, (Path, str)):
        with open(output, mode="w", encoding="utf-8") as out_file:
            out_file.writelines("%s\n" % line for line in pdb_lines)
    else:
        output.writelines("%s\n" % line for line in pdb_lines)


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
        structure = structure[bts.filter_amino_acids(structure)]

    if not utils.KEEP_HYDROGENS:
        structure = structure[~np.isin(structure.element, ["H"])]

    return structure


def get_structure_coords(structure: bts.AtomArray) -> pd.DataFrame:
    """
    Return the coordinates of an atom array as a dataframe.

    This is done because the PyntCloud package takes a dataframe as input to generate the cloud.
    """
    coordinates = structure.coord
    elements = structure.element

    coordinates = pd.DataFrame(coordinates, columns=["x", "y", "z"])
    coordinates["element"] = elements

    return coordinates


def structure_to_pdb(structure: bts.AtomArray) -> list[str]:
    """
    Convert an AtomArray into individual PDB lines
    """
    try:
        return [
            make_atom_line(
                atom.atom_name,
                atom.res_name,
                atom.res_id,
                atom.chain_id,
                idx,
                atom.coord,
                atom.element,
                beta=atom.b_factor,
            )
            for idx, atom in enumerate(structure)
        ]
    except AttributeError:
        return [
            make_atom_line(
                atom.atom_name,
                atom.res_name,
                atom.res_id,
                atom.chain_id,
                idx,
                atom.coord,
                atom.element,
            )
            for idx, atom in enumerate(structure)
        ]


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
    Make a single atom line in a PDB
    """
    return f"ATOM  {atom_index:>5d} {format_atom_name(atom_name)} {res_name} {chain_id}{res_num:>4d}    {coord[0]:>8.3f}{coord[1]:>8.3f}{coord[2]:>8.3f}  1.00{beta:>6.2f}           {element}"


def volume_to_structure(
    voxel_type: str,
    voxel_group_index: int,
    voxel_indices: set[int] | np.ndarray,
    surface_indices: set[int] | np.ndarray,
    voxel_grid_centers: np.ndarray,
) -> bts.AtomArray:
    """
    Convert one volume into a set of atoms in a biotite AtomArray.
    """
    if isinstance(voxel_indices, set):
        voxel_index_array = np.fromiter(
            voxel_indices,
            dtype=np.int64,
            count=len(voxel_indices),
        )
    else:
        voxel_index_array = np.asarray(voxel_indices, dtype=np.int64).reshape(-1)

    num_voxels = int(voxel_index_array.size)
    volume_structure = bts.AtomArray(num_voxels)

    if isinstance(surface_indices, set):
        surface_index_set = surface_indices
    else:
        surface_index_set = set(np.asarray(surface_indices, dtype=np.int64).tolist())

    for i, voxel_index in enumerate(voxel_index_array):
        b_factor = 50.0 if int(voxel_index) in surface_index_set else 0.0
        volume_structure[i] = bts.Atom(
            voxel_grid_centers[int(voxel_index)],
            atom_id=int(voxel_index),
            b_factor=b_factor,
        )

    volume_structure.atom_name = [VOXEL_TYPE_ATOM_MAP[voxel_type]] * num_voxels
    volume_structure.res_name = [voxel_type] * num_voxels
    volume_structure.res_id = [voxel_group_index] * num_voxels
    volume_structure.chain_id = [VOXEL_TYPE_CHAIN_MAP[voxel_type]] * num_voxels
    volume_structure.element = [VOXEL_TYPE_ELEMENT_MAP[voxel_type]] * num_voxels

    return volume_structure



def volumes_to_structure(
    voxel_grid: VoxelGrid,
    hubs: dict[int, VoxelGroup],
    pores: dict[int, VoxelGroup],
    pockets: dict[int, VoxelGroup],
    cavities: dict[int, VoxelGroup],
    occluded: dict[int, VoxelGroup],
) -> bts.AtomArray:
    """
    Convert the voxels of all volumes into a set atoms in a biotite AtomArray.
    """
    volume_structure = bts.AtomArray(0)

    for i, voxel_group in hubs.items():
        volume_structure += volume_to_structure(
            "HUB",
            i,
            voxel_group.indices,
            voxel_group.surface_indices,
            voxel_grid.voxel_centers,
        )
    for i, voxel_group in pores.items():
        volume_structure += volume_to_structure(
            "POR",
            i,
            voxel_group.indices,
            voxel_group.surface_indices,
            voxel_grid.voxel_centers,
        )
    for i, voxel_group in pockets.items():
        volume_structure += volume_to_structure(
            "POK",
            i,
            voxel_group.indices,
            voxel_group.surface_indices,
            voxel_grid.voxel_centers,
        )
    for i, voxel_group in cavities.items():
        volume_structure += volume_to_structure(
            "CAV",
            i,
            voxel_group.indices,
            voxel_group.surface_indices,
            voxel_grid.voxel_centers,
        )
    for i, voxel_group in occluded.items():
        volume_structure += volume_to_structure(
            "OCC",
            i,
            voxel_group.indices,
            voxel_group.surface_indices,
            voxel_grid.voxel_centers,
        )

    return volume_structure


def make_volumized_pdb_lines(
    structures: list[bts.AtomArray],
    deliminator: str = "END",
) -> list[str]:
    """
    Convenience function to produce the final set of PDB lines for the annotated PDB.
    """
    pdb_lines = []
    for structure in structures:
        pdb_lines.extend(structure_to_pdb(structure))
        pdb_lines.append(deliminator)

    return pdb_lines


def coordinates_to_structure(
    coords: np.ndarray, res_num: int = 0, chain_id: str = "A", element: str = "C"
) -> bts.AtomArray:
    """
    Convert a set of 3D coordinates into a biotite structure.

    By default these are carbon atoms.
    """
    # TODO: this needs testing alone and with pipeline
    #   left here for future improvement/use
    coordinate_structure = bts.AtomArray(len(coords))
    for i, coord in enumerate(coords):
        coordinate_structure[i] = bts.Atom(coord, atom_id=i, b_factor=0.0)

    coordinate_structure.atom_name = [element] * len(coords)
    coordinate_structure.res_name = ["CRD"] * len(coords)
    coordinate_structure.res_num = [res_num] * len(coords)
    coordinate_structure.chain_id = [chain_id] * len(coords)
    coordinate_structure.element = [element] * len(coords)

    return coordinate_structure
