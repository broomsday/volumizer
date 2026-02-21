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


def _index_values_to_array(index_values: set[int] | np.ndarray) -> np.ndarray:
    """
    Normalize voxel index containers into a flat int64 numpy array.
    """
    if isinstance(index_values, set):
        return np.fromiter(
            index_values,
            dtype=np.int64,
            count=len(index_values),
        )

    return np.asarray(index_values, dtype=np.int64).reshape(-1)


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
    voxel_index_array = _index_values_to_array(voxel_indices)
    surface_index_array = _index_values_to_array(surface_indices)

    num_voxels = int(voxel_index_array.size)
    volume_structure = bts.AtomArray(num_voxels)

    is_surface = (
        np.isin(voxel_index_array, surface_index_array)
        if surface_index_array.size > 0
        else np.zeros(num_voxels, dtype=bool)
    )

    volume_structure.coord = voxel_grid_centers[voxel_index_array]
    volume_structure.atom_id = voxel_index_array
    volume_structure.b_factor = np.where(is_surface, 50.0, 0.0)

    volume_structure.atom_name = np.full(num_voxels, VOXEL_TYPE_ATOM_MAP[voxel_type])
    volume_structure.res_name = np.full(num_voxels, voxel_type)
    volume_structure.res_id = np.full(num_voxels, voxel_group_index, dtype=np.int64)
    volume_structure.chain_id = np.full(num_voxels, VOXEL_TYPE_CHAIN_MAP[voxel_type])
    volume_structure.element = np.full(num_voxels, VOXEL_TYPE_ELEMENT_MAP[voxel_type])

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
    grouped_volumes = [
        ("HUB", hubs),
        ("POR", pores),
        ("POK", pockets),
        ("CAV", cavities),
        ("OCC", occluded),
    ]

    volume_chunks: list[tuple[str, int, np.ndarray, np.ndarray]] = []
    total_voxels = 0

    for voxel_type, groups in grouped_volumes:
        for group_index, voxel_group in groups.items():
            voxel_index_array = _index_values_to_array(voxel_group.indices)
            if voxel_index_array.size == 0:
                continue

            surface_index_array = _index_values_to_array(voxel_group.surface_indices)
            volume_chunks.append(
                (
                    voxel_type,
                    int(group_index),
                    voxel_index_array,
                    surface_index_array,
                )
            )
            total_voxels += int(voxel_index_array.size)

    if total_voxels == 0:
        return bts.AtomArray(0)

    all_voxel_indices = np.empty(total_voxels, dtype=np.int64)
    all_b_factors = np.zeros(total_voxels, dtype=np.float64)
    all_atom_names = np.empty(total_voxels, dtype=object)
    all_res_names = np.empty(total_voxels, dtype=object)
    all_res_ids = np.empty(total_voxels, dtype=np.int64)
    all_chain_ids = np.empty(total_voxels, dtype=object)
    all_elements = np.empty(total_voxels, dtype=object)

    cursor = 0
    for (
        voxel_type,
        voxel_group_index,
        voxel_index_array,
        surface_index_array,
    ) in volume_chunks:
        chunk_size = int(voxel_index_array.size)
        next_cursor = cursor + chunk_size

        all_voxel_indices[cursor:next_cursor] = voxel_index_array
        if surface_index_array.size > 0:
            all_b_factors[cursor:next_cursor] = np.where(
                np.isin(voxel_index_array, surface_index_array),
                50.0,
                0.0,
            )

        all_atom_names[cursor:next_cursor] = VOXEL_TYPE_ATOM_MAP[voxel_type]
        all_res_names[cursor:next_cursor] = voxel_type
        all_res_ids[cursor:next_cursor] = voxel_group_index
        all_chain_ids[cursor:next_cursor] = VOXEL_TYPE_CHAIN_MAP[voxel_type]
        all_elements[cursor:next_cursor] = VOXEL_TYPE_ELEMENT_MAP[voxel_type]

        cursor = next_cursor

    volume_structure = bts.AtomArray(total_voxels)
    volume_structure.coord = voxel_grid.voxel_centers[all_voxel_indices]
    volume_structure.atom_id = all_voxel_indices
    volume_structure.b_factor = all_b_factors
    volume_structure.atom_name = all_atom_names
    volume_structure.res_name = all_res_names
    volume_structure.res_id = all_res_ids
    volume_structure.chain_id = all_chain_ids
    volume_structure.element = all_elements

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
