"""
Functions for parsing, cleaning, and modifying PDBs.
"""

from pathlib import Path
from time import perf_counter

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


DEFAULT_ASSEMBLY_POLICY = "biological"
VALID_ASSEMBLY_POLICIES = (
    "biological",
    "asymmetric",
    "auto",
)

_LOAD_STRUCTURE_PARSE_STAGE = "load_structure_parse_decode"
_LOAD_STRUCTURE_ASSEMBLY_STAGE = "load_structure_assembly_expand"
_LOAD_STRUCTURE_FALLBACK_STAGE = "load_structure_fallback"
_IDENTITY_OPERATION_EXPRESSIONS = {"1", "(1)"}


def _accumulate_stage_timing(
    stage_timings: dict[str, float] | None,
    stage_name: str,
    elapsed_seconds: float,
) -> None:
    if stage_timings is None:
        return

    stage_timings[stage_name] = stage_timings.get(stage_name, 0.0) + float(
        elapsed_seconds
    )


def _safe_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _is_identity_operation_row(
    oper_list_category: dict[str, list[str]],
    row_index: int,
) -> bool:
    expected_fields = {
        "matrix[1][1]": 1.0,
        "matrix[1][2]": 0.0,
        "matrix[1][3]": 0.0,
        "matrix[2][1]": 0.0,
        "matrix[2][2]": 1.0,
        "matrix[2][3]": 0.0,
        "matrix[3][1]": 0.0,
        "matrix[3][2]": 0.0,
        "matrix[3][3]": 1.0,
        "vector[1]": 0.0,
        "vector[2]": 0.0,
        "vector[3]": 0.0,
    }

    for field_name, expected_value in expected_fields.items():
        field_values = oper_list_category.get(field_name)
        if field_values is None or row_index >= len(field_values):
            return False

        parsed_value = _safe_float(field_values[row_index])
        if parsed_value is None:
            return False

        if abs(parsed_value - expected_value) > 1e-6:
            return False

    return True


def _can_use_identity_assembly_shortcut_cif(file: pdbx.PDBxFile) -> bool:
    """
    Return True only when the default biological assembly is provably identity-only.

    This is intentionally strict: if required categories are missing or ambiguous,
    return False and fall back to explicit biological assembly expansion.
    """
    try:
        assembly_map = pdbx.list_assemblies(file)
    except (InvalidFileError, KeyError, ValueError, TypeError):
        return False

    if not assembly_map:
        return False

    selected_assembly_id = str(next(iter(assembly_map.keys())))

    try:
        assembly_gen = file.get_category("pdbx_struct_assembly_gen", expect_looped=True)
        oper_list = file.get_category("pdbx_struct_oper_list", expect_looped=True)
        atom_site = file.get_category("atom_site", expect_looped=True)
    except (KeyError, ValueError, TypeError):
        return False

    if assembly_gen is None or oper_list is None or atom_site is None:
        return False

    assembly_ids = assembly_gen.get("assembly_id")
    oper_expressions = assembly_gen.get("oper_expression")
    asym_id_lists = assembly_gen.get("asym_id_list")
    if (
        assembly_ids is None
        or oper_expressions is None
        or asym_id_lists is None
        or len(assembly_ids) != len(oper_expressions)
        or len(assembly_ids) != len(asym_id_lists)
    ):
        return False

    assembly_asym_ids: set[str] = set()
    selected_rows = 0
    for row_index, assembly_id in enumerate(assembly_ids):
        if str(assembly_id) != selected_assembly_id:
            continue

        selected_rows += 1
        expression = oper_expressions[row_index].replace(" ", "")
        if expression not in _IDENTITY_OPERATION_EXPRESSIONS:
            return False

        asym_ids = {
            token.strip()
            for token in asym_id_lists[row_index].split(",")
            if len(token.strip()) > 0
        }
        assembly_asym_ids |= asym_ids

    if selected_rows == 0 or len(assembly_asym_ids) == 0:
        return False

    oper_ids = oper_list.get("id")
    if oper_ids is None:
        return False

    identity_row_index = None
    for row_index, oper_id in enumerate(oper_ids):
        if str(oper_id).strip() == "1":
            identity_row_index = row_index
            break

    if identity_row_index is None:
        return False

    if not _is_identity_operation_row(oper_list, identity_row_index):
        return False

    asym_column = atom_site.get("label_asym_id")
    if asym_column is None:
        asym_column = atom_site.get("auth_asym_id")
    if asym_column is None:
        return False

    all_asym_ids = {
        str(asym_id).strip()
        for asym_id in asym_column
        if str(asym_id).strip() not in {"", "?", "."}
    }
    if len(all_asym_ids) == 0:
        return False

    return assembly_asym_ids == all_asym_ids


def load_structure(
    file_path: Path,
    assembly_policy: str = DEFAULT_ASSEMBLY_POLICY,
    stage_timings: dict[str, float] | None = None,
) -> bts.AtomArray:
    """
    Load various structure formats with explicit assembly policy controls.
    """
    normalized_policy = str(assembly_policy).strip().lower()
    if normalized_policy not in VALID_ASSEMBLY_POLICIES:
        raise ValueError(
            f"Unsupported assembly policy: {assembly_policy}. "
            f"Expected one of: {', '.join(VALID_ASSEMBLY_POLICIES)}"
        )

    suffix = file_path.suffix.lower()

    parse_start = perf_counter()
    parse_recorded = False
    assembly_start: float | None = None
    assembly_recorded = False

    try:
        if suffix == ".pdb":
            file = pdb.PDBFile.read(file_path)
            _accumulate_stage_timing(
                stage_timings,
                _LOAD_STRUCTURE_PARSE_STAGE,
                perf_counter() - parse_start,
            )
            parse_recorded = True

            if normalized_policy == "asymmetric":
                return pdb.get_structure(file, model=1)

            assembly_start = perf_counter()
            structure = pdb.get_assembly(file, model=1)
            _accumulate_stage_timing(
                stage_timings,
                _LOAD_STRUCTURE_ASSEMBLY_STAGE,
                perf_counter() - assembly_start,
            )
            assembly_recorded = True
            return structure

        if suffix in {".cif", ".mmcif"}:
            file = pdbx.PDBxFile.read(file_path)
            _accumulate_stage_timing(
                stage_timings,
                _LOAD_STRUCTURE_PARSE_STAGE,
                perf_counter() - parse_start,
            )
            parse_recorded = True

            if normalized_policy == "asymmetric":
                return pdbx.get_structure(file, model=1)

            if normalized_policy == "auto" and _can_use_identity_assembly_shortcut_cif(file):
                return pdbx.get_structure(file, model=1)

            assembly_start = perf_counter()
            structure = pdbx.get_assembly(file, model=1)
            _accumulate_stage_timing(
                stage_timings,
                _LOAD_STRUCTURE_ASSEMBLY_STAGE,
                perf_counter() - assembly_start,
            )
            assembly_recorded = True
            return structure

        if suffix == ".mmtf":
            file = mmtf.MMTFFile.read(file_path)
            _accumulate_stage_timing(
                stage_timings,
                _LOAD_STRUCTURE_PARSE_STAGE,
                perf_counter() - parse_start,
            )
            parse_recorded = True

            if normalized_policy == "asymmetric":
                return mmtf.get_structure(file, model=1)

            assembly_start = perf_counter()
            structure = mmtf.get_assembly(file, model=1)
            _accumulate_stage_timing(
                stage_timings,
                _LOAD_STRUCTURE_ASSEMBLY_STAGE,
                perf_counter() - assembly_start,
            )
            assembly_recorded = True
            return structure

        structure = biotite_load_structure(file_path)
        _accumulate_stage_timing(
            stage_timings,
            _LOAD_STRUCTURE_PARSE_STAGE,
            perf_counter() - parse_start,
        )
        return structure

    except (InvalidFileError, NotImplementedError):
        now = perf_counter()
        if not parse_recorded:
            _accumulate_stage_timing(
                stage_timings,
                _LOAD_STRUCTURE_PARSE_STAGE,
                now - parse_start,
            )
        elif assembly_start is not None and not assembly_recorded:
            _accumulate_stage_timing(
                stage_timings,
                _LOAD_STRUCTURE_ASSEMBLY_STAGE,
                now - assembly_start,
            )

        fallback_start = perf_counter()
        structure = biotite_load_structure(file_path)
        _accumulate_stage_timing(
            stage_timings,
            _LOAD_STRUCTURE_FALLBACK_STAGE,
            perf_counter() - fallback_start,
        )
        return structure


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
