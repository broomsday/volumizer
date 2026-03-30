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


def _generate_chain_id_pool() -> list[str]:
    """
    Return a deterministic pool of unique chain IDs up to Biotite's `U4` width.

    The biological-assembly path can create hundreds or thousands of chain
    copies. The older `A-Z`, `a-z`, `0-9`, `AA-ZZ` pool exhausted on larger
    assemblies such as 8B12/2BBV. Stay within 4 characters so assignments fit
    the fixed-width chain-id dtype Biotite uses for structure arrays.
    """
    import itertools
    import string

    alphabet = string.ascii_uppercase + string.ascii_lowercase + string.digits
    pool: list[str] = []
    for width in range(1, 5):
        for chars in itertools.product(alphabet, repeat=width):
            pool.append("".join(chars))
    return pool


def _deduplicate_assembly_chain_ids(structure: bts.AtomArray) -> bts.AtomArray:
    """
    Assign unique chain IDs to symmetry copies produced by assembly expansion.

    Biotite's get_assembly() reuses the original chain IDs for each copy.
    This function detects repeated chain ID blocks and renames them so every
    physical chain has a distinct ID.
    """
    chain_ids = structure.chain_id
    n = len(chain_ids)
    if n == 0:
        return structure

    # Identify contiguous chain blocks: [(chain_id, start, end), ...]
    blocks: list[tuple[str, int, int]] = []
    block_start = 0
    for i in range(1, n):
        if chain_ids[i] != chain_ids[i - 1]:
            blocks.append((chain_ids[block_start], block_start, i))
            block_start = i
    blocks.append((chain_ids[block_start], block_start, n))

    # Detect copy boundaries: a chain ID reappearing after we've moved past it
    # within the current copy means a new copy started.
    seen_in_copy: set[str] = set()
    copy_index = 0
    block_copy: list[int] = []

    for block_id, _, _ in blocks:
        if block_id in seen_in_copy:
            copy_index += 1
            seen_in_copy.clear()
        seen_in_copy.add(block_id)
        block_copy.append(copy_index)

    num_copies = copy_index + 1
    if num_copies <= 1:
        # Chain-block detection failed — happens when the asymmetric unit has
        # only one chain so all copies share a single contiguous chain ID.
        # Fall back to splitting by residue-ID resets (res_id jumping backwards
        # by more than a small gap marks a copy boundary).
        if len(blocks) == 1:
            block_id, block_start, block_end = blocks[0]
            res_ids = structure.res_id[block_start:block_end]
            drops = np.where(np.diff(res_ids) < -10)[0]
            if len(drops) > 0:
                # Re-build blocks and block_copy from the detected boundaries
                split_points = [block_start] + [block_start + int(d) + 1 for d in drops] + [block_end]
                blocks = [
                    (block_id, split_points[i], split_points[i + 1])
                    for i in range(len(split_points) - 1)
                ]
                block_copy = list(range(len(blocks)))
                num_copies = len(blocks)

        if num_copies <= 1:
            return structure

    # Build mapping: (original_chain_id, copy_index) -> new_chain_id
    original_ids_ordered: list[str] = []
    seen_originals: set[str] = set()
    for block_id, _, _ in blocks:
        if block_id not in seen_originals:
            seen_originals.add(block_id)
            original_ids_ordered.append(block_id)

    pool = _generate_chain_id_pool()
    # Reserve the original IDs for copy 0, assign new ones for copies 1+
    used: set[str] = set(original_ids_ordered)
    pool_iter = iter(cid for cid in pool if cid not in used)

    rename_map: dict[tuple[str, int], str] = {}
    for orig_id in original_ids_ordered:
        rename_map[(orig_id, 0)] = orig_id
        for ci in range(1, num_copies):
            try:
                rename_map[(orig_id, ci)] = next(pool_iter)
            except StopIteration as error:
                raise RuntimeError(
                    "Exhausted generated chain IDs while deduplicating biological "
                    f"assembly copies: original_chains={len(original_ids_ordered)}, "
                    f"copies={num_copies}. Increase chain-ID generation capacity."
                ) from error

    # Apply renaming
    new_chain_ids = chain_ids.copy()
    for (block_id, start, end), ci in zip(blocks, block_copy):
        new_id = rename_map[(block_id, ci)]
        new_chain_ids[start:end] = new_id

    structure = structure.copy()
    structure.chain_id = new_chain_ids
    return structure


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
            return _deduplicate_assembly_chain_ids(structure)

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
            return _deduplicate_assembly_chain_ids(structure)

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
            return _deduplicate_assembly_chain_ids(structure)

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


_VOLUME_RES_NAMES = frozenset(VOXEL_TYPE_CHAIN_MAP.keys())


def _add_entity_table(
    pdbx_file: pdbx.PDBxFile,
    data_block: str = "structure",
) -> None:
    """
    Write an ``_entity`` category so that Mol* can distinguish polymer
    chains from volume pseudo-atoms.  Without this table Mol* must guess
    entity types, and structures starting with modified residues (e.g.
    FME) are misclassified as non-polymer.
    """
    try:
        atom_site = pdbx_file.get_category("atom_site", block=data_block)
    except KeyError:
        return

    if "label_entity_id" not in atom_site or "label_comp_id" not in atom_site:
        return

    entity_ids = np.asarray(atom_site["label_entity_id"])
    comp_ids = np.asarray(atom_site["label_comp_id"])

    unique_eids = sorted(set(entity_ids), key=lambda x: (int(x) if x.isdigit() else 0, x))
    entity_types = []

    for eid in unique_eids:
        names = set(comp_ids[entity_ids == eid])
        # If any residue is NOT a volume pseudo-atom → polymer
        if names - _VOLUME_RES_NAMES:
            entity_types.append("polymer")
        else:
            entity_types.append("non-polymer")

    pdbx_file.set_category(
        "entity",
        {"id": np.array(unique_eids), "type": np.array(entity_types)},
        block=data_block,
    )


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
        _add_entity_table(pdbx_file, data_block="structure")
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


def get_structure_residue_count(structure: bts.AtomArray) -> int:
    """
    Return the number of residues represented in a structure.
    """
    return int(bts.get_residue_count(structure))


def compute_sse_fractions(
    structure: bts.AtomArray,
) -> dict[str, float | None]:
    """
    Compute secondary structure element fractions using biotite P-SEA.

    Returns dict with keys ``frac_alpha``, ``frac_beta``, ``frac_coil``,
    each a float in [0, 1] or None if annotation is not possible.
    """
    try:
        sse = bts.annotate_sse(structure)
    except Exception:
        return {"frac_alpha": None, "frac_beta": None, "frac_coil": None}

    assigned = sse[sse != ""]
    if len(assigned) == 0:
        return {"frac_alpha": None, "frac_beta": None, "frac_coil": None}

    total = len(assigned)
    return {
        "frac_alpha": float(np.sum(assigned == "a")) / total,
        "frac_beta": float(np.sum(assigned == "b")) / total,
        "frac_coil": float(np.sum(assigned == "c")) / total,
    }


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


def _build_volume_chain_map(
    existing_chain_ids: set[str],
) -> dict[str, str]:
    """
    Return a voxel-type → chain-ID mapping that avoids *existing_chain_ids*.

    Falls back to the default ``VOXEL_TYPE_CHAIN_MAP`` when there is no
    collision.  When a collision exists the next available single-character
    ID is chosen from ``A-Z``, ``a-z``, ``0-9``.
    """
    candidates = (
        [chr(c) for c in range(ord("A"), ord("Z") + 1)]
        + [chr(c) for c in range(ord("a"), ord("z") + 1)]
        + [chr(c) for c in range(ord("0"), ord("9") + 1)]
    )

    used = set(existing_chain_ids)
    chain_map: dict[str, str] = {}

    for vtype, default_id in VOXEL_TYPE_CHAIN_MAP.items():
        if default_id not in used:
            chain_map[vtype] = default_id
        else:
            for cand in candidates:
                if cand not in used:
                    chain_map[vtype] = cand
                    break
            else:
                chain_map[vtype] = default_id  # last resort
        used.add(chain_map[vtype])

    return chain_map


def volumes_to_structure(
    voxel_grid: VoxelGrid,
    hubs: dict[int, VoxelGroup],
    pores: dict[int, VoxelGroup],
    pockets: dict[int, VoxelGroup],
    cavities: dict[int, VoxelGroup],
    occluded: dict[int, VoxelGroup],
    existing_chain_ids: set[str] | None = None,
) -> bts.AtomArray:
    """
    Convert the voxels of all volumes into a set atoms in a biotite AtomArray.

    *existing_chain_ids*, when provided, is the set of chain IDs already
    used by the protein structure.  Volume pseudo-atoms will be assigned
    chain IDs that do not collide with these.
    """
    chain_map = (
        _build_volume_chain_map(existing_chain_ids)
        if existing_chain_ids
        else VOXEL_TYPE_CHAIN_MAP
    )

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
        all_chain_ids[cursor:next_cursor] = chain_map[voxel_type]
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
