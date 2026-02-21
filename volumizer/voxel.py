"""
Functions to manipulate and analyze voxels.
"""


from copy import deepcopy
import ctypes
from time import perf_counter

import numpy as np
import pandas as pd

from pyntcloud import PyntCloud
from pyntcloud.structures.voxelgrid import VoxelGrid

from volumizer import utils, native_backend
from volumizer.constants import OCCLUDED_DIMENSION_LIMIT, MIN_NUM_VOXELS
from volumizer.types import VoxelGroup
from volumizer.paths import C_CODE_DIR


if utils.using_performant():
    VOXEL_C_PATH = C_CODE_DIR / "voxel.so"
    VOXEL_C = ctypes.CDLL(str(VOXEL_C_PATH.absolute()))


NATIVE_COMPONENT_TYPE_CODE_MAP = {
    0: "occluded",
    1: "cavity",
    2: "pocket",
    3: "pore",
    4: "hub",
}


def _accumulate_stage_timing(
    stage_timings: dict[str, float] | None,
    stage_name: str,
    elapsed_seconds: float,
) -> None:
    """
    Accumulate elapsed time under a stage key when profiling is enabled.
    """
    if stage_timings is None:
        return

    stage_timings[stage_name] = stage_timings.get(stage_name, 0.0) + float(
        elapsed_seconds
    )


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
        "voxelgrid",
        size_x=utils.VOXEL_SIZE,
        size_y=utils.VOXEL_SIZE,
        size_z=utils.VOXEL_SIZE,
        regular_bounding_box=False,
    )
    return cloud, voxel_grid_id


def get_voxel_grid(cloud: PyntCloud, voxel_grid_id: str) -> VoxelGrid:
    """
    Generate an array representing a binary voxel grid where zero represents no protein
    atoms in that voxel (e.g. solvent) and a non-zero represents a protein atom in that voxel.
    """
    return cloud.structures[voxel_grid_id]


def get_protein_solvent_voxel_array(voxel_grid: VoxelGrid) -> np.ndarray:
    """
    Generate a 3D array of x[y[z]] as the voxel coordinate and the value of that array entry
    being zero for no protein and one when a protein atom is in that voxel.
    """
    return voxel_grid.get_feature_vector(mode="binary")


def get_protein_and_solvent_voxels(
    binary_voxel_array: np.ndarray,
    voxel_grid_dimensions: np.ndarray,
) -> tuple[VoxelGroup, VoxelGroup]:
    """
    From the voxel grid, return the indices of the protein containing voxels,
    and those not containing protein (e.g. solvent).
    """
    protein_voxels = np.nonzero(binary_voxel_array)
    solvent_voxels = np.nonzero(binary_voxel_array == 0)

    protein_voxel_indices = compute_voxel_indices(protein_voxels, voxel_grid_dimensions)
    solvent_voxel_indices = compute_voxel_indices(solvent_voxels, voxel_grid_dimensions)

    return (
        VoxelGroup(
            voxels=protein_voxels,
            indices=protein_voxel_indices,
            num_voxels=len(protein_voxel_indices),
            voxel_type="protein",
            volume=compute_voxel_group_volume(len(protein_voxel_indices)),
        ),
        VoxelGroup(
            voxels=solvent_voxels,
            indices=solvent_voxel_indices,
            num_voxels=len(solvent_voxel_indices),
            voxel_type="solvent",
            volume=compute_voxel_group_volume(len(solvent_voxel_indices)),
        ),
    )


def get_occluded_dimensions(
    query_voxel: tuple[int, int, int],
    occluding_x: list[int],
    occluding_y: list[int],
    occluding_z: list[int],
) -> list[int]:
    """
    Determine how many ordinal axes are occluded.

    We simply compare our query voxel with the coordinates of the `occluding` voxels,
    which are e.g. protein voxels.

    For performance reasons rather than looking across ALL occluding voxels, we look at planar arrays
    of coordinates.
    """
    occluded_dimensions = [0, 0, 0, 0, 0, 0]

    # here we assume that e.g. occluding_z[0] is the min value and occluding_z[-1] is the max
    #   if this becomes untrue we would need to explicitly take the min/max
    #   this is less performant for large structures
    if len(occluding_z) != 0:
        if occluding_z[0] < query_voxel[2]:
            occluded_dimensions[4] = 1
        if occluding_z[-1] > query_voxel[2]:
            occluded_dimensions[5] = 1
    if len(occluding_y) != 0:
        if occluding_y[0] < query_voxel[1]:
            occluded_dimensions[2] = 1
        if occluding_y[-1] > query_voxel[1]:
            occluded_dimensions[3] = 1
    if len(occluding_x) != 0:
        if occluding_x[0] < query_voxel[0]:
            occluded_dimensions[0] = 1
        if occluding_x[-1] > query_voxel[0]:
            occluded_dimensions[1] = 1

    return occluded_dimensions


def is_buried(occluded_dimensions: list[int]) -> bool:
    """
    If 5 or 6 dimensions are occluded, return True.
    If less than 4 dimensions are occluded, return False.
    If exactly 4 dimensions are occluded, return True if the two unoccluded dimensions are the same axis,
    False otherwise
    """

    if occluded_dimensions.count(1) > OCCLUDED_DIMENSION_LIMIT:
        return True
    elif occluded_dimensions.count(1) == OCCLUDED_DIMENSION_LIMIT:
        if occluded_dimensions[0] == 0 and occluded_dimensions[1] == 0:
            return True
        elif occluded_dimensions[2] == 0 and occluded_dimensions[3] == 0:
            return True
        elif occluded_dimensions[4] == 0 and occluded_dimensions[5] == 0:
            return True

    return False


def build_planar_voxel_coordinate_arrays(
    voxels: tuple[np.ndarray, ...], voxel_grid_dimensions: np.ndarray
) -> tuple[list[list[list[int]]], list[list[list[int]]], list[list[list[int]]],]:
    """
    For each two-dimension pair, construct a tuple of the min/max coordinates in the 3rd dimension.
    """
    z_array = [
        [[] for y in range(voxel_grid_dimensions[1])]
        for x in range(voxel_grid_dimensions[0])
    ]
    y_array = [
        [[] for z in range(voxel_grid_dimensions[2])]
        for x in range(voxel_grid_dimensions[0])
    ]
    x_array = [
        [[] for z in range(voxel_grid_dimensions[2])]
        for y in range(voxel_grid_dimensions[1])
    ]

    for i in range(voxels[0].size):
        z_array[voxels[0][i]][voxels[1][i]].append(voxels[2][i])
        y_array[voxels[0][i]][voxels[2][i]].append(voxels[1][i])
        x_array[voxels[1][i]][voxels[2][i]].append(voxels[0][i])

    return z_array, y_array, x_array


def _build_exposed_and_buried_voxel_groups(
    exposed_voxels: tuple[np.ndarray, np.ndarray, np.ndarray],
    buried_voxels: tuple[np.ndarray, np.ndarray, np.ndarray],
    voxel_grid_dimensions: np.ndarray,
) -> tuple[VoxelGroup, VoxelGroup]:
    """
    Build typed voxel groups for exposed/buried solvent subsets.
    """
    exposed_arrays = (
        np.asarray(exposed_voxels[0]),
        np.asarray(exposed_voxels[1]),
        np.asarray(exposed_voxels[2]),
    )
    buried_arrays = (
        np.asarray(buried_voxels[0]),
        np.asarray(buried_voxels[1]),
        np.asarray(buried_voxels[2]),
    )

    buried_voxel_indices = compute_voxel_indices(buried_arrays, voxel_grid_dimensions)
    exposed_voxel_indices = compute_voxel_indices(exposed_arrays, voxel_grid_dimensions)

    return (
        VoxelGroup(
            voxels=exposed_arrays,
            indices=exposed_voxel_indices,
            num_voxels=len(exposed_voxel_indices),
            voxel_type="exposed",
            volume=compute_voxel_group_volume(len(exposed_voxel_indices)),
        ),
        VoxelGroup(
            voxels=buried_arrays,
            indices=buried_voxel_indices,
            num_voxels=len(buried_voxel_indices),
            voxel_type="buried",
            volume=compute_voxel_group_volume(len(buried_voxel_indices)),
        ),
    )


def get_exposed_and_buried_voxels_python(
    solvent_voxels: VoxelGroup,
    protein_voxels: VoxelGroup,
    voxel_grid_dimensions: np.ndarray,
) -> tuple[VoxelGroup, VoxelGroup]:
    """
    Use simple geometric heuristics to determine if a group of solvent voxels is buried or exposed.
    """
    buried_voxels = ([], [], [])
    exposed_voxels = ([], [], [])

    z_array, y_array, x_array = build_planar_voxel_coordinate_arrays(
        protein_voxels.voxels, voxel_grid_dimensions
    )

    for i in range(solvent_voxels.voxels[0].size):
        query_voxel = (
            solvent_voxels.voxels[0][i],
            solvent_voxels.voxels[1][i],
            solvent_voxels.voxels[2][i],
        )

        occluded_dimensions = get_occluded_dimensions(
            query_voxel,
            x_array[query_voxel[1]][query_voxel[2]],
            y_array[query_voxel[0]][query_voxel[2]],
            z_array[query_voxel[0]][query_voxel[1]],
        )

        if is_buried(occluded_dimensions):
            buried_voxels[0].append(query_voxel[0])
            buried_voxels[1].append(query_voxel[1])
            buried_voxels[2].append(query_voxel[2])
        else:
            exposed_voxels[0].append(query_voxel[0])
            exposed_voxels[1].append(query_voxel[1])
            exposed_voxels[2].append(query_voxel[2])

    return _build_exposed_and_buried_voxel_groups(
        exposed_voxels,
        buried_voxels,
        voxel_grid_dimensions,
    )

def get_exposed_and_buried_voxels_native(
    solvent_voxels: VoxelGroup,
    protein_voxels: VoxelGroup,
    voxel_grid_dimensions: np.ndarray,
) -> tuple[VoxelGroup, VoxelGroup]:
    """
    Native backend version of exposed/buried solvent split.
    """
    native_module = native_backend.get_native_module_for_mode("native")
    if native_module is None:
        raise RuntimeError(
            "Native backend requested but `volumizer_native` is not importable."
        )

    if not hasattr(native_module, "get_exposed_and_buried_voxel_indices"):
        raise RuntimeError(
            "Native backend does not provide `get_exposed_and_buried_voxel_indices`."
        )

    solvent_array = _voxel_tuple_to_native_array(solvent_voxels.voxels)
    protein_array = _voxel_tuple_to_native_array(protein_voxels.voxels)
    native_output = native_module.get_exposed_and_buried_voxel_indices(
        solvent_array,
        protein_array,
        np.asarray(voxel_grid_dimensions, dtype=np.int32),
    )
    if not isinstance(native_output, dict):
        raise RuntimeError(
            "Unexpected native exposed/buried output, expected dict with index arrays."
        )

    exposed_indices = np.asarray(native_output["exposed_indices"], dtype=np.int64)
    buried_indices = np.asarray(native_output["buried_indices"], dtype=np.int64)

    exposed_voxels = (
        solvent_voxels.voxels[0][exposed_indices],
        solvent_voxels.voxels[1][exposed_indices],
        solvent_voxels.voxels[2][exposed_indices],
    )
    buried_voxels = (
        solvent_voxels.voxels[0][buried_indices],
        solvent_voxels.voxels[1][buried_indices],
        solvent_voxels.voxels[2][buried_indices],
    )

    return _build_exposed_and_buried_voxel_groups(
        exposed_voxels, buried_voxels, voxel_grid_dimensions
    )


def get_exposed_and_buried_voxels(
    solvent_voxels: VoxelGroup,
    protein_voxels: VoxelGroup,
    voxel_grid_dimensions: np.ndarray,
    backend: str | None = None,
) -> tuple[VoxelGroup, VoxelGroup]:
    """
    Use simple geometric heuristics to determine if a group of solvent voxels is buried or exposed.
    """
    requested_backend = backend if backend is not None else utils.get_active_backend()
    if requested_backend == "native":
        native_module = native_backend.get_native_module_for_mode("native")
        if native_module is not None and hasattr(
            native_module, "get_exposed_and_buried_voxel_indices"
        ):
            return get_exposed_and_buried_voxels_native(
                solvent_voxels, protein_voxels, voxel_grid_dimensions
            )

    return get_exposed_and_buried_voxels_python(
        solvent_voxels, protein_voxels, voxel_grid_dimensions
    )


def get_neighbor_voxels_python(
    query_voxels: tuple[np.ndarray, np.ndarray, np.ndarray],
    reference_voxels: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return the voxels from a query group that neighbor a reference group.
    """
    neighbor_indices = []
    for query_index in range(len(query_voxels[0])):
        query_voxel = get_single_voxel(query_voxels, query_index)
        for reference_index in range(len(reference_voxels[0])):
            reference_voxel = get_single_voxel(reference_voxels, reference_index)
            if is_neighbor_voxel(query_voxel, reference_voxel):
                neighbor_indices.append(query_index)
                break

    return (
        np.array([query_voxels[0][i] for i in neighbor_indices]),
        np.array([query_voxels[1][i] for i in neighbor_indices]),
        np.array([query_voxels[2][i] for i in neighbor_indices]),
    )


def get_neighbor_voxels_c(
    query_voxels: tuple[np.ndarray, np.ndarray, np.ndarray],
    reference_voxels: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return the voxels from a query group that neighbor a reference group.
    """
    # make the voxels compatible with C
    num_query = len(query_voxels[0])
    query_x = (ctypes.c_int * num_query)(*query_voxels[0])
    query_y = (ctypes.c_int * num_query)(*query_voxels[1])
    query_z = (ctypes.c_int * num_query)(*query_voxels[2])

    num_reference = len(reference_voxels[0])
    reference_x = (ctypes.c_int * num_reference)(*reference_voxels[0])
    reference_y = (ctypes.c_int * num_reference)(*reference_voxels[1])
    reference_z = (ctypes.c_int * num_reference)(*reference_voxels[2])

    # setup return and call the function
    initial_indices = (ctypes.c_int * num_query)(*([-1] * num_query))
    VOXEL_C.get_neighbor_voxels.restype = ctypes.POINTER(ctypes.c_int * num_query)
    first_shell_indices_c = VOXEL_C.get_neighbor_voxels(
        ctypes.byref(query_x),
        ctypes.byref(query_y),
        ctypes.byref(query_z),
        ctypes.byref(reference_x),
        ctypes.byref(reference_y),
        ctypes.byref(reference_z),
        num_query,
        num_reference,
        ctypes.byref(initial_indices),
    )
    first_shell_indices = [i for i in first_shell_indices_c.contents if i != -1]

    return (
        np.array([query_voxels[0][i] for i in first_shell_indices]),
        np.array([query_voxels[1][i] for i in first_shell_indices]),
        np.array([query_voxels[2][i] for i in first_shell_indices]),
    )


def get_neighbor_voxels_native(
    query_voxels: tuple[np.ndarray, np.ndarray, np.ndarray],
    reference_voxels: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return the voxels from a query group that neighbor a reference group using native backend.
    """
    native_module = native_backend.get_native_module_for_mode("native")
    if native_module is None:
        raise RuntimeError(
            "Native backend requested but `volumizer_native` is not importable."
        )

    query_array = _voxel_tuple_to_native_array(query_voxels)
    reference_array = _voxel_tuple_to_native_array(reference_voxels)
    neighbor_indices = np.asarray(
        native_module.get_neighbor_voxel_indices(query_array, reference_array),
        dtype=np.int64,
    )
    if neighbor_indices.size == 0:
        return np.array([]), np.array([]), np.array([])

    return (
        query_voxels[0][neighbor_indices],
        query_voxels[1][neighbor_indices],
        query_voxels[2][neighbor_indices],
    )


def get_first_shell_exposed_voxels(
    exposed_voxels: VoxelGroup,
    buried_voxels: VoxelGroup,
    voxel_grid: VoxelGrid,
    performant: bool = utils.using_performant(),
    backend: str | None = None,
) -> VoxelGroup:
    """
    Subset exposed voxels into only those neighboring one or more buried voxels.
    """
    requested_backend = backend if backend is not None else utils.get_active_backend()

    if requested_backend == "native":
        native_module = native_backend.get_native_module_for_mode("native")
        if native_module is None:
            raise RuntimeError(
                "Native backend requested but `volumizer_native` is not importable."
            )

        if hasattr(native_module, "get_first_shell_exposed_indices"):
            exposed_array = _voxel_tuple_to_native_array(exposed_voxels.voxels)
            buried_array = _voxel_tuple_to_native_array(buried_voxels.voxels)
            grid_dimensions = np.asarray(voxel_grid.x_y_z, dtype=np.int32)
            first_shell_neighbor_indices = np.asarray(
                native_module.get_first_shell_exposed_indices(
                    exposed_array,
                    buried_array,
                    grid_dimensions,
                ),
                dtype=np.int64,
            )
            first_shell_arrays = (
                exposed_voxels.voxels[0][first_shell_neighbor_indices],
                exposed_voxels.voxels[1][first_shell_neighbor_indices],
                exposed_voxels.voxels[2][first_shell_neighbor_indices],
            )
        else:
            first_shell_voxels = get_neighbor_voxels_native(
                exposed_voxels.voxels, buried_voxels.voxels
            )
            first_shell_arrays = (
                np.asarray(first_shell_voxels[0]),
                np.asarray(first_shell_voxels[1]),
                np.asarray(first_shell_voxels[2]),
            )
    elif performant:
        first_shell_voxels = get_neighbor_voxels_c(
            exposed_voxels.voxels, buried_voxels.voxels
        )
        first_shell_arrays = (
            np.asarray(first_shell_voxels[0]),
            np.asarray(first_shell_voxels[1]),
            np.asarray(first_shell_voxels[2]),
        )
    else:
        first_shell_voxels = get_neighbor_voxels_python(
            exposed_voxels.voxels, buried_voxels.voxels
        )
        first_shell_arrays = (
            np.asarray(first_shell_voxels[0]),
            np.asarray(first_shell_voxels[1]),
            np.asarray(first_shell_voxels[2]),
        )
    first_shell_indices = compute_voxel_indices(first_shell_arrays, voxel_grid.x_y_z)

    return VoxelGroup(
        voxels=first_shell_arrays,
        indices=first_shell_indices,
        num_voxels=len(first_shell_indices),
        voxel_type="exposed",
        volume=compute_voxel_group_volume(len(first_shell_indices)),
    )


def is_neighbor_voxel(
    voxel_one: tuple[np.int64, ...], voxel_two: tuple[np.int64, ...]
) -> bool:
    """
    Given two voxels return True if they are ordinal neighbors.
    """
    # below approach is ~2x faster than using a list comprehension, presumably because of early exit
    sum_differences = 0
    for dimension in range(3):
        difference = abs(voxel_one[dimension] - voxel_two[dimension])

        if difference > 1:
            return False

        sum_differences += difference

    # a voxel is an ordinal neighbor when the sum of the absolute differences in axes indices is 1
    if sum_differences == 1:
        return True

    return False


def get_single_voxel(
    voxels: tuple[np.ndarray, ...], index: int
) -> tuple[np.int64, ...]:
    """
    Given a set of voxels return just one voxel at index.
    """
    return (
        voxels[0][index],
        voxels[1][index],
        voxels[2][index],
    )


def _voxel_tuple_to_native_array(voxels: tuple[np.ndarray, ...]) -> np.ndarray:
    """
    Convert `(x,y,z)` voxel arrays into contiguous int32 `(N,3)` for native kernels.
    """
    num_voxels = len(voxels[0])
    native_voxels = np.empty((num_voxels, 3), dtype=np.int32)
    if num_voxels == 0:
        return native_voxels

    native_voxels[:, 0] = voxels[0]
    native_voxels[:, 1] = voxels[1]
    native_voxels[:, 2] = voxels[2]
    return native_voxels


def _compute_voxel_indices_array(
    voxels: tuple[np.ndarray, ...], grid_dimensions: np.ndarray
) -> np.ndarray:
    """
    Given a 3D array of voxels, compute flat 1D voxel-grid indices.
    """
    if len(voxels[0]) == 0:
        return np.array([], dtype=np.int64)

    grid_dimension_1_2 = int(grid_dimensions[1]) * int(grid_dimensions[2])
    grid_dimension_2 = int(grid_dimensions[2])
    x_coords = np.asarray(voxels[0], dtype=np.int64)
    y_coords = np.asarray(voxels[1], dtype=np.int64)
    z_coords = np.asarray(voxels[2], dtype=np.int64)

    return x_coords * grid_dimension_1_2 + y_coords * grid_dimension_2 + z_coords


def compute_voxel_indices(
    voxels: tuple[np.ndarray, ...], grid_dimensions: np.ndarray
) -> set[int]:
    """
    Given a 3D array of voxels, compute the 1D set of their indices.
    """
    return set(_compute_voxel_indices_array(voxels, grid_dimensions).tolist())


def compute_voxel_group_volume(num_voxels: int) -> float:
    """
    Return the volume of a number of voxels
    """
    return num_voxels * utils.VOXEL_VOLUME


def breadth_first_search_python(
    voxels: tuple[np.ndarray, ...], searchable_indices: set[int]
) -> set[int]:
    """
    Given a set of voxels and list of possible indices to add,
    add indices for all ordinal neighbors iteratively until no more such neighbors exist.
    """
    searchable_indices = deepcopy(searchable_indices)
    start_index = searchable_indices.pop()
    queue_indices = set([start_index])
    neighbor_indices = set([start_index])

    while len(queue_indices) > 0:
        current_index = queue_indices.pop()
        current_voxel = get_single_voxel(voxels, current_index)
        for searched_index in searchable_indices:
            searched_voxel = get_single_voxel(voxels, searched_index)
            if is_neighbor_voxel(current_voxel, searched_voxel):
                queue_indices.add(searched_index)
                neighbor_indices.add(searched_index)
        searchable_indices -= neighbor_indices

    return neighbor_indices


def breadth_first_search_c(
    voxels: tuple[np.ndarray, ...], searchable_indices: set[int]
) -> set[int]:
    """
    Given a set of voxels and list of possible indices to add,
    add indices for all ordinal neighbors iteratively until no more such neighbors exist.
    """
    c_searchable_indices = list(deepcopy(searchable_indices))

    # make the voxels compatible with C
    num_voxels = len(voxels[0])
    voxels_x = (ctypes.c_int * num_voxels)(*voxels[0])
    voxels_y = (ctypes.c_int * num_voxels)(*voxels[1])
    voxels_z = (ctypes.c_int * num_voxels)(*voxels[2])

    # generate other inputs needed for C function
    c_searchable_indices = (ctypes.c_int * len(c_searchable_indices))(
        *c_searchable_indices
    )
    initial_indices = (ctypes.c_int * num_voxels)(*([-1] * num_voxels))

    # run the breadth-first search
    VOXEL_C.breadth_first_search.restype = ctypes.POINTER(ctypes.c_int * num_voxels)
    neighbor_indices_c = VOXEL_C.breadth_first_search(
        ctypes.byref(voxels_x),
        ctypes.byref(voxels_y),
        ctypes.byref(voxels_z),
        ctypes.byref(c_searchable_indices),
        len(searchable_indices),
        ctypes.byref(initial_indices),
    )

    return set([i for i in neighbor_indices_c.contents if i != -1])


def breadth_first_search_native(
    voxels: tuple[np.ndarray, ...], searchable_indices: set[int]
) -> set[int]:
    """
    Native backend version of breadth-first search over voxel indices.
    """
    if len(searchable_indices) == 0:
        return set()

    native_module = native_backend.get_native_module_for_mode("native")
    if native_module is None:
        raise RuntimeError(
            "Native backend requested but `volumizer_native` is not importable."
        )

    voxel_array = _voxel_tuple_to_native_array(voxels)
    searchable_index_array = np.fromiter(
        searchable_indices,
        dtype=np.int32,
        count=len(searchable_indices),
    )
    neighbor_indices = native_module.bfs_component_indices(
        voxel_array, searchable_index_array
    )
    return set(np.asarray(neighbor_indices, dtype=np.int64).tolist())


def breadth_first_search(
    voxels: tuple[np.ndarray, ...],
    searchable_indices: set[int],
    performant: bool = utils.using_performant(),
    backend: str | None = None,
) -> set[int]:
    """
    Given a set of voxels and list of possible indices to add,
    add indices for all ordinal neighbors iteratively until no more such neighbors exist.
    """
    requested_backend = backend if backend is not None else utils.get_active_backend()
    if requested_backend == "native":
        return breadth_first_search_native(voxels, searchable_indices)

    if performant:
        return breadth_first_search_c(voxels, searchable_indices)

    return breadth_first_search_python(voxels, searchable_indices)


def is_edge_voxel(voxel: np.ndarray, voxel_grid_dimensions: np.ndarray) -> bool:
    """
    Is this voxel at the edge of our grid on any dimension?
    """
    if 0 in voxel:
        return True
    elif voxel[0] == voxel_grid_dimensions[0] - 1:
        return True
    elif voxel[1] == voxel_grid_dimensions[1] - 1:
        return True
    elif voxel[2] == voxel_grid_dimensions[2] - 1:
        return True

    return False


def get_agglomerated_type(
    query_indices: set[int],
    buried_voxels: tuple[np.ndarray, ...],
    exposed_voxels: tuple[np.ndarray, ...],
    voxel_grid_dimensions: np.ndarray,
    backend: str | None = None,
) -> tuple[set[int], str]:
    """
    Find "surface" voxels, being buried voxel in direct contact with an exposed voxel. Four possibilites:

    a) there are no "surface" voxels -> this is a cavity
    b) all "surface" voxels can be agglomerated (BFS) into a single group -> this is a pocket
    c) the "surface" voxels can be agglomerated (BFS) into exactly two groups -> this is a pore
    d) the "surface" voxels cannot be agglomerated (BFS) into less than 3 groups -> this is a hub
    """
    direct_surface_indices = set()
    for query_index in query_indices:
        query_voxel = get_single_voxel(buried_voxels, query_index)
        if is_edge_voxel(query_voxel, voxel_grid_dimensions):
            direct_surface_indices.add(query_index)
            continue
        for exposed_voxel_index in range(len(exposed_voxels[0])):
            exposed_voxel = get_single_voxel(exposed_voxels, exposed_voxel_index)
            if is_neighbor_voxel(query_voxel, exposed_voxel):
                direct_surface_indices.add(query_index)
                break

    # if there are no surface contacts, this must be a cavity
    if len(direct_surface_indices) == 0:
        return direct_surface_indices, "cavity"

    # add all voxels that are neighbours to the direct surface voxels
    neighbor_surface_indices = set()
    for query_index in query_indices - direct_surface_indices:
        query_voxel = get_single_voxel(buried_voxels, query_index)
        for surface_index in direct_surface_indices:
            surface_voxel = get_single_voxel(buried_voxels, surface_index)
            if is_neighbor_voxel(surface_voxel, query_voxel):
                neighbor_surface_indices.add(query_index)
                break

    # we pass all surface indices and the buried voxels for BFS
    # If there is more than one distinct surface (e.g. a pore) then
    #   BFS will terminate after agglomerating just one of the surfaces
    #   and `single_surface_indices` will be less than `surface_indices`
    # NOTE: we have to union the direct and neighbor surfaces, otherwise small discritization
    #   on the surface would look like a distinct surface
    surface_indices = direct_surface_indices.union(neighbor_surface_indices)
    single_surface_indices = breadth_first_search(
        buried_voxels, surface_indices, backend=backend
    )
    if len(single_surface_indices) < len(surface_indices):
        # run the agglomeration one more time on what's left
        remaining_surface_indices = surface_indices - single_surface_indices
        single_surface_indices = breadth_first_search(
            buried_voxels, remaining_surface_indices, backend=backend
        )

        # if there are stil indices left, there were more than 2 surfaces, hence a hub
        if len(single_surface_indices) < len(remaining_surface_indices):
            return direct_surface_indices.union(neighbor_surface_indices), "hub"

        # othewise, exactly 2 and a pore
        return direct_surface_indices.union(neighbor_surface_indices), "pore"

    # if the first BFS went to completion and consumed all surfaces, there was only a single surface, hence pocket
    return direct_surface_indices.union(neighbor_surface_indices), "pocket"


def _get_voxel_group_center_from_indices(
    voxel_indices: np.ndarray,
    voxel_grid: VoxelGrid,
) -> np.ndarray:
    """
    Compute the center of geometry for flat voxel-grid indices.
    """
    if voxel_indices.size == 0:
        return np.array([0.0, 0.0, 0.0])

    voxel_coords = np.asarray(voxel_grid.voxel_centers[voxel_indices], dtype=float)
    return np.mean(voxel_coords, axis=0)


def _get_voxel_group_axial_lengths_from_indices(
    voxel_indices: np.ndarray,
    voxel_grid: VoxelGrid,
) -> list[float]:
    """
    Compute principal-axis lengths for flat voxel-grid indices.
    """
    if voxel_indices.size < 3:
        return [0, 0, 0]

    voxel_coords = np.asarray(voxel_grid.voxel_centers[voxel_indices], dtype=np.float64)
    centered_coords = voxel_coords - np.mean(voxel_coords, axis=0)

    covariance = centered_coords.T @ centered_coords
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    axis_order = np.argsort(eigenvalues)[::-1]
    principal_axes = eigenvectors[:, axis_order]

    aligned_coords = centered_coords @ principal_axes
    axial_lengths = np.round(np.ptp(aligned_coords, axis=0), 3)
    return sorted(axial_lengths.tolist(), reverse=True)


def get_voxel_group_center(
    voxel_indices: set[int], voxel_grid: VoxelGrid
) -> np.ndarray:
    """
    Compute the center of geometry for the voxel group.
    """
    voxel_index_array = np.fromiter(
        voxel_indices,
        dtype=np.int64,
        count=len(voxel_indices),
    )
    return _get_voxel_group_center_from_indices(voxel_index_array, voxel_grid)


def get_voxel_group_axial_lengths(
    voxel_indices: set[int], voxel_grid: VoxelGrid
) -> list[float]:
    """
    Align the voxel group to it's principal axes, then compute the maximum length along each axis.
    """
    voxel_index_array = np.fromiter(
        voxel_indices,
        dtype=np.int64,
        count=len(voxel_indices),
    )
    return _get_voxel_group_axial_lengths_from_indices(voxel_index_array, voxel_grid)


def _extract_component_voxels(
    buried_voxels: VoxelGroup, component_indices: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert indices into a component voxel tuple.
    """
    component_indices = np.asarray(component_indices, dtype=np.int64)
    if component_indices.size == 0:
        return (
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
        )

    return (
        buried_voxels.voxels[0][component_indices],
        buried_voxels.voxels[1][component_indices],
        buried_voxels.voxels[2][component_indices],
    )


def _build_voxel_group_from_component_indices(
    buried_voxels: VoxelGroup,
    voxel_grid: VoxelGrid,
    component_indices: np.ndarray,
    surface_component_indices: np.ndarray,
    agglomerated_type: str,
    buried_flat_indices: np.ndarray | None = None,
) -> VoxelGroup:
    """
    Construct a VoxelGroup from native component-index outputs.
    """
    component_indices = np.asarray(component_indices, dtype=np.int64)
    surface_component_indices = np.asarray(surface_component_indices, dtype=np.int64)

    volume_voxels = _extract_component_voxels(buried_voxels, component_indices)

    if buried_flat_indices is None:
        surface_voxels = _extract_component_voxels(
            buried_voxels, surface_component_indices
        )
        volume_index_values = compute_voxel_indices(volume_voxels, voxel_grid.x_y_z)
        surface_indices = compute_voxel_indices(surface_voxels, voxel_grid.x_y_z)
        voxel_index_array = np.fromiter(
            volume_index_values,
            dtype=np.int64,
            count=len(volume_index_values),
        )
        volume_indices: set[int] | np.ndarray = volume_index_values
    else:
        voxel_index_array = np.asarray(
            buried_flat_indices[component_indices],
            dtype=np.int64,
        )
        surface_index_array = np.asarray(
            buried_flat_indices[surface_component_indices],
            dtype=np.int64,
        )
        # Keep component indices as a NumPy array to avoid large set materialization overhead.
        volume_indices = voxel_index_array
        surface_indices = surface_index_array

    if voxel_index_array.size == 0:
        center = np.array([0.0, 0.0, 0.0])
        axial_lengths = [0.0, 0.0, 0.0]
    else:
        center = _get_voxel_group_center_from_indices(voxel_index_array, voxel_grid)
        axial_lengths = _get_voxel_group_axial_lengths_from_indices(
            voxel_index_array,
            voxel_grid,
        )

    num_voxels = int(len(voxel_index_array))
    return VoxelGroup(
        voxels=volume_voxels,
        indices=volume_indices,
        surface_indices=surface_indices,
        num_voxels=num_voxels,
        voxel_type=agglomerated_type,
        volume=compute_voxel_group_volume(num_voxels),
        center=center,
        axial_lengths=axial_lengths,
    )


def classify_buried_components_native(
    buried_voxels: VoxelGroup,
    exposed_voxels: VoxelGroup,
    voxel_grid: VoxelGrid,
    stage_timings: dict[str, float] | None = None,
) -> tuple[
    dict[int, VoxelGroup],
    dict[int, VoxelGroup],
    dict[int, VoxelGroup],
    dict[int, VoxelGroup],
    dict[int, VoxelGroup],
]:
    """
    Use native classifier output (if exposed by the native backend) and map into VoxelGroups.
    """
    native_module = native_backend.get_native_module_for_mode("native")
    if native_module is None or not hasattr(
        native_module, "classify_buried_components"
    ):
        raise RuntimeError(
            "Native backend does not provide `classify_buried_components`."
        )

    total_start = perf_counter() if stage_timings is not None else 0.0
    kernel_seconds = 0.0
    mapping_seconds = 0.0

    buried_array = _voxel_tuple_to_native_array(buried_voxels.voxels)
    exposed_array = _voxel_tuple_to_native_array(exposed_voxels.voxels)
    grid_dimensions = np.asarray(voxel_grid.x_y_z, dtype=np.int32)

    kernel_start = perf_counter() if stage_timings is not None else 0.0
    native_output = native_module.classify_buried_components(
        buried_array,
        exposed_array,
        grid_dimensions,
        int(MIN_NUM_VOXELS),
        float(utils.VOXEL_SIZE),
    )
    if stage_timings is not None:
        kernel_seconds = perf_counter() - kernel_start
        _accumulate_stage_timing(
            stage_timings,
            "classify_components_native_kernel",
            kernel_seconds,
        )

    if not isinstance(native_output, dict):
        raise RuntimeError(
            "Unexpected native classifier output, expected dict with flattened arrays."
        )

    if stage_timings is not None:
        kernel_stage_timings = native_output.get("kernel_stage_timings_seconds")
        if isinstance(kernel_stage_timings, dict):
            for substage_name, substage_seconds in kernel_stage_timings.items():
                try:
                    parsed_seconds = float(substage_seconds)
                except (TypeError, ValueError):
                    continue
                if parsed_seconds <= 0.0:
                    continue
                _accumulate_stage_timing(
                    stage_timings,
                    f"classify_components_native_kernel_{substage_name}",
                    parsed_seconds,
                )

    mapping_start = perf_counter() if stage_timings is not None else 0.0

    component_type_codes = np.asarray(
        native_output["component_type_codes"], dtype=np.int64
    )
    component_offsets = np.asarray(native_output["component_offsets"], dtype=np.int64)
    voxel_indices_flat = np.asarray(
        native_output["component_voxel_indices_flat"], dtype=np.int64
    )
    surface_indices_flat = np.asarray(
        native_output["component_surface_indices_flat"], dtype=np.int64
    )
    surface_offsets = np.asarray(
        native_output.get("surface_offsets", component_offsets), dtype=np.int64
    )

    if (
        len(component_offsets) != len(component_type_codes) + 1
        or len(surface_offsets) != len(component_type_codes) + 1
    ):
        raise RuntimeError(
            "Native classifier output has inconsistent component offsets."
        )

    buried_flat_indices = _compute_voxel_indices_array(
        buried_voxels.voxels, voxel_grid.x_y_z
    )

    hubs, pores, pockets, cavities, occluded = {}, {}, {}, {}, {}
    counters = {"hub": 0, "pore": 0, "pocket": 0, "cavity": 0, "occluded": 0}
    bucket_map = {
        "hub": hubs,
        "pore": pores,
        "pocket": pockets,
        "cavity": cavities,
        "occluded": occluded,
    }

    for component_id, type_code in enumerate(component_type_codes):
        agglomerated_type = NATIVE_COMPONENT_TYPE_CODE_MAP.get(int(type_code))
        if agglomerated_type is None:
            raise RuntimeError(f"Unknown native component type code: {type_code}")

        volume_slice = slice(
            component_offsets[component_id],
            component_offsets[component_id + 1],
        )
        surface_slice = slice(
            surface_offsets[component_id],
            surface_offsets[component_id + 1],
        )
        component_indices = voxel_indices_flat[volume_slice]
        surface_component_indices = surface_indices_flat[surface_slice]

        voxel_group = _build_voxel_group_from_component_indices(
            buried_voxels,
            voxel_grid,
            component_indices,
            surface_component_indices,
            agglomerated_type,
            buried_flat_indices=buried_flat_indices,
        )

        bucket = bucket_map[agglomerated_type]
        type_index = counters[agglomerated_type]
        bucket[type_index] = voxel_group
        counters[agglomerated_type] += 1

    if stage_timings is not None:
        mapping_seconds = perf_counter() - mapping_start
        _accumulate_stage_timing(
            stage_timings,
            "classify_components_native_mapping",
            mapping_seconds,
        )

        other_seconds = (perf_counter() - total_start) - kernel_seconds - mapping_seconds
        if other_seconds > 0:
            _accumulate_stage_timing(
                stage_timings,
                "classify_components_native_other",
                other_seconds,
            )

    return hubs, pores, pockets, cavities, occluded


def get_pores_pockets_cavities_occluded(
    buried_voxels: VoxelGroup,
    exposed_voxels: VoxelGroup,
    voxel_grid: VoxelGrid,
    backend: str | None = None,
    stage_timings: dict[str, float] | None = None,
) -> tuple[
    dict[int, VoxelGroup],
    dict[int, VoxelGroup],
    dict[int, VoxelGroup],
    dict[int, VoxelGroup],
    dict[int, VoxelGroup],
]:
    """
    Agglomerate buried solvent voxels into hubs, pores, pockets, cavities, and simply occluded.
    """
    requested_backend = backend if backend is not None else utils.get_active_backend()
    if requested_backend == "native":
        native_module = native_backend.get_native_module_for_mode("native")
        if native_module is not None and hasattr(
            native_module, "classify_buried_components"
        ):
            if stage_timings is None:
                return classify_buried_components_native(
                    buried_voxels,
                    exposed_voxels,
                    voxel_grid,
                )
            return classify_buried_components_native(
                buried_voxels,
                exposed_voxels,
                voxel_grid,
                stage_timings=stage_timings,
            )

    python_start = perf_counter() if stage_timings is not None else 0.0

    buried_indices = set(range(buried_voxels.voxels[0].size))
    agglomerated_indices = set()

    hubs, pores, pockets, cavities, occluded = {}, {}, {}, {}, {}
    hub_id, pore_id, pocket_id, cavity_id, occluded_id = 0, 0, 0, 0, 0

    while len(agglomerated_indices) < len(buried_indices):
        remaining_indices = buried_indices - agglomerated_indices

        # perform BFS over the remaining indices
        agglomerable_indices = breadth_first_search(
            buried_voxels.voxels, remaining_indices, backend=requested_backend
        )
        # iterate our counter of finished indices
        agglomerated_indices = agglomerated_indices.union(agglomerable_indices)

        # if too small, don't assign direct surface indices and assign type "occluded"
        if len(agglomerable_indices) <= MIN_NUM_VOXELS:
            direct_surface_indices = set()
            agglomerated_type = "occluded"
        else:
            # identify what these agglomerated voxels are
            direct_surface_indices, agglomerated_type = get_agglomerated_type(
                agglomerable_indices,
                buried_voxels.voxels,
                exposed_voxels.voxels,
                voxel_grid.x_y_z,
                backend=requested_backend,
            )

        # get the surface voxels for use in getting their voxel-grid indices
        surface_voxels = (
            np.array(
                [buried_voxels.voxels[0][index] for index in direct_surface_indices]
            ),
            np.array(
                [buried_voxels.voxels[1][index] for index in direct_surface_indices]
            ),
            np.array(
                [buried_voxels.voxels[2][index] for index in direct_surface_indices]
            ),
        )

        # get the voxels
        volume_voxels = (
            np.array(
                [buried_voxels.voxels[0][index] for index in agglomerable_indices]
            ),
            np.array(
                [buried_voxels.voxels[1][index] for index in agglomerable_indices]
            ),
            np.array(
                [buried_voxels.voxels[2][index] for index in agglomerable_indices]
            ),
        )
        # get the voxel indices
        volume_indices = compute_voxel_indices(volume_voxels, voxel_grid.x_y_z)
        # create the voxelgroup
        voxel_group = VoxelGroup(
            voxels=volume_voxels,
            indices=volume_indices,
            surface_indices=compute_voxel_indices(surface_voxels, voxel_grid.x_y_z),
            num_voxels=len(volume_indices),
            voxel_type=agglomerated_type,
            volume=compute_voxel_group_volume(len(volume_indices)),
            center=get_voxel_group_center(volume_indices, voxel_grid),
            axial_lengths=get_voxel_group_axial_lengths(volume_indices, voxel_grid),
        )

        # add the voxelgroup depending on type and increment the type counter
        if agglomerated_type == "hub":
            hubs[hub_id] = voxel_group
            hub_id += 1
        elif agglomerated_type == "pore":
            pores[pore_id] = voxel_group
            pore_id += 1
        elif agglomerated_type == "pocket":
            pockets[pocket_id] = voxel_group
            pocket_id += 1
        elif agglomerated_type == "cavity":
            cavities[cavity_id] = voxel_group
            cavity_id += 1
        elif agglomerated_type == "occluded":
            occluded[occluded_id] = voxel_group
            occluded_id += 1

    if stage_timings is not None:
        _accumulate_stage_timing(
            stage_timings,
            "classify_components_python",
            perf_counter() - python_start,
        )

    return hubs, pores, pockets, cavities, occluded


