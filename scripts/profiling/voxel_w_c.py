from copy import deepcopy
import ctypes

import numpy as np

from pore.constants import DIAGONAL_NEIGHBORS
from pore.paths import C_CODE_DIR


voxel_compute_path = C_CODE_DIR / "voxel_compute.so"
voxel_compute = ctypes.CDLL(voxel_compute_path.absolute())
# voxel_compute.get_single_voxel.restype = ctypes.c_void_p
voxel_compute.get_single_voxel.restype = ctypes.POINTER(ctypes.c_int * 3)


class VoxelIndices(ctypes.Structure):
    _fields_ = [("x", ctypes.c_int), ("y", ctypes.c_int), ("z", ctypes.c_int)]


def is_neighbor_voxel(voxel_one, voxel_two, diagonal_neighbors: bool = DIAGONAL_NEIGHBORS) -> bool:
    """
    Given two voxels return True if they are ordinal neighbors.
    """
    differences = []
    for dimension in range(3):
        differences.append(abs(voxel_one[dimension] - voxel_two[dimension]))

    # a voxel is an ordinal neighbor when the sum of the absolute differences in axes indices is 1
    if sum(differences) == 1:
        return True
    # alternatively a voxel is a diagonal neighbour when the max distance on any axis is 1
    # here we only count being diagonal on a plane, not in all 3 dimensions
    #   which would be a longer distance between voxels
    elif (diagonal_neighbors) and (sum(differences) == 2) and (max(differences) == 1):
        return True

    return False


def get_single_voxel(voxels: tuple[np.ndarray, ...], index: int) -> tuple[np.int64, ...]:
    """
    Given a set of voxels return just one voxel at index.
    """
    return (
        voxels[0][index],
        voxels[1][index],
        voxels[2][index],
    )


def breadth_first_search(voxels: tuple[np.ndarray, ...], searchable_indices: set[int]) -> set[int]:
    """
    Given a set of voxels and list of possible indices to add,
    add indices for all ordinal neighbors iteratively until no more such neighbors exist.
    """
    searchable_indices = deepcopy(searchable_indices)
    start_index = searchable_indices.pop()
    queue_indices = set([start_index])
    neighbor_indices = set([start_index])

    # make the voxels compatible with C
    num_voxels = len(voxels[0])
    voxels_x = (ctypes.c_int * num_voxels)(*voxels[0])
    voxels_y = (ctypes.c_int * num_voxels)(*voxels[1])
    voxels_z = (ctypes.c_int * num_voxels)(*voxels[2])

    # TODO: make the searchable indices and resulting neighbor indices compatible with C

    # NOTE: below block is just test-code, have not written the C version of breadth_first_search, just working on sub-parts
    voxel_compute.breadth_first_search.restype = ctypes.POINTER(ctypes.c_int * num_voxels)
    initial_indices = (ctypes.c_int * num_voxels)(*([-1] * num_voxels))
    test_indices_c = voxel_compute.breadth_first_search(ctypes.byref(initial_indices))
    test_indices = [i for i in test_indices_c.contents]

    initial_single_voxel_indices = (ctypes.c_int * 3)(*([-1] * 3))
    while len(queue_indices) > 0:
        current_index = queue_indices.pop()
        current_voxel_c = voxel_compute.get_single_voxel(
            ctypes.byref(voxels_x),
            ctypes.byref(voxels_y),
            ctypes.byref(voxels_z),
            current_index,
            ctypes.byref(initial_single_voxel_indices),
        )
        current_voxel = [i for i in current_voxel_c.contents]

        for searched_index in searchable_indices:
            searched_voxel_c = voxel_compute.get_single_voxel(
                ctypes.byref(voxels_x),
                ctypes.byref(voxels_y),
                ctypes.byref(voxels_z),
                searched_index,
                ctypes.byref(initial_single_voxel_indices),
            )
            searched_voxel = [i for i in searched_voxel_c.contents]

            if is_neighbor_voxel(current_voxel, searched_voxel):
                queue_indices.add(searched_index)
                neighbor_indices.add(searched_index)
        searchable_indices -= neighbor_indices

    return neighbor_indices
