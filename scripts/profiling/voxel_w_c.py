from copy import deepcopy
import ctypes

import numpy as np

from pore.constants import DIAGONAL_NEIGHBORS
from pore.paths import C_CODE_DIR


voxel_compute_path = C_CODE_DIR / "voxel.so"
voxel_compute = ctypes.CDLL(str(voxel_compute_path.absolute()))
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
    c_searchable_indices = list(deepcopy(searchable_indices))

    # make the voxels compatible with C
    num_voxels = len(voxels[0])
    voxels_x = (ctypes.c_int * num_voxels)(*voxels[0])
    voxels_y = (ctypes.c_int * num_voxels)(*voxels[1])
    voxels_z = (ctypes.c_int * num_voxels)(*voxels[2])

    # generate other inputs needed for C function
    c_searchable_indices = (ctypes.c_int * len(c_searchable_indices))(*c_searchable_indices)
    initial_indices = (ctypes.c_int * num_voxels)(*([-1] * num_voxels))

    # run the breadth-first search
    voxel_compute.breadth_first_search.restype = ctypes.POINTER(ctypes.c_int * num_voxels)
    neighbor_indices_c = voxel_compute.breadth_first_search(
        ctypes.byref(voxels_x),
        ctypes.byref(voxels_y),
        ctypes.byref(voxels_z),
        ctypes.byref(c_searchable_indices),
        len(searchable_indices),
        ctypes.byref(initial_indices),
    )

    return set([i for i in neighbor_indices_c.contents if i != -1])
