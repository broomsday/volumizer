import pytest
import numpy as np
import ctypes

from pore.voxel import get_single_voxel, is_neighbor_voxel, breadth_first_search
from pore.paths import C_CODE_DIR


voxel_compute_path = C_CODE_DIR / "voxel_compute.so"
voxel_compute = ctypes.CDLL(str(voxel_compute_path.absolute()))


@pytest.mark.parametrize(
    "voxels, index, voxel",
    [
        (
            (np.array([0, 9, 12]), np.array([1, 2, 7]), np.array([21, 0, 0])),
            0,
            (0, 1, 21),
        ),
        (
            (np.array([0, 9, 12]), np.array([1, 2, 7]), np.array([21, 0, 0])),
            1,
            (9, 2, 0),
        ),
    ],
)
def test_get_single_voxel(voxels, index, voxel):
    assert get_single_voxel(voxels, index) == voxel


@pytest.mark.parametrize(
    "voxels, index, voxel",
    [
        (
            (np.array([0, 9, 12]), np.array([1, 2, 7]), np.array([21, 0, 0])),
            0,
            [0, 1, 21],
        ),
        (
            (np.array([0, 9, 12]), np.array([1, 2, 7]), np.array([21, 0, 0])),
            1,
            [9, 2, 0],
        ),
    ],
)
def test_get_single_voxel_c(voxels, index, voxel):
    num_voxels = len(voxels[0])
    voxels_x = (ctypes.c_int * num_voxels)(*voxels[0])
    voxels_y = (ctypes.c_int * num_voxels)(*voxels[1])
    voxels_z = (ctypes.c_int * num_voxels)(*voxels[2])

    initial_single_voxel_indices = (ctypes.c_int * 3)(*([-1] * 3))

    voxel_compute.get_single_voxel.restype = ctypes.POINTER(ctypes.c_int * 3)
    voxel_indices_c = voxel_compute.get_single_voxel(
        ctypes.byref(voxels_x),
        ctypes.byref(voxels_y),
        ctypes.byref(voxels_z),
        index,
        ctypes.byref(initial_single_voxel_indices),
    )

    assert [i for i in voxel_indices_c.contents] == voxel


@pytest.mark.parametrize(
    "voxel_one, voxel_two, diagonal, is_neighbor",
    [
        (
            (1, 7, 9),
            (2, 7, 9),
            True,
            True,
        ),
        (
            (1, 7, 9),
            (2, 8, 9),
            True,
            True,
        ),
        (
            (1, 7, 9),
            (2, 8, 8),
            True,
            False,
        ),
        (
            (1, 7, 9),
            (2, 7, 9),
            False,
            True,
        ),
        (
            (1, 7, 9),
            (2, 8, 9),
            False,
            False,
        ),
    ],
)
def test_is_neighbor_voxel(voxel_one, voxel_two, diagonal, is_neighbor):
    assert is_neighbor_voxel(voxel_one, voxel_two, diagonal) == is_neighbor


@pytest.mark.parametrize(
    "voxel_one, voxel_two, diagonal, is_neighbor",
    [
        (
            (1, 7, 9),
            (2, 7, 9),
            True,
            True,
        ),
        (
            (1, 7, 9),
            (2, 8, 9),
            True,
            True,
        ),
        (
            (1, 7, 9),
            (2, 8, 8),
            True,
            False,
        ),
        (
            (1, 7, 9),
            (2, 7, 9),
            False,
            True,
        ),
        (
            (1, 7, 9),
            (2, 8, 9),
            False,
            False,
        ),
    ],
)
def test_is_neighbor_voxel_c(voxel_one, voxel_two, diagonal, is_neighbor):
    voxel_one = (ctypes.c_int * 3)(*voxel_one)
    voxel_two = (ctypes.c_int * 3)(*voxel_two)
    diagonal_c = int(diagonal)

    is_neighbor_c = voxel_compute.is_neighbor_voxel(
        ctypes.byref(voxel_one),
        ctypes.byref(voxel_two),
        diagonal_c,
    )

    assert bool(is_neighbor_c) == is_neighbor


@pytest.mark.parametrize(
    "voxels, searchable_indices, neighbor_indices",
    [
        (
            (
                np.array([0, 1, 1, 1, 1, 2]),
                np.array([0, 0, 1, 1, 2, 2]),
                np.array([0, 0, 0, 1, 0, 2]),
            ),
            set([0, 1, 2, 3, 4, 5]),
            set([0, 1, 2, 3, 4]),
        ),
        (
            (
                np.array([0, 1, 1, 1, 1, 2]),
                np.array([0, 0, 1, 1, 2, 2]),
                np.array([0, 0, 0, 1, 0, 2]),
            ),
            set([0, 1, 4, 5]),
            set([0, 1]),
        ),
    ],
)
def test_breadth_first_search(voxels, searchable_indices, neighbor_indices):
    assert breadth_first_search(voxels, searchable_indices) == neighbor_indices


@pytest.mark.parametrize(
    "voxels, searchable_indices, neighbor_indices",
    [
        (
            (
                np.array([0, 1, 1, 1, 1, 2]),
                np.array([0, 0, 1, 1, 2, 2]),
                np.array([0, 0, 0, 1, 0, 2]),
            ),
            [0, 1, 2, 3, 4, 5],
            [0, 1, 2, 3, 4],
        ),
        (
            (
                np.array([0, 1, 1, 1, 1, 2]),
                np.array([0, 0, 1, 1, 2, 2]),
                np.array([0, 0, 0, 1, 0, 2]),
            ),
            [0, 1, 4, 5],
            [0, 1],
        ),
    ],
)
def test_breadth_first_search_c(voxels, searchable_indices, neighbor_indices):
    num_voxels = len(voxels[0])
    voxels_x = (ctypes.c_int * num_voxels)(*voxels[0])
    voxels_y = (ctypes.c_int * num_voxels)(*voxels[1])
    voxels_z = (ctypes.c_int * num_voxels)(*voxels[2])

    searchable_indices = (ctypes.c_int * len(searchable_indices))(*searchable_indices)

    initial_indices = (ctypes.c_int * num_voxels)(*([-1] * num_voxels))

    voxel_compute.breadth_first_search.restype = ctypes.POINTER(ctypes.c_int * num_voxels)
    neighbor_indices_c = voxel_compute.breadth_first_search(
        ctypes.byref(voxels_x),
        ctypes.byref(voxels_y),
        ctypes.byref(voxels_z),
        ctypes.byref(searchable_indices),
        len(searchable_indices),
        ctypes.byref(initial_indices),
    )

    assert [i for i in neighbor_indices_c.contents if i != -1] == neighbor_indices
