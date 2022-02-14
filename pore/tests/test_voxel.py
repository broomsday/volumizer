import pytest
import numpy as np
import ctypes

from pore.voxel import (
    get_single_voxel,
    is_neighbor_voxel,
    breadth_first_search_python,
    breadth_first_search_c,
    get_neighbor_voxels_python,
    get_neighbor_voxels_c,
)
from pore.paths import C_CODE_DIR


VOXEL_C_PATH = C_CODE_DIR / "voxel.so"
VOXEL_C = ctypes.CDLL(str(VOXEL_C_PATH.absolute()))


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

    VOXEL_C.get_single_voxel.restype = ctypes.POINTER(ctypes.c_int * 3)
    voxel_indices_c = VOXEL_C.get_single_voxel(
        ctypes.byref(voxels_x),
        ctypes.byref(voxels_y),
        ctypes.byref(voxels_z),
        index,
        ctypes.byref(initial_single_voxel_indices),
    )

    assert [i for i in voxel_indices_c.contents] == voxel


@pytest.mark.parametrize(
    "voxel_one, voxel_two, is_neighbor",
    [
        (
            (1, 7, 9),
            (2, 7, 9),
            True,
        ),
        (
            (1, 7, 9),
            (2, 8, 9),
            False,
        ),
    ],
)
def test_is_neighbor_voxel(voxel_one, voxel_two, is_neighbor):
    assert is_neighbor_voxel(voxel_one, voxel_two) == is_neighbor


@pytest.mark.parametrize(
    "voxel_one, voxel_two, is_neighbor",
    [
        (
            (1, 7, 9),
            (2, 7, 9),
            True,
        ),
        (
            (1, 7, 9),
            (2, 8, 9),
            False,
        ),
    ],
)
def test_is_neighbor_voxel_c(voxel_one, voxel_two, is_neighbor):
    voxel_one = (ctypes.c_int * 3)(*voxel_one)
    voxel_two = (ctypes.c_int * 3)(*voxel_two)

    is_neighbor_c = VOXEL_C.is_neighbor_voxel(
        ctypes.byref(voxel_one),
        ctypes.byref(voxel_two),
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
def test_breadth_first_search_python(voxels, searchable_indices, neighbor_indices):
    assert breadth_first_search_python(voxels, searchable_indices) == neighbor_indices


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
def test_breadth_first_search_c(voxels, searchable_indices, neighbor_indices):
    assert breadth_first_search_c(voxels, searchable_indices) == neighbor_indices


@pytest.mark.parametrize(
    "query_voxels, reference_voxels, neighbor_voxels",
    [
        (
            (
                np.array([1, 1, 1, 1, 2]),
                np.array([0, 1, 1, 2, 2]),
                np.array([0, 0, 1, 0, 2]),
            ),
            (
                np.array([0, 2]),
                np.array([0, 2]),
                np.array([0, 3]),
            ),
            (
                np.array([1, 2]),
                np.array([0, 2]),
                np.array([0, 2]),
            ),
        ),
    ],
)
def test_get_neighbor_voxels_python(query_voxels, reference_voxels, neighbor_voxels):
    computed_neighbor_voxels = get_neighbor_voxels_python(query_voxels, reference_voxels)
    assert (
        np.allclose(computed_neighbor_voxels[0], neighbor_voxels[0])
        and np.allclose(computed_neighbor_voxels[0], neighbor_voxels[0])
        and np.allclose(computed_neighbor_voxels[0], neighbor_voxels[0])
    )


@pytest.mark.parametrize(
    "query_voxels, reference_voxels, neighbor_voxels",
    [
        (
            (
                np.array([1, 1, 1, 1, 2]),
                np.array([0, 1, 1, 2, 2]),
                np.array([0, 0, 1, 0, 2]),
            ),
            (
                np.array([0, 2]),
                np.array([0, 2]),
                np.array([0, 3]),
            ),
            (
                np.array([1, 2]),
                np.array([0, 2]),
                np.array([0, 2]),
            ),
        ),
    ],
)
def test_get_neighbor_voxels_c(query_voxels, reference_voxels, neighbor_voxels):
    computed_neighbor_voxels = get_neighbor_voxels_c(query_voxels, reference_voxels)
    assert (
        np.allclose(computed_neighbor_voxels[0], neighbor_voxels[0])
        and np.allclose(computed_neighbor_voxels[0], neighbor_voxels[0])
        and np.allclose(computed_neighbor_voxels[0], neighbor_voxels[0])
    )
