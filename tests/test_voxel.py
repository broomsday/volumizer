import pytest
import numpy as np
import ctypes

from volumizer.voxel import (
    _get_surface_components,
    _is_wrapped_hub_direction_spread,
    get_single_voxel,
    is_neighbor_voxel,
    breadth_first_search_python,
    breadth_first_search_c,
    get_neighbor_voxels_python,
    get_neighbor_voxels_c,
)
from volumizer.paths import C_CODE_DIR


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
    "sorted_direction_bucket_counts, expected",
    [
        ([2011, 1972, 1893, 1706, 1552, 1551], True),
        ([95, 65, 60, 38, 23, 22], False),
        ([107, 63, 0, 0, 0, 0], False),
    ],
)
def test_is_wrapped_hub_direction_spread(
    sorted_direction_bucket_counts, expected
):
    assert _is_wrapped_hub_direction_spread(sorted_direction_bucket_counts) is expected


def test_surface_component_connectivity_modes_distinguish_6_18_and_custom18():
    buried_voxels = (
        np.array([0, 1, 1]),
        np.array([0, 1, 0]),
        np.array([0, 0, 0]),
    )
    direct_surface_indices = set([0, 1])
    support_surface_indices = set([0, 1])

    assert len(_get_surface_components(
        direct_surface_indices,
        buried_voxels,
        support_indices=support_surface_indices,
        connectivity_mode="6",
    )) == 2
    assert len(_get_surface_components(
        direct_surface_indices,
        buried_voxels,
        support_indices=support_surface_indices,
        connectivity_mode="18",
    )) == 1
    assert len(_get_surface_components(
        direct_surface_indices,
        buried_voxels,
        support_indices=support_surface_indices,
        connectivity_mode="26",
    )) == 1
    assert len(_get_surface_components(
        direct_surface_indices,
        buried_voxels,
        support_indices=support_surface_indices,
        connectivity_mode="custom18",
    )) == 2


def test_surface_component_custom18_uses_supported_shell_diagonal():
    buried_voxels = (
        np.array([0, 1, 1]),
        np.array([0, 1, 0]),
        np.array([0, 0, 0]),
    )
    direct_surface_indices = set([0, 1])
    support_surface_indices = set([0, 1, 2])

    assert len(_get_surface_components(
        direct_surface_indices,
        buried_voxels,
        support_indices=support_surface_indices,
    )) == 1


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
