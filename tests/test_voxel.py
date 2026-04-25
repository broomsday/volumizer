import pytest
import numpy as np
import ctypes

from volumizer import utils
from volumizer.voxel import (
    _get_surface_components,
    _is_wrapped_hub_direction_spread,
    compute_voxel_indices,
    get_agglomerated_type,
    get_single_voxel,
    is_neighbor_voxel,
    breadth_first_search_python,
    breadth_first_search_c,
    get_neighbor_voxels_python,
    get_neighbor_voxels_c,
    refine_voxel_group_annotation,
)
from volumizer.paths import C_CODE_DIR
from volumizer.types import VoxelGroup


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


def _make_linear_two_mouth_component(num_buried_voxels: int):
    buried_voxels = (
        np.arange(1, num_buried_voxels + 1, dtype=np.int64),
        np.ones(num_buried_voxels, dtype=np.int64),
        np.ones(num_buried_voxels, dtype=np.int64),
    )
    exposed_voxels = (
        np.array([0, num_buried_voxels + 1], dtype=np.int64),
        np.array([1, 1], dtype=np.int64),
        np.array([1, 1], dtype=np.int64),
    )
    grid_dimensions = np.array([num_buried_voxels + 3, 3, 3], dtype=np.int64)
    query_indices = set(range(num_buried_voxels))
    return query_indices, buried_voxels, exposed_voxels, grid_dimensions


def _make_necked_pocket_voxel_group(
    core_edge: int,
) -> tuple[VoxelGroup, tuple[np.ndarray, np.ndarray, np.ndarray], np.ndarray]:
    shell_coords = [(1, 1, 1), (2, 1, 1), (3, 1, 1)]
    neck_coords = [(4, 1, 1), (5, 1, 1), (6, 1, 1)]
    core_coords = [
        (x_coord, y_coord, z_coord)
        for x_coord in range(7, 7 + core_edge)
        for y_coord in range(1, 1 + core_edge)
        for z_coord in range(1, 1 + core_edge)
    ]
    all_coords = shell_coords + neck_coords + core_coords

    voxels = (
        np.array([coord[0] for coord in all_coords], dtype=np.int64),
        np.array([coord[1] for coord in all_coords], dtype=np.int64),
        np.array([coord[2] for coord in all_coords], dtype=np.int64),
    )
    shell_voxels = (
        np.array([coord[0] for coord in shell_coords], dtype=np.int64),
        np.array([coord[1] for coord in shell_coords], dtype=np.int64),
        np.array([coord[2] for coord in shell_coords], dtype=np.int64),
    )
    exposed_voxels = (
        np.array([1, 2, 3], dtype=np.int64),
        np.array([0, 0, 0], dtype=np.int64),
        np.array([1, 1, 1], dtype=np.int64),
    )
    grid_dimensions = np.array(
        [7 + core_edge + 2, core_edge + 3, core_edge + 3],
        dtype=np.int64,
    )

    return (
        VoxelGroup(
            voxels=voxels,
            indices=compute_voxel_indices(voxels, grid_dimensions),
            num_voxels=len(all_coords),
            surface_indices=compute_voxel_indices(shell_voxels, grid_dimensions),
            voxel_type="pocket",
        ),
        exposed_voxels,
        grid_dimensions,
    )


def _empty_exposed_voxels() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    empty = np.array([], dtype=np.int64)
    return (empty, empty, empty)


def _make_voxel_group_from_coords(
    coords: list[tuple[int, int, int]],
    grid_dimensions: np.ndarray,
    voxel_type: str,
) -> VoxelGroup:
    voxels = (
        np.array([coord[0] for coord in coords], dtype=np.int64),
        np.array([coord[1] for coord in coords], dtype=np.int64),
        np.array([coord[2] for coord in coords], dtype=np.int64),
    )
    return VoxelGroup(
        voxels=voxels,
        indices=compute_voxel_indices(voxels, grid_dimensions),
        num_voxels=len(coords),
        voxel_type=voxel_type,
    )


def _make_single_surface_pore_voxel_group() -> tuple[VoxelGroup, np.ndarray]:
    grid_dimensions = np.array([6, 5, 5], dtype=np.int64)
    coords = [
        (x_coord, y_coord, z_coord)
        for x_coord in range(4)
        for y_coord in range(1, 3)
        for z_coord in range(1, 3)
    ]
    return _make_voxel_group_from_coords(coords, grid_dimensions, "pore"), grid_dimensions


def _make_one_narrow_one_wide_pore_voxel_group() -> tuple[VoxelGroup, np.ndarray]:
    grid_dimensions = np.array([8, 5, 4], dtype=np.int64)
    chamber_coords = [
        (x_coord, y_coord, z_coord)
        for x_coord in range(4)
        for y_coord in range(1, 4)
        for z_coord in range(1, 3)
    ]
    neck_coords = [(x_coord, 2, 1) for x_coord in range(4, 8)]
    coords = chamber_coords + neck_coords
    return _make_voxel_group_from_coords(coords, grid_dimensions, "pore"), grid_dimensions


def _make_two_mouth_pore_voxel_group(
    num_buried_voxels: int,
) -> tuple[VoxelGroup, tuple[np.ndarray, np.ndarray, np.ndarray], np.ndarray]:
    _, buried_voxels, exposed_voxels, grid_dimensions = _make_linear_two_mouth_component(
        num_buried_voxels
    )
    voxel_group = VoxelGroup(
        voxels=buried_voxels,
        indices=compute_voxel_indices(buried_voxels, grid_dimensions),
        num_voxels=num_buried_voxels,
        voxel_type="pore",
    )
    return voxel_group, exposed_voxels, grid_dimensions


def test_refine_voxel_group_annotation_promotes_large_necked_pocket_to_cavity():
    voxel_group, exposed_voxels, grid_dimensions = _make_necked_pocket_voxel_group(
        core_edge=13
    )
    refined = refine_voxel_group_annotation(
        voxel_group,
        exposed_voxels,
        grid_dimensions,
    )

    assert refined.voxel_type == "cavity"
    assert refined.surface_indices == set()


def test_refine_voxel_group_annotation_keeps_low_ratio_necked_pocket_as_pocket():
    voxel_group, exposed_voxels, grid_dimensions = _make_necked_pocket_voxel_group(
        core_edge=6
    )
    refined = refine_voxel_group_annotation(
        voxel_group,
        exposed_voxels,
        grid_dimensions,
    )

    assert refined.voxel_type == "pocket"
    assert len(refined.surface_indices) > 0


def test_refine_voxel_group_annotation_demotes_single_surface_pore_to_pocket():
    voxel_group, grid_dimensions = _make_single_surface_pore_voxel_group()

    refined = refine_voxel_group_annotation(
        voxel_group,
        _empty_exposed_voxels(),
        grid_dimensions,
    )

    assert refined.voxel_type == "pocket"
    assert len(refined.surface_indices) > 0


def test_refine_voxel_group_annotation_demotes_one_narrow_one_wide_pore_to_pocket():
    voxel_group, grid_dimensions = _make_one_narrow_one_wide_pore_voxel_group()

    refined = refine_voxel_group_annotation(
        voxel_group,
        _empty_exposed_voxels(),
        grid_dimensions,
    )

    assert refined.voxel_type == "pocket"
    assert len(refined.surface_indices) > 0


def test_refine_voxel_group_annotation_promotes_two_narrow_mouth_pore_to_cavity():
    voxel_group, exposed_voxels, grid_dimensions = _make_two_mouth_pore_voxel_group(6)

    refined = refine_voxel_group_annotation(
        voxel_group,
        exposed_voxels,
        grid_dimensions,
    )

    assert refined.voxel_type == "cavity"
    assert refined.surface_indices == set()


@pytest.mark.parametrize(
    "num_buried_voxels, gap_threshold, expected_type",
    [
        (4, -1, "pore"),
        (4, 0, "pocket"),
        (5, 0, "pore"),
        (5, 1, "pocket"),
        (6, 1, "pore"),
        (6, 2, "pocket"),
    ],
)
def test_get_agglomerated_type_merges_close_two_mouth_pores(
    monkeypatch,
    num_buried_voxels,
    gap_threshold,
    expected_type,
):
    monkeypatch.setattr(
        utils,
        "DIRECT_SURFACE_MOUTH_MERGE_GAP_VOXELS",
        gap_threshold,
    )
    (
        query_indices,
        buried_voxels,
        exposed_voxels,
        grid_dimensions,
    ) = _make_linear_two_mouth_component(num_buried_voxels)

    _, agglomerated_type = get_agglomerated_type(
        query_indices,
        buried_voxels,
        exposed_voxels,
        grid_dimensions,
    )

    assert agglomerated_type == expected_type


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
