import numpy as np
import pytest

from volumizer import native_backend
from volumizer.voxel import (
    breadth_first_search_native,
    breadth_first_search_python,
    get_exposed_and_buried_voxels_native,
    get_exposed_and_buried_voxels_python,
    get_protein_and_solvent_voxels,
    get_neighbor_voxels_native,
    get_neighbor_voxels_python,
)


@pytest.fixture(scope="module", autouse=True)
def _require_native_backend():
    native_backend.clear_backend_cache()
    if native_backend.get_native_module_for_mode("auto") is None:
        pytest.skip("native backend artifact/module not available")


def test_get_neighbor_voxels_native_matches_python():
    query_voxels = (
        np.array([1, 1, 1, 1, 2]),
        np.array([0, 1, 1, 2, 2]),
        np.array([0, 0, 1, 0, 2]),
    )
    reference_voxels = (
        np.array([0, 2]),
        np.array([0, 2]),
        np.array([0, 3]),
    )

    python_neighbor_voxels = get_neighbor_voxels_python(query_voxels, reference_voxels)
    native_neighbor_voxels = get_neighbor_voxels_native(query_voxels, reference_voxels)
    assert np.array_equal(python_neighbor_voxels[0], native_neighbor_voxels[0])
    assert np.array_equal(python_neighbor_voxels[1], native_neighbor_voxels[1])
    assert np.array_equal(python_neighbor_voxels[2], native_neighbor_voxels[2])


def test_breadth_first_search_native_matches_python():
    voxels = (
        np.array([0, 1, 1, 1, 1, 2]),
        np.array([0, 0, 1, 1, 2, 2]),
        np.array([0, 0, 0, 1, 0, 2]),
    )
    searchable_indices = set([0, 1, 2, 3, 4, 5])

    python_neighbor_indices = breadth_first_search_python(voxels, searchable_indices)
    native_neighbor_indices = breadth_first_search_native(voxels, searchable_indices)
    assert python_neighbor_indices == native_neighbor_indices


def test_get_exposed_and_buried_voxels_native_matches_python():
    native_module = native_backend.get_native_module_for_mode("auto")
    if native_module is None or not hasattr(
        native_module, "get_exposed_and_buried_voxel_indices"
    ):
        pytest.skip("native exposed/buried split kernel not available")

    voxel_grid_dimensions = np.array([6, 6, 6], dtype=np.int32)
    binary_grid = np.zeros(voxel_grid_dimensions, dtype=np.int32)
    binary_grid[2, 2, 2] = 1
    binary_grid[2, 2, 3] = 1
    binary_grid[2, 3, 2] = 1
    binary_grid[3, 2, 2] = 1
    binary_grid[3, 3, 3] = 1

    protein_voxels, solvent_voxels = get_protein_and_solvent_voxels(
        binary_grid, voxel_grid_dimensions
    )

    python_exposed, python_buried = get_exposed_and_buried_voxels_python(
        solvent_voxels, protein_voxels, voxel_grid_dimensions
    )
    native_exposed, native_buried = get_exposed_and_buried_voxels_native(
        solvent_voxels, protein_voxels, voxel_grid_dimensions
    )

    assert native_exposed.indices == python_exposed.indices
    assert native_buried.indices == python_buried.indices
    assert native_exposed.num_voxels == python_exposed.num_voxels
    assert native_buried.num_voxels == python_buried.num_voxels
