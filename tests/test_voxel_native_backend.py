import importlib

import numpy as np
import pytest

from volumizer import native_backend, voxel
from volumizer.types import VoxelGroup


@pytest.fixture(autouse=True)
def _clear_backend_cache():
    native_backend.clear_backend_cache()
    yield
    native_backend.clear_backend_cache()


def test_get_neighbor_voxels_native_uses_native_module(monkeypatch):
    class FakeNativeModule:
        @staticmethod
        def get_neighbor_voxel_indices(query, reference):
            assert query.shape == (3, 3)
            assert reference.shape == (2, 3)
            return np.array([0, 2], dtype=np.int32)

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    monkeypatch.setattr(
        native_backend.importlib, "import_module", lambda module_name: FakeNativeModule
    )

    query_voxels = (
        np.array([1, 2, 3]),
        np.array([10, 20, 30]),
        np.array([100, 200, 300]),
    )
    reference_voxels = (
        np.array([0, 9]),
        np.array([0, 9]),
        np.array([0, 9]),
    )

    computed = voxel.get_neighbor_voxels_native(query_voxels, reference_voxels)
    assert np.array_equal(computed[0], np.array([1, 3]))
    assert np.array_equal(computed[1], np.array([10, 30]))
    assert np.array_equal(computed[2], np.array([100, 300]))


def test_breadth_first_search_prefers_native_when_requested(monkeypatch):
    monkeypatch.setattr(
        voxel,
        "breadth_first_search_native",
        lambda voxels, searchable_indices: {42},
    )
    monkeypatch.setattr(
        voxel,
        "breadth_first_search_c",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            RuntimeError("C path should not be used")
        ),
    )

    voxels = (
        np.array([0, 1, 2]),
        np.array([0, 0, 0]),
        np.array([0, 0, 0]),
    )
    out = voxel.breadth_first_search(
        voxels, set([0, 1, 2]), performant=True, backend="native"
    )
    assert out == {42}


def test_native_mode_raises_when_native_missing(monkeypatch):
    real_import = importlib.import_module

    def fake_import(module_name):
        if module_name == "volumizer_native":
            raise ImportError("missing native module")
        return real_import(module_name)

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    monkeypatch.setattr(native_backend.importlib, "import_module", fake_import)
    monkeypatch.setattr(
        native_backend, "_load_native_module_from_local_artifact", lambda: None
    )

    with pytest.raises(RuntimeError):
        voxel.get_neighbor_voxels_native(
            (
                np.array([0]),
                np.array([0]),
                np.array([0]),
            ),
            (
                np.array([1]),
                np.array([1]),
                np.array([1]),
            ),
        )


def _make_test_buried_and_exposed_voxels():
    buried_voxels = VoxelGroup(
        voxels=(
            np.array([1, 1]),
            np.array([1, 1]),
            np.array([1, 2]),
        ),
        indices=set([21, 22]),
        num_voxels=2,
        voxel_type="buried",
        volume=16.0,
    )
    exposed_voxels = VoxelGroup(
        voxels=(np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)),
        indices=set(),
        num_voxels=0,
        voxel_type="exposed",
        volume=0.0,
    )
    return buried_voxels, exposed_voxels


def _make_dummy_voxel_grid():
    class DummyGrid:
        x_y_z = np.array([4, 4, 4], dtype=int)
        voxel_centers = np.arange(64 * 3, dtype=float).reshape(64, 3)

    return DummyGrid()


def test_classify_buried_components_native_maps_native_output(monkeypatch):
    class FakeNativeModule:
        @staticmethod
        def classify_buried_components(
            buried_array, exposed_array, grid_dimensions, min_num_voxels, voxel_size
        ):
            assert buried_array.shape == (2, 3)
            assert exposed_array.shape == (0, 3)
            assert np.array_equal(grid_dimensions, np.array([4, 4, 4], dtype=np.int32))
            return {
                "component_type_codes": np.array([2], dtype=np.int32),  # pocket
                "component_offsets": np.array([0, 2], dtype=np.int32),
                "surface_offsets": np.array([0, 1], dtype=np.int32),
                "component_voxel_indices_flat": np.array([0, 1], dtype=np.int32),
                "component_surface_indices_flat": np.array([0], dtype=np.int32),
            }

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    monkeypatch.setattr(
        native_backend.importlib, "import_module", lambda module_name: FakeNativeModule
    )

    buried_voxels, exposed_voxels = _make_test_buried_and_exposed_voxels()
    dummy_grid = _make_dummy_voxel_grid()
    hubs, pores, pockets, cavities, occluded = voxel.classify_buried_components_native(
        buried_voxels, exposed_voxels, dummy_grid
    )

    assert len(hubs) == 0
    assert len(pores) == 0
    assert len(cavities) == 0
    assert len(occluded) == 0
    assert len(pockets) == 1
    assert pockets[0].voxel_type == "pocket"
    assert pockets[0].num_voxels == 2


def test_get_pores_dispatches_to_native_classifier(monkeypatch):
    class FakeNativeModule:
        @staticmethod
        def classify_buried_components(*args, **kwargs):
            return {}

    expected_output = ({1: "hub"}, {2: "pore"}, {3: "pocket"}, {}, {})

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    monkeypatch.setattr(
        native_backend.importlib, "import_module", lambda module_name: FakeNativeModule
    )
    monkeypatch.setattr(
        voxel,
        "classify_buried_components_native",
        lambda buried_voxels, exposed_voxels, voxel_grid: expected_output,
    )
    monkeypatch.setattr(
        voxel,
        "breadth_first_search",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            RuntimeError("Python/C path should not be used")
        ),
    )

    buried_voxels, exposed_voxels = _make_test_buried_and_exposed_voxels()
    dummy_grid = _make_dummy_voxel_grid()
    result = voxel.get_pores_pockets_cavities_occluded(
        buried_voxels, exposed_voxels, dummy_grid, backend="native"
    )

    assert result == expected_output
