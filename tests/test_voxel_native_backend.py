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


def _make_first_shell_test_inputs() -> tuple[VoxelGroup, VoxelGroup]:
    exposed_voxels = VoxelGroup(
        voxels=(np.array([1, 2, 3]), np.array([1, 1, 1]), np.array([1, 1, 1])),
        indices=set([21, 37, 53]),
        num_voxels=3,
        voxel_type="exposed",
        volume=24.0,
    )
    buried_voxels = VoxelGroup(
        voxels=(np.array([0, 3]), np.array([1, 1]), np.array([1, 1])),
        indices=set([5, 53]),
        num_voxels=2,
        voxel_type="buried",
        volume=16.0,
    )

    return exposed_voxels, buried_voxels


def test_get_first_shell_exposed_voxels_native_prefers_specialized_selection(monkeypatch):
    class FakeNativeModule:
        @staticmethod
        def get_first_shell_exposed_selection(
            exposed_x,
            exposed_y,
            exposed_z,
            buried_x,
            buried_y,
            buried_z,
            grid_dimensions,
        ):
            assert exposed_x.dtype == np.int64
            assert exposed_y.dtype == np.int64
            assert exposed_z.dtype == np.int64
            assert buried_x.dtype == np.int64
            assert buried_y.dtype == np.int64
            assert buried_z.dtype == np.int64
            assert np.array_equal(grid_dimensions, np.array([4, 4, 4], dtype=np.int32))
            return {
                "neighbor_indices": np.array([0, 2], dtype=np.int32),
                "flat_indices": np.array([21, 53], dtype=np.int64),
            }

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    monkeypatch.setattr(
        native_backend.importlib, "import_module", lambda module_name: FakeNativeModule
    )

    exposed_voxels, buried_voxels = _make_first_shell_test_inputs()

    class DummyGrid:
        x_y_z = np.array([4, 4, 4], dtype=int)

    first_shell = voxel.get_first_shell_exposed_voxels(
        exposed_voxels,
        buried_voxels,
        DummyGrid(),
        backend="native",
    )

    assert first_shell.num_voxels == 2
    assert np.array_equal(first_shell.voxels[0], np.array([1, 3]))
    assert np.array_equal(first_shell.voxels[1], np.array([1, 1]))
    assert np.array_equal(first_shell.voxels[2], np.array([1, 1]))
    assert first_shell.indices == set([21, 53])


def test_get_first_shell_exposed_voxels_native_falls_back_to_legacy_kernel(monkeypatch):
    class FakeNativeModule:
        @staticmethod
        def get_first_shell_exposed_indices(exposed, buried, grid_dimensions):
            assert exposed.shape == (3, 3)
            assert buried.shape == (2, 3)
            assert np.array_equal(grid_dimensions, np.array([4, 4, 4], dtype=np.int32))
            return np.array([0, 2], dtype=np.int32)

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    monkeypatch.setattr(
        native_backend.importlib, "import_module", lambda module_name: FakeNativeModule
    )

    exposed_voxels, buried_voxels = _make_first_shell_test_inputs()

    class DummyGrid:
        x_y_z = np.array([4, 4, 4], dtype=int)

    first_shell = voxel.get_first_shell_exposed_voxels(
        exposed_voxels,
        buried_voxels,
        DummyGrid(),
        backend="native",
    )

    assert first_shell.num_voxels == 2
    assert np.array_equal(first_shell.voxels[0], np.array([1, 3]))
    assert np.array_equal(first_shell.voxels[1], np.array([1, 1]))
    assert np.array_equal(first_shell.voxels[2], np.array([1, 1]))
    assert first_shell.indices == set([21, 53])


def test_get_exposed_and_buried_voxels_native_uses_native_module(monkeypatch):
    class FakeNativeModule:
        @staticmethod
        def get_exposed_and_buried_voxel_indices(solvent, protein, grid_dimensions):
            assert solvent.shape == (3, 3)
            assert protein.shape == (2, 3)
            assert np.array_equal(grid_dimensions, np.array([4, 4, 4], dtype=np.int32))
            return {
                "exposed_indices": np.array([0, 2], dtype=np.int32),
                "buried_indices": np.array([1], dtype=np.int32),
            }

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    monkeypatch.setattr(
        native_backend.importlib, "import_module", lambda module_name: FakeNativeModule
    )

    solvent_voxels = VoxelGroup(
        voxels=(np.array([0, 1, 2]), np.array([0, 1, 2]), np.array([1, 1, 1])),
        indices=set([1, 21, 41]),
        num_voxels=3,
        voxel_type="solvent",
        volume=24.0,
    )
    protein_voxels = VoxelGroup(
        voxels=(np.array([1, 2]), np.array([0, 3]), np.array([0, 3])),
        indices=set([16, 47]),
        num_voxels=2,
        voxel_type="protein",
        volume=16.0,
    )

    exposed, buried = voxel.get_exposed_and_buried_voxels(
        solvent_voxels,
        protein_voxels,
        np.array([4, 4, 4], dtype=np.int32),
        backend="native",
    )

    assert exposed.num_voxels == 2
    assert buried.num_voxels == 1
    assert exposed.indices == set([1, 41])
    assert buried.indices == set([21])
    assert np.array_equal(exposed.voxels[0], np.array([0, 2]))
    assert np.array_equal(buried.voxels[0], np.array([1]))


def test_get_exposed_and_buried_voxels_native_falls_back_when_missing(monkeypatch):
    class FakeNativeModule:
        @staticmethod
        def get_neighbor_voxel_indices(*args, **kwargs):
            return np.array([], dtype=np.int32)

    expected_output = ("exposed", "buried")

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    monkeypatch.setattr(
        native_backend.importlib, "import_module", lambda module_name: FakeNativeModule
    )
    monkeypatch.setattr(
        voxel,
        "get_exposed_and_buried_voxels_python",
        lambda solvent_voxels, protein_voxels, voxel_grid_dimensions: expected_output,
    )

    solvent_voxels = VoxelGroup(
        voxels=(np.array([0]), np.array([0]), np.array([0])),
        indices=set([0]),
        num_voxels=1,
        voxel_type="solvent",
        volume=8.0,
    )
    protein_voxels = VoxelGroup(
        voxels=(np.array([1]), np.array([1]), np.array([1])),
        indices=set([21]),
        num_voxels=1,
        voxel_type="protein",
        volume=8.0,
    )

    result = voxel.get_exposed_and_buried_voxels(
        solvent_voxels,
        protein_voxels,
        np.array([4, 4, 4], dtype=np.int32),
        backend="native",
    )
    assert result == expected_output


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


def test_classify_buried_components_native_records_kernel_substage_timings(monkeypatch):
    class FakeNativeModule:
        @staticmethod
        def classify_buried_components(
            buried_array, exposed_array, grid_dimensions, min_num_voxels, voxel_size
        ):
            assert buried_array.shape == (2, 3)
            assert exposed_array.shape == (0, 3)
            assert np.array_equal(grid_dimensions, np.array([4, 4, 4], dtype=np.int32))
            return {
                "component_type_codes": np.array([2], dtype=np.int32),
                "component_offsets": np.array([0, 2], dtype=np.int32),
                "surface_offsets": np.array([0, 1], dtype=np.int32),
                "component_voxel_indices_flat": np.array([0, 1], dtype=np.int32),
                "component_surface_indices_flat": np.array([0], dtype=np.int32),
                "kernel_stage_timings_seconds": {
                    "bfs_component_expansion": 0.123,
                    "classify_component_type": 0.045,
                },
            }

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    monkeypatch.setattr(
        native_backend.importlib, "import_module", lambda module_name: FakeNativeModule
    )

    buried_voxels, exposed_voxels = _make_test_buried_and_exposed_voxels()
    dummy_grid = _make_dummy_voxel_grid()
    stage_timings = {}

    voxel.classify_buried_components_native(
        buried_voxels,
        exposed_voxels,
        dummy_grid,
        stage_timings=stage_timings,
    )

    assert stage_timings["classify_components_native_kernel"] > 0.0
    assert stage_timings["classify_components_native_mapping"] > 0.0
    assert stage_timings["classify_components_native_kernel_bfs_component_expansion"] == pytest.approx(
        0.123
    )
    assert stage_timings["classify_components_native_kernel_classify_component_type"] == pytest.approx(
        0.045
    )


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
