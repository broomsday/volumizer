import importlib

import numpy as np
import pandas as pd
import pytest

from volumizer import native_backend, fib_sphere


@pytest.fixture(autouse=True)
def _clear_backend_cache():
    native_backend.clear_backend_cache()
    yield
    native_backend.clear_backend_cache()


def test_auto_backend_falls_back_to_python_when_native_missing(monkeypatch):
    real_import = importlib.import_module

    def fake_import(module_name):
        if module_name == "volumizer_native":
            raise ImportError("missing native module")
        return real_import(module_name)

    monkeypatch.setenv("VOLUMIZER_BACKEND", "auto")
    monkeypatch.setattr(native_backend.importlib, "import_module", fake_import)
    monkeypatch.setattr(
        native_backend, "_load_native_module_from_local_artifact", lambda: None
    )

    assert native_backend.get_native_module() is None
    assert native_backend.active_backend() == "python"


def test_forced_native_backend_raises_when_missing(monkeypatch):
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
        native_backend.get_native_module()


def test_fibonacci_sphere_native_path_uses_native_module(monkeypatch):
    class FakeNativeModule:
        @staticmethod
        def fibonacci_sphere_points(radius, x, y, z, samples):
            return np.array(
                [
                    [x + radius, y, z],
                    [x, y + radius, z],
                ][:samples]
            )

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    monkeypatch.setattr(
        native_backend.importlib, "import_module", lambda module_name: FakeNativeModule
    )

    coords = fib_sphere.fibonacci_sphere(2.0, 1.0, 3.0, 5.0, 2, backend="native")
    assert coords == {"x": [3.0, 1.0], "y": [3.0, 5.0], "z": [5.0, 5.0]}


def test_add_extra_points_prefers_native_when_requested(monkeypatch):
    class FakeNativeModule:
        @staticmethod
        def fibonacci_sphere_points(radius, x, y, z, samples):
            # return repeated center points, shape (samples, 3)
            return np.array([[x, y, z] for _ in range(samples)], dtype=float)

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    monkeypatch.setattr(
        native_backend.importlib, "import_module", lambda module_name: FakeNativeModule
    )
    monkeypatch.setattr(
        fib_sphere, "add_extra_points_c", lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("C path should not be used"))
    )

    coords = pd.DataFrame.from_dict(
        {"x": [0.0], "y": [0.0], "z": [0.0], "element": ["C"]}
    )
    out_df = fib_sphere.add_extra_points(
        coords,
        voxel_size=2.0,
        performant=True,
        backend="native",
    )

    assert len(out_df) > len(coords)
