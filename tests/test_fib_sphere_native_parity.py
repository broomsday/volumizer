import numpy as np
import pandas as pd
import pytest

from volumizer import fib_sphere


pytest.importorskip("volumizer_native")


@pytest.mark.parametrize(
    "radius, x, y, z, samples",
    [
        (1.4, 0.0, 0.0, 0.0, 20),
        (1.8, 1.0, -10.0, 5.5, 40),
    ],
)
def test_fibonacci_sphere_native_matches_python(radius, x, y, z, samples):
    python_coords = fib_sphere.fibonacci_sphere(
        radius, x, y, z, samples, backend="python"
    )
    native_coords = fib_sphere.fibonacci_sphere(
        radius, x, y, z, samples, backend="native"
    )

    assert np.allclose(python_coords["x"], native_coords["x"], atol=1e-5)
    assert np.allclose(python_coords["y"], native_coords["y"], atol=1e-5)
    assert np.allclose(python_coords["z"], native_coords["z"], atol=1e-5)


def test_add_extra_points_native_matches_python():
    coords = pd.DataFrame.from_dict(
        {
            "x": [0.0, 1.5, -2.0],
            "y": [0.0, -3.5, 4.0],
            "z": [0.0, 2.5, -1.0],
            "element": ["C", "N", "O"],
        }
    )
    voxel_size = 2.0

    python_df = fib_sphere.add_extra_points(coords, voxel_size, backend="python")
    native_df = fib_sphere.add_extra_points(coords, voxel_size, backend="native")

    python_df = python_df.sort_values(
        by=["element", "x", "y", "z"], ignore_index=True
    )
    native_df = native_df.sort_values(
        by=["element", "x", "y", "z"], ignore_index=True
    )

    assert len(python_df) == len(native_df)
    assert np.array_equal(python_df["element"].to_numpy(), native_df["element"].to_numpy())
    assert np.allclose(
        python_df[["x", "y", "z"]].to_numpy(),
        native_df[["x", "y", "z"]].to_numpy(),
        atol=1e-5,
    )
