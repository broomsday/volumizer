import pytest
import numpy as np
import ctypes

from pore.fib_sphere import fibonacci_sphere, estimate_fibonacci_sphere_samples
from pore.paths import C_CODE_DIR


fib_sphere_c_path = C_CODE_DIR / "fib_sphere.so"
fib_sphere_c = ctypes.CDLL(str(fib_sphere_c_path.absolute()))


@pytest.mark.parametrize(
    "radius, voxel_size, samples",
    [
        (
            1.4,
            1.0,
            40,
        ),
        (
            1.4,
            3.0,
            20,
        ),
        (
            1.8,
            2.0,
            20,
        ),
    ],
)
def test_estimate_fibonacci_sphere_samples(radius, voxel_size, samples):
    assert estimate_fibonacci_sphere_samples(radius, voxel_size) == samples


@pytest.mark.parametrize(
    "radius, x, y, z, samples, coords",
    [
        (
            1.4,
            0.0,
            0.0,
            0.0,
            5,
            {
                "x": [0.84, -0.946, 0.122, 0.781, -0.827],
                "y": [0.0, -0.867, 1.395, -1.018, 0.146],
                "z": [1.12, 0.56, 0.0, -0.56, -1.12],
            },
        ),
        (
            1.8,
            1.0,
            -10.0,
            5.5,
            10,
            {
                "x": [1.785, 0.052, 1.136, 2.045, -0.764, 2.511, 0.554, 0.282, 2.207, 0.275],
                "y": [-10.0, -10.868, -8.447, -11.363, -9.688, -9.039, -11.658, -8.617, -10.441, -10.299],
                "z": [7.12, 6.76, 6.4, 6.04, 5.68, 5.32, 4.96, 4.6, 4.24, 3.88],
            },
        ),
    ],
)
def test_fibonacci_sphere(radius, x, y, z, samples, coords):
    computed_coords = fibonacci_sphere(radius, x, y, z, samples)
    assert (
        np.allclose(coords["x"], computed_coords["x"], atol=0.001)
        and np.allclose(coords["y"], computed_coords["y"], atol=0.001)
        and np.allclose(coords["z"], computed_coords["z"], atol=0.001)
    )


@pytest.mark.parametrize(
    "radius, x, y, z, samples, coords",
    [
        (
            1.4,
            0.0,
            0.0,
            0.0,
            5,
            {
                "x": [0.84, -0.946, 0.122, 0.781, -0.827],
                "y": [0.0, -0.867, 1.395, -1.018, 0.146],
                "z": [1.12, 0.56, 0.0, -0.56, -1.12],
            },
        ),
        (
            1.8,
            1.0,
            -10.0,
            5.5,
            10,
            {
                "x": [1.785, 0.052, 1.136, 2.045, -0.764, 2.511, 0.554, 0.282, 2.207, 0.275],
                "y": [-10.0, -10.868, -8.447, -11.363, -9.688, -9.039, -11.658, -8.617, -10.441, -10.299],
                "z": [7.12, 6.76, 6.4, 6.04, 5.68, 5.32, 4.96, 4.6, 4.24, 3.88],
            },
        ),
    ],
)
def test_fibonacci_sphere_c(radius, x, y, z, samples, coords):
    num_floats = samples * 3
    initial_values = (ctypes.c_float * num_floats)(*([0.0] * num_floats))

    fib_sphere_c.fibonacci_sphere.restype = ctypes.POINTER(ctypes.c_float * num_floats)
    c_coords = fib_sphere_c.fibonacci_sphere(
        ctypes.c_float(radius),
        ctypes.c_float(x),
        ctypes.c_float(y),
        ctypes.c_float(z),
        samples,
        ctypes.byref(initial_values),
    )
    computed_coords_flat = [coord for coord in c_coords.contents]

    computed_coords = {
        "x": computed_coords_flat[:samples],
        "y": computed_coords_flat[samples : samples * 2],
        "z": computed_coords_flat[samples * 2 :],
    }

    assert (
        np.allclose(coords["x"], computed_coords["x"], atol=0.001)
        and np.allclose(coords["y"], computed_coords["y"], atol=0.001)
        and np.allclose(coords["z"], computed_coords["z"], atol=0.001)
    )
