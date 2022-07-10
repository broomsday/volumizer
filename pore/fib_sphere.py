"""
Functions for creating points on a sphere surface.
"""


import ctypes

import pandas as pd
import numpy as np

from pore import utils
from pore.constants import ATOMIC_RADII, BASE_ATOMIC_RADIUS
from pore.paths import C_CODE_DIR


GOLDEN_RATIO = (1 + np.sqrt(5.0)) / 4

if utils.using_performant():
    FIB_SPHERE_C_PATH = C_CODE_DIR / "fib_sphere.so"
    FIB_SPHERE_C = ctypes.CDLL(str(FIB_SPHERE_C_PATH.absolute()))


def fibonacci_sphere(radius: float, x: float, y: float, z: float, samples: int) -> dict[str, list[float]]:
    """
    Generate `samples` points at `radius` from `x,y,z` such that points are roughly equally spaced.
    """
    # this function could likely be made more performant by using numpy more intelligently
    x_coords = []
    y_coords = []
    z_coords = []

    for i in range(samples):
        phi = np.arccos(1 - 2 * (i + 0.5) / samples)
        theta = (np.pi * i) / GOLDEN_RATIO

        x_coords.append(x + (np.cos(theta) * np.sin(phi)) * radius)
        y_coords.append(y + (np.sin(theta) * np.sin(phi)) * radius)
        z_coords.append(z + (np.cos(phi)) * radius)

    return {"x": x_coords, "y": y_coords, "z": z_coords}


def get_example_point_distance(coords: dict[str, list[float]]) -> float:
    """
    Get the distance between two adjacent points on a fibonacci sphere.
    Due to the way the spiral traverses the sphere, for an unknown number of points on the sphere,
    we can only be certain that the first and second, or last and second last points are "adjacent"
    """
    coord_diff = [
        coords["x"][0] - coords["x"][1],
        coords["y"][0] - coords["y"][1],
        coords["z"][0] - coords["z"][1],
    ]
    return np.linalg.norm(coord_diff)


def estimate_fibonacci_sphere_samples(radius: float, voxel_size: float) -> int:
    """
    Based on the voxel-size and radius being used, estimate how many samples will be needed.
    Smaller voxel-size requires more samples
    """
    num_samples = 20
    coords = fibonacci_sphere(radius, 0, 0, 0, num_samples)
    while get_example_point_distance(coords) > voxel_size:
        num_samples *= 2
        coords = fibonacci_sphere(radius, 0, 0, 0, num_samples)

    return num_samples


def get_fibonacci_sphere_radii(element: str, voxel_size: float) -> list[float]:
    """
    Based on the voxel-size being used, estimate the size and number of radii needed.
    Smaller voxel-size requires more radii.

    We always place a single radii at the VDW radius for the element corresponding to this
    atom, and only use additional, smaller radii when needed.
    """
    vdw_radius = ATOMIC_RADII.get(element, BASE_ATOMIC_RADIUS)
    radii = [vdw_radius]

    scale_factor = 1
    while min(radii) > voxel_size:
        radii.append(vdw_radius - (voxel_size * scale_factor))
        scale_factor += 1

    return radii


def add_extra_points_python(coords: pd.DataFrame, voxel_size: float = utils.VOXEL_SIZE) -> pd.DataFrame:
    """
    For each given point which represents the center of a heavy atom,
    add additional points on a surface around that point at one or more radii.
    """
    elemental_radii = {element: get_fibonacci_sphere_radii(element, voxel_size) for element in set(coords["element"])}

    extra_x, extra_y, extra_z, extra_element = [], [], [], []
    for element, radii in elemental_radii.items():
        for radius in radii:
            num_samples = estimate_fibonacci_sphere_samples(radius, voxel_size)
            for _, coord in coords[coords["element"] == element].iterrows():
                extra_points = fibonacci_sphere(radius, coord.x, coord.y, coord.z, num_samples)
                extra_x += extra_points["x"]
                extra_y += extra_points["y"]
                extra_z += extra_points["z"]
                extra_element += [element] * num_samples

    extra_coords = pd.DataFrame.from_dict({"x": extra_x, "y": extra_y, "z": extra_z, "element": extra_element})
    coords = pd.concat([coords, extra_coords], ignore_index=True)

    return coords


def add_extra_points_c(coords: pd.DataFrame, voxel_size: float = utils.VOXEL_SIZE) -> pd.DataFrame:
    """
    For each given point which represents the center of a heavy atom,
    add additional points on a surface around that point at one or more radii.
    """
    elemental_radii = {element: get_fibonacci_sphere_radii(element, voxel_size) for element in set(coords["element"])}

    extra_x, extra_y, extra_z, extra_element = [], [], [], []
    for element, radii in elemental_radii.items():
        for radius in radii:
            num_samples = estimate_fibonacci_sphere_samples(radius, voxel_size)
            num_floats = num_samples * 3
            initial_values = (ctypes.c_float * num_floats)(*([0.0] * num_floats))
            FIB_SPHERE_C.fibonacci_sphere.restype = ctypes.POINTER(ctypes.c_float * num_floats)
            for _, coord in coords[coords["element"] == element].iterrows():
                c_coords = FIB_SPHERE_C.fibonacci_sphere(
                    ctypes.c_float(radius),
                    ctypes.c_float(coord.x),
                    ctypes.c_float(coord.y),
                    ctypes.c_float(coord.z),
                    num_samples,
                    ctypes.byref(initial_values),
                )
                computed_coords_flat = [coord for coord in c_coords.contents]

                extra_x += computed_coords_flat[:num_samples]
                extra_y += computed_coords_flat[num_samples : num_samples * 2]
                extra_z += computed_coords_flat[num_samples * 2 :]
                extra_element += [element] * num_samples

    extra_coords = pd.DataFrame.from_dict({"x": extra_x, "y": extra_y, "z": extra_z, "element": extra_element})
    coords = pd.concat([coords, extra_coords])

    return coords


def add_extra_points(
    coords: pd.DataFrame, voxel_size: float = utils.VOXEL_SIZE, performant: bool = utils.using_performant()
) -> pd.DataFrame:
    """
    For each given point which represents the center of a heavy atom,
    add additional points on a surface around that point at one or more radii.
    """
    if performant:
        return add_extra_points_c(coords, voxel_size)

    return add_extra_points_python(coords, voxel_size)
