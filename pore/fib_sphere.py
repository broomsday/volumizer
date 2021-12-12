"""
Functions for creating points on a sphere surface.
"""


import pandas as pd
import numpy as np

from pore import utils
from pore.constants import ATOMIC_RADII, BASE_ATOMIC_RADII


GOLDEN_RATIO = (1 + np.sqrt(5.0)) / 4


def fibonacci_sphere(radius: float, x: float, y: float, z: float, samples: int) -> dict[str, list[float]]:
    """
    Generate `samples` points at `radius` from `x,y,z` such that points are roughly equally spaced.
    """
    x_coords = []
    y_coords = []
    z_coords = []

    # TODO convert this to a dictionary comprehension right from the start once we know it works
    # TODO we could be more performant by generating phi, theta as np matrices/vectors and then the same for x,y,z
    for i in range(samples):
        phi = np.arccos(1 - 2 * (i + 0.5) / samples)
        theta = np.pi * i / GOLDEN_RATIO

        x_coords.append(x + (np.cos(theta) * np.sin(phi)) * radius)
        y_coords.append(y + (np.sin(theta) * np.sin(phi)) * radius)
        z_coords.append(z + (np.cos(phi)) * radius)

    return {"x": x_coords, "y": y_coords, "z": z_coords}


def get_example_point_distance(coords: dict[str, list[float]]) -> float:
    """
    Get the distance between two adjacent points on a fibonacci sphere.
    """
    #i = np.random.randint(0, len(coords["x"]))
    i = 0
    coord_diff = [
        coords["x"][i] - coords["x"][i + 1],
        coords["y"][i] - coords["y"][i + 1],
        coords["z"][i] - coords["z"][i + 1],
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
    vdw_radius = ATOMIC_RADII.get(element, BASE_ATOMIC_RADII)
    radii = [vdw_radius]

    scale_factor = 1
    while min(radii) > voxel_size:
        radii.append(vdw_radius - (voxel_size * scale_factor))
        scale_factor += 1

    return radii


def add_extra_points(coords: pd.DataFrame, voxel_size: float = utils.VOXEL_SIZE) -> pd.DataFrame:
    """
    For each given point which represents the center of a heavy atom,
    add additional points on a surface around that point at one or more radii.
    """
    # TODO elements should be added to the coords dataframe (done before this function is called)
    # TODO use the actual element to get the radii rather than assuming carbon
    radii = get_fibonacci_sphere_radii("C", voxel_size)

    extra_coords = pd.DataFrame.from_dict({"x": [], "y": [], "z": []})
    for radius in radii:
        num_samples = estimate_fibonacci_sphere_samples(radius, voxel_size)
        for _, coord in coords.iterrows():
            extra_points = pd.DataFrame.from_dict(fibonacci_sphere(radius, coord.x, coord.y, coord.z, num_samples))
            extra_coords = extra_coords.append(extra_points, ignore_index=True)

    coords = coords.append(extra_coords, ignore_index=True)
    return coords
