"""
Various utility functions
"""


from pathlib import Path

import gzip
import pandas as pd
import numpy as np

from pore.paths import DOWNLOADED_PDB_DIR, PREPARED_PDB_DIR, ANNOTATED_PDB_DIR, PROTEIN_COMPONENTS_FILE
from pore.constants import VOXEL_SIZE, ATOMIC_RADII, BASE_ATOMIC_RADII
from pore import rcsb


GOLDEN_RATIO = (1 + np.sqrt(5.0)) / 4


def get_downloaded_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb"


def get_prepared_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return PREPARED_PDB_DIR / f"{pdb_id}.pdb"


def get_annotated_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return ANNOTATED_PDB_DIR / f"{pdb_id}.pdb"


def is_pdb_downloaded(pdb_id: str) -> bool:
    """
    Does the downloaded PDB archive exist.
    """
    return get_downloaded_pdb_path(pdb_id).is_file()


def is_pdb_prepared(pdb_id: str) -> bool:
    """
    Has the PDB been fully processed.
    """
    return get_prepared_pdb_path(pdb_id).is_file()


def is_pdb_annotated(pdb_id: str) -> bool:
    """
    Has the PDB been fully processed.
    """
    return get_annotated_pdb_path(pdb_id).is_file()


def decompress_pdb(pdb_id: str) -> None:
    """
    Decompress the gzipped PDB file in the PDB_DIR and save
    the decompressed file into the PROCESSED_PDB_DIR
    """
    zipped_path = DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb1.gz"
    unzipped_path = DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb"

    with gzip.open(zipped_path, mode="rb") as fi:
        pdb_data = fi.read()
    with open(unzipped_path, mode="wb") as fo:
        fo.write(pdb_data)

    zipped_path.unlink()


def setup_dirs():
    """
    Make sure the base data directories exist.
    """
    DOWNLOADED_PDB_DIR.mkdir(parents=True, exist_ok=True)
    PREPARED_PDB_DIR.mkdir(parents=True, exist_ok=True)
    ANNOTATED_PDB_DIR.mkdir(parents=True, exist_ok=True)


def load_protein_components() -> set[str]:
    """
    Load the processed protein components file.
    """
    with open(PROTEIN_COMPONENTS_FILE, mode="r", encoding="utf-8") as fi:
        return set([component.rstrip("\n") for component in fi.readlines()])


def ensure_protein_components_file() -> None:
    """
    Make sure the processed components file exists.
    If not already present, generate it.
    """
    if not PROTEIN_COMPONENTS_FILE.is_file():
        components = rcsb.get_components()
        protein_components = rcsb.build_protein_component_set(components)
        with open(PROTEIN_COMPONENTS_FILE, mode="w", encoding="utf-8") as fo:
            fo.write("\n".join([component for component in protein_components]))


# TODO all this fibonacci extra points stuff should actually be it's own module
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
    Get the distance between a random two points in a set of coords
    """
    i_one = 0  # TODO make random
    i_two = 1  # TODO make random
    coord_diff = [
        coords["x"][i_one] - coords["x"][i_two],
        coords["y"][i_one] - coords["y"][i_two],
        coords["z"][i_one] - coords["z"][i_two],
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

    scale_factor = 2
    while min(radii) > voxel_size:
        radii.append(vdw_radius / scale_factor)
        scale_factor *= 2

    return radii


def add_extra_points(coords: pd.DataFrame, voxel_size: float = VOXEL_SIZE) -> pd.DataFrame:
    """
    For each given point which represents the center of a heavy atom,
    add additional points on a surface around that point at one or more radii.
    """
    # TODO elements should be added to the coords dataframe (done before this function is called)
    # TODO use the actual element to get the radii rather than assuming carbon
    radii = get_fibonacci_sphere_radii("C", voxel_size)

    # TODO when there is an extra radius we start adding coords to the existing additions rather than just the originals
    extra_coords = pd.DataFrame.from_dict({"x": [], "y": [], "z": []})
    for radius in radii:
        num_samples = estimate_fibonacci_sphere_samples(radius, voxel_size)
        for _, coord in coords.iterrows():
            extra_points = pd.DataFrame.from_dict(fibonacci_sphere(radius, coord.x, coord.y, coord.z, num_samples))
            extra_coords = extra_coords.append(extra_points, ignore_index=True)

    coords = coords.append(extra_coords, ignore_index=True)
    return coords