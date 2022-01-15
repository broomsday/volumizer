"""
Module containing functions for alignment of coordinates, e.g. along principal axes.
"""


import pandas as pd
import numpy as np
from Bio.PDB import Structure

from pore import pdb


def get_centering_vector(array_coords: np.ndarray) -> np.ndarray:
    """
    Compute the translation vector to place the coordinates center-of-geometry at [0,0,0]
    """
    return np.mean(array_coords, 0)


def order_principal_axis_matrix(eigen_values: list, eigen_vectors: list) -> tuple[np.ndarray, np.ndarray]:
    """
    Put the principal axis eigen values and vectors in order, largest to smallest.
    """
    order = np.argsort(eigen_values)
    eigen_values = eigen_values[order]
    eigen_vectors = eigen_vectors[:, order].transpose() # TODO why transposing here?

    return eigen_values, eigen_vectors


def compute_principal_axis_matrix(array_coords: np.ndarray) -> np.ndarray:
    """
    Compute the eigen values and eigen vectors of the points making up the system
    """
    inertia = np.dot(array_coords.transpose(), array_coords)

    eigen_values, eigen_vectors = np.linalg.eig(inertia)
    eigen_values, eigen_vectors = order_principal_axis_matrix(eigen_values, eigen_vectors)

    return eigen_vectors


def get_numpy_coords(coords: pd.DataFrame) -> np.ndarray:
    """
    Convert the coordinates from a pandas dataframe to a numpy array
    """
    return np.array([
        [row["x"], row["y"], row["z"]] for _, row in coords.iterrows()
    ], float)


def align_coords_to_principal_axes(coords: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """
    Translate the center-of-geometry of a set of coordinates to [0,0,0].
    Then align to the principal axes.
    """
    # convert coordinates into a numpy array for math operations
    array_coords = get_numpy_coords(coords)

    # get rotation and translation matrices for alignment to principal axes
    translation = get_centering_vector(array_coords)
    eigen_vectors = compute_principal_axis_matrix(array_coords)
    print(eigen_vectors)
    # TODO need module to compute rotation matrix (TODO: do this last)
    rotation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # TODO

    return rotation, translation


def apply_rotation_translation(structure: Structure, rotation: np.ndarray, translation: np.ndarray) -> Structure:
    """
    Apply a rotation matrix and translation vector to the coordinates of the supplied structure.
    """
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.transform(rotation, translation)

    return structure

def align_structure(structure: Structure) -> Structure:
    """
    Translate the center-of-geometry of a protein to [0,0,0].
    Then align to the principal axes.
    """

    coords = pdb.get_structure_coords(structure)
    rotation, translation = align_coords_to_principal_axes(coords)
    structure = apply_rotation_translation(structure, rotation, translation)

    return structure