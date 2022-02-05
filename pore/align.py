"""
Module containing functions for alignment of coordinates, e.g. along principal axes.
"""


import pandas as pd
import numpy as np
from scipy.spatial.transform import Rotation
from Bio.PDB import Structure

from pore import pdb


def get_centering_vector(array_coords: np.ndarray) -> np.ndarray:
    """
    Compute the translation vector to place the coordinates center-of-geometry at [0,0,0]
    """
    return np.mean(array_coords, axis=0) * -1


def compute_principal_axis_matrix(array_coords: np.ndarray) -> np.ndarray:
    """
    Compute the eigen values and eigen vectors of the points making up the system
    """
    inertia = np.dot(array_coords.transpose(), array_coords)
    _, eigen_vectors = np.linalg.eig(inertia)

    return eigen_vectors


def df_to_numpy_coords(coords: pd.DataFrame) -> np.ndarray:
    """
    Convert the coordinates from a pandas dataframe to a numpy array
    """
    return np.array([[row["x"], row["y"], row["z"]] for _, row in coords.iterrows()], float)


def get_principal_axis_alignment_translation_rotation(array_coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Get the translation vector and rotation matrix to align the coords to their principal axis
    and place the center-of-geometry at the origin.
    """
    coords = array_coords.copy()
    translation = get_centering_vector(coords)
    coords += translation

    # get rotation matrix for alignment to principal axes
    eigen_vectors = compute_principal_axis_matrix(coords)
    identity = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    rotation = Rotation.align_vectors(identity, eigen_vectors)[0].as_matrix()

    return rotation, translation


def rotate_and_translate_structure(
    structure: Structure.Structure, rotation: np.ndarray, translation: np.ndarray
) -> Structure.Structure:
    """
    Apply a rotation matrix and translation vector to the coordinates of the supplied structure.
    """
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.transform(rotation, translation)

    return structure


def rotate_and_translate_coords(coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    """
    Apply a rotation matrix and translation vector to the supplied coordinates.
    """
    rotation = Rotation.from_matrix(rotation)
    coords = rotation.apply(coords)
    coords += translation

    return coords


def align_structure(structure: Structure.Structure) -> tuple[Structure.Structure, np.ndarray, np.ndarray]:
    """
    Translate the center-of-geometry of a protein to [0,0,0].
    Then align to the principal axes.

    Return the aligned structure as well as the rotation matrix and translation vector used.
    """

    coords = df_to_numpy_coords(pdb.get_structure_coords(structure))
    rotation, translation = get_principal_axis_alignment_translation_rotation(coords)
    structure = rotate_and_translate_structure(structure, rotation, translation)

    return structure, rotation, translation
