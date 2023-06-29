"""
Module containing functions for alignment of coordinates, e.g. along principal axes.
"""


import numpy as np
import biotite.structure as bts


def align_structure(structure: bts.AtomArray) -> tuple[bts.AtomArray, np.ndarray, np.ndarray]:
    """
    Translate the center-of-geometry of a protein to [0,0,0].
    Then align to the principal axes.

    Return the aligned structure as well as the rotation matrix and translation vector used.
    """
    aligned_structure = bts.orient_principal_components(structure)
    _, transformation = bts.superimpose(structure, aligned_structure)

    return aligned_structure, transformation[1], transformation[0]
