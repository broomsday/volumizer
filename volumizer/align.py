"""
Module containing functions for alignment of coordinates, e.g. along principal axes.
"""


import biotite.structure as bts
import numpy as np


def align_structure(
    structure: bts.AtomArray,
) -> bts.AtomArray:
    """
    Translate the center-of-geometry of a protein to [0,0,0] and align to the principal axes.
    """
    coords = bts.coord(structure)
    if coords.ndim != 2:
        raise ValueError(f"Expected input shape of (n, 3), got {coords.shape}.")

    row, col = coords.shape
    if (row < 3) or (col != 3):
        raise ValueError(
            f"Expected at least 3 entries, {row} given, and 3 dimensions, {col} given."
        )

    centered = coords - coords.mean(axis=0)
    order = np.array([0, 1, 2], dtype=int)
    identity = np.eye(3)

    for _ in range(50):
        # We only use singular values and right-singular vectors; avoid allocating an NxN U matrix.
        _, sigma, components = np.linalg.svd(centered, full_matrices=False)
        idx = sigma.argsort()[::-1][order]
        ident = np.eye(3)[:, idx]

        v, _, wt = np.linalg.svd(components.T @ ident, full_matrices=False)
        if np.linalg.det(v) * np.linalg.det(wt) < -0:
            v[:, -1] *= -1

        rotation = v @ wt
        if np.isclose(rotation, identity, atol=1e-5).all():
            break
        centered = centered @ rotation

    aligned_structure = structure.copy()
    aligned_structure.coord = centered
    return aligned_structure
