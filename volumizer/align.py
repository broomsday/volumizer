"""
Module containing functions for alignment of coordinates, e.g. along principal axes.
"""


import biotite.structure as bts


def align_structure(
    structure: bts.AtomArray,
) -> bts.AtomArray:
    """
    Translate the center-of-geometry of a protein to [0,0,0] and align to the principal axes.
    """
    return bts.orient_principal_components(structure)
