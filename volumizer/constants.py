"""
Constants used to avoid magic numbers.
"""

BASE_ATOMIC_RADIUS = 2.0
ATOMIC_RADII = {
    "H": 1.11,
    "He": 1.40,
    "Li": 1.82,
    "Be": 1.53,
    "B": 1.92,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "Na": 2.27,
    "Mg": 1.73,
    "Al": 1.84,
    "Si": 2.10,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "K": 2.75,
    "Ca": 2.31,
    "Ni": 1.63,
    "Cu": 1.40,
    "Zn": 1.39,
    "Se": 1.90,
    "Br": 1.85,
    "I": 1.98,
    "X": BASE_ATOMIC_RADIUS,
}
VOXEL_SIZE = 3.0
OCCLUDED_DIMENSION_LIMIT = 4
MIN_NUM_VOXELS = 4

VOXEL_TYPE_CHAIN_MAP = {"HUB": "A", "POR": "B", "POK": "C", "CAV": "D", "OCC": "E"}
VOXEL_TYPE_ATOM_MAP = {"HUB": "F", "POR": "O", "POK": "N", "CAV": "S", "OCC": "H"}
VOXEL_TYPE_ELEMENT_MAP = {"HUB": "F", "POR": "O", "POK": "N", "CAV": "S", "OCC": "H"}
