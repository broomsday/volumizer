"""
Constants used to avoid magic numbers.
"""

PDB_ID_LENGTH = 4

RCSB_CLUSTER_URL = "https://cdn.rcsb.org/resources/sequence/clusters/bc-90.out"
RCSB_CCD_URL = "https://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif"
RCSB_BIOUNIT_URL = "https://ftp.wwpdb.org/pub/pdb/data/biounit/PDB/divided/"
RCSB_STRUCTURE_URL = "https://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/"

VOXEL_ATOM_NAMES = frozenset(
    [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CG1",
        "CG2",
        "CD",
        "CD1",
        "CD2",
        "CE",
        "CE1",
        "CE2",
        "CE3",
        "CZ",
        "CZ2",
        "CZ3",
        "CH2",
        "ND1",
        "ND2",
        "NE",
        "NE1",
        "NE2",
        "NZ",
        "OG",
        "OG1",
        "OD1",
        "OD2",
        "OE1",
        "OE2",
        "OH",
        "SG",
        "SD",
    ]
)
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
VOXEL_SIZE = 2.0
OCCLUDED_DIMENSION_LIMIT = 4
POCKET_VOLUME_THRESHOLD = 20
DIAGONAL_NEIGHBORS = True

VOXEL_TYPE_CHAIN_MAP = {"POK": "A", "POR": "B", "CAV": "C", "OCC": "D"}
VOXEL_TYPE_ATOM_MAP = {"POK": "N", "POR": "O", "CAV": "S", "OCC": "H"}
