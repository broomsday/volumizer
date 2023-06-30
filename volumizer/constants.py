"""
Constants used to avoid magic numbers.
"""

PDB_ID_LENGTH = 4

RCSB_CLUSTER_URL = "https://cdn.rcsb.org/resources/sequence/clusters/bc-90.out"
RCSB_CCD_URL = "https://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif"
RCSB_BIOUNIT_URL = "https://ftp.wwpdb.org/pub/pdb/data/biounit/PDB/divided/"
RCSB_BUNDLE_URL = "https://files.rcsb.org/pub/pdb/compatible/pdb_bundle/"
RCSB_STRUCTURE_URL = "https://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/"
RCSB_GENERAL_INFO_URL = "https://data.rcsb.org/rest/v1/core/entry/"
RCSB_ASSEMBLY_INFO_URL = "https://data.rcsb.org/rest/v1/core/assembly/"
RCSB_CONTACT_RETRIES = 10

# TODO: this is only used in some optional downstream coordinate cleaning, but I suspect we should just ditch it
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
VOXEL_SIZE = 3.0
OCCLUDED_DIMENSION_LIMIT = 4
MIN_NUM_VOXELS = 4

VOXEL_TYPE_CHAIN_MAP = {"HUB": "A", "POR": "B", "POK": "C", "CAV": "D", "OCC": "E"}
VOXEL_TYPE_ATOM_MAP = {"HUB": "F", "POR": "O", "POK": "N", "CAV": "S", "OCC": "H"}

RESIDUE_LETTER_CONVERSION = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}
SEQUENCE_IDENTITY_CUTOFF = 0.90

STRIDE_CODES = {
    "helix": frozenset(["H", "G", "I"]),
    "strand": frozenset(["E", "B"]),
}
