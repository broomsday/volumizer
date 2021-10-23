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
VOXEL_SIZE = 3.0
OCCLUDED_DIMENSION_LIMIT = 4

MONGO_CONNECTION_STRING = "mongodb://127.0.0.1:27017/"