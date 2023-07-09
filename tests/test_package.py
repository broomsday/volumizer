from pathlib import Path
import pytest

from volumizer import volumizer, pdb, utils
from volumizer.paths import TEST_DIR

TEST_PDB_DIR = TEST_DIR / "pdbs"
VOLUME_TEST_DATA = [
    {
        "pdb_file": TEST_PDB_DIR / "cavity.pdb",
        "cavity_volume": 184.0,
        "pocket_volume": None,
        "pore_volume": None,
        "hub_volume": None,
    },
    {
        "pdb_file": TEST_PDB_DIR / "pocket.pdb",
        "cavity_volume": None,
        "pocket_volume": 840.0,
        "pore_volume": None,
        "hub_volume": None,
    },
    {
        "pdb_file": TEST_PDB_DIR / "pore.pdb",
        "cavity_volume": None,
        "pocket_volume": None,
        "pore_volume": 1904.0,
        "hub_volume": None,
    },
    {
        "pdb_file": TEST_PDB_DIR / "hub.pdb",
        "cavity_volume": None,
        "pocket_volume": None,
        "pore_volume": None,
        "hub_volume": 2176.0,
    },
]
COMPONENT_TEST_PDB = TEST_PDB_DIR / "4jpn.pdb"


@pytest.mark.parametrize(
    "pdb_file, cavity_volume, pocket_volume, pore_volume, hub_volume",
    [(test_data["pdb_file"], test_data["cavity_volume"], test_data["pocket_volume"],test_data["pore_volume"],test_data["hub_volume"]) for test_data in VOLUME_TEST_DATA
    ],
)
def test_volume_annotations(pdb_file: Path, cavity_volume: float, pocket_volume: float, pore_volume: float, hub_volume: float):
    utils.set_resolution(2.0)

    pdb_structure = pdb.load_structure(pdb_file)
    prepared_structure, _ = volumizer.prepare_pdb_structure(pdb_structure)
    annotation, _ = volumizer.annotate_pdb_structure(prepared_structure)

    if cavity_volume is not None:
        assert annotation.cavity_volumes[0] == cavity_volume
    if pocket_volume is not None:
        assert annotation.pocket_volumes[0] == pocket_volume
    if pore_volume is not None:
        assert annotation.pore_volumes[0] == pore_volume
    if hub_volume is not None:
        assert annotation.hub_volumes[0] == hub_volume


@pytest.mark.parametrize(
    "pdb_file, removed_components, pore_volume, hub_volume",
    [
        (COMPONENT_TEST_PDB, set(), 38286.0, None),
        (COMPONENT_TEST_PDB, {"GLY", "ALA", "LEU", "GLU", "ASP", "LYS", "ILE"}, None, 35964.0)
    ]
)
def test_component_removal(pdb_file: Path, removed_components: set[str], pore_volume: float, hub_volume: float):
    utils.set_resolution(3.0)
    utils.reset_protein_components()
    utils.remove_protein_components(removed_components)
    pdb_structure = pdb.load_structure(pdb_file)
    prepared_structure, _ = volumizer.prepare_pdb_structure(pdb_structure)
    annotation, _ = volumizer.annotate_pdb_structure(prepared_structure)

    if pore_volume is not None:
        assert annotation.pore_volumes[0] == pore_volume
    if hub_volume is not None:
        assert annotation.hub_volumes[0] == hub_volume
