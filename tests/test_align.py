import numpy as np
import biotite.structure as bts
import pytest

from volumizer import align, pdb
from volumizer.paths import TEST_DIR


@pytest.mark.parametrize(
    "pdb_file",
    [
        TEST_DIR / "pdbs" / "cavity.pdb",
        TEST_DIR / "pdbs" / "4jpn.pdb",
        TEST_DIR / "pdbs" / "4jpp.cif",
    ],
)
def test_align_structure_matches_biotite_orient_principal_components(pdb_file):
    structure = pdb.clean_structure(pdb.load_structure(pdb_file))

    expected = bts.orient_principal_components(structure)
    actual = align.align_structure(structure)

    assert np.array_equal(actual.coord, expected.coord)
    assert expected == actual
