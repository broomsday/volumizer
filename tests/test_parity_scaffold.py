from volumizer import utils, volumizer
from volumizer.paths import TEST_DIR

from tests.parity_helpers import (
    assert_annotation_dataframes_close,
    load_annotation_json,
)


INPUT_PDB = TEST_DIR / "pdbs" / "4jpn.pdb"
GOLDEN_ANNOTATION = TEST_DIR / "golden" / "4jpn.annotation.json"


def test_golden_fixture_exists():
    assert GOLDEN_ANNOTATION.is_file()


def test_4jpn_annotation_matches_golden_output():
    utils.set_resolution(2.0)

    annotation_df, _, _ = volumizer.volumize_pdb(INPUT_PDB)
    expected_df = load_annotation_json(GOLDEN_ANNOTATION)

    assert_annotation_dataframes_close(
        annotation_df,
        expected_df,
        volume_tolerance=1e-6,
        dimension_tolerance=1e-3,
    )
