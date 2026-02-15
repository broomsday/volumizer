import pytest

from volumizer import native_backend, pdb, utils, volumizer
from volumizer.paths import TEST_DIR
from tests.parity_helpers import assert_annotation_dataframes_close


@pytest.fixture(scope="module", autouse=True)
def _require_native_backend():
    native_backend.clear_backend_cache()
    if native_backend.get_native_module_for_mode("auto") is None:
        pytest.skip("native backend artifact/module not available")


INPUT_PDB = TEST_DIR / "pdbs" / "pore.pdb"


def _run_annotation():
    utils.set_resolution(2.0)
    native_backend.clear_backend_cache()

    structure = pdb.load_structure(INPUT_PDB)
    prepared_structure = volumizer.prepare_pdb_structure(structure)
    annotation_df, _ = volumizer.annotate_structure_volumes(prepared_structure)
    return annotation_df


def test_native_pipeline_matches_python(monkeypatch):
    monkeypatch.setenv("VOLUMIZER_BACKEND", "python")
    python_df = _run_annotation()

    monkeypatch.setenv("VOLUMIZER_BACKEND", "native")
    native_df = _run_annotation()

    assert_annotation_dataframes_close(
        native_df,
        python_df,
        volume_tolerance=1e-6,
        dimension_tolerance=1e-2,
    )
