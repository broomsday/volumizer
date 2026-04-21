import pytest
from pathlib import Path
import tempfile
import inspect

import pandas as pd

from volumizer import volumizer, pdb, utils
from volumizer.paths import TEST_DIR

TEST_PDB_DIR = TEST_DIR / "pdbs"
TEST_DF_DIR = TEST_DIR / "dfs"
RCSB70_RUN_DIR = Path("data/runs/rcsb70/downloads")

GENERAL_TEST_PDB = TEST_PDB_DIR / "4jpn.pdb"
GENERAL_TEST_CIF = TEST_PDB_DIR / "4jpn.cif"
GENERAL_TEST_MMTF = TEST_PDB_DIR / "4jpn.mmtf"
GENERAL_TEST_DF = TEST_DF_DIR / "4jpn.json"
ASSEMBLY_TEST_CIF = TEST_PDB_DIR / "4jpp.cif"


def test_annotate_structure_volumes_defaults_min_voxels_to_four():
    signature = inspect.signature(volumizer.annotate_structure_volumes)
    assert signature.parameters["min_voxels"].default == 4


@pytest.mark.parametrize(
    "pdb_file, largest_type, largest_volume",
    [
        (
            TEST_PDB_DIR / "cavity.pdb",
            "cavity",
            184.0,
        ),
        (
            TEST_PDB_DIR / "pocket.pdb",
            "pocket",
            840.0,
        ),
        (
            TEST_PDB_DIR / "pore.pdb",
            "pore",
            1904.0,
        ),
        (
            TEST_PDB_DIR / "hub.pdb",
            "hub",
            2176.0,
        ),
        (
            ASSEMBLY_TEST_CIF,
            "pore",
            71160.0,
        ),
    ],
)
def test_volume_annotations(
    pdb_file: Path,
    largest_type: str,
    largest_volume: float,
):
    utils.set_resolution(2.0)

    pdb_structure = pdb.load_structure(pdb_file)
    prepared_structure = volumizer.prepare_pdb_structure(pdb_structure)
    annotation_df, _ = volumizer.annotate_structure_volumes(prepared_structure)

    assert (annotation_df.iloc[0].type == largest_type) and (
        annotation_df.iloc[0].volume == largest_volume
    )


@pytest.mark.parametrize(
    "pdb_file, expected_volume, expected_type",
    [
        (
            RCSB70_RUN_DIR / "5J3V.cif",
            19008.0,
            "pocket",
        ),
        (
            RCSB70_RUN_DIR / "9LJ9.cif",
            18117.0,
            "hub",
        ),
    ],
)
def test_reported_surface_regressions(
    pdb_file: Path,
    expected_volume: float,
    expected_type: str,
):
    utils.set_resolution(3.0)

    pdb_structure = pdb.load_structure(pdb_file)
    prepared_structure = volumizer.prepare_pdb_structure(pdb_structure)
    annotation_df, _ = volumizer.annotate_structure_volumes(prepared_structure)

    matching_rows = annotation_df[
        (annotation_df["type"] == expected_type)
        & (annotation_df["volume"] == expected_volume)
    ]
    assert len(matching_rows) >= 1


def test_9lj9_promotes_both_symmetric_large_volumes_to_hubs():
    utils.set_resolution(3.0)

    pdb_structure = pdb.load_structure(RCSB70_RUN_DIR / "9LJ9.cif")
    prepared_structure = volumizer.prepare_pdb_structure(pdb_structure)
    annotation_df, _ = volumizer.annotate_structure_volumes(prepared_structure)

    hub_volumes = set(annotation_df.loc[annotation_df["type"] == "hub", "volume"].tolist())
    assert 17091.0 in hub_volumes
    assert 18117.0 in hub_volumes


@pytest.mark.parametrize(
    "pdb_file, expected_annotation_df",
    [
        (
            GENERAL_TEST_PDB,
            pd.DataFrame.from_dict(
                {
                    "id": [0],
                    "type": ["pore"],
                    "volume": [38286.0],
                    "x": [108.221],
                    "y": [38.574],
                    "z": [36.310],
                    "cross_section_circularity": [0.8453],
                    "cross_section_uniformity": [0.8535],
                }
            ),
        ),
    ],
)
def test_volumize_pdb(pdb_file: Path, expected_annotation_df: pd.DataFrame):
    utils.set_resolution(3.0)

    annotation_df, _, _ = volumizer.volumize_pdb(pdb_file)
    assert annotation_df.iloc[0].to_dict() == expected_annotation_df.iloc[0].to_dict()


@pytest.mark.parametrize(
    "in_pdb_file, out_pdb_file, out_annotation_df",
    [
        (
            GENERAL_TEST_PDB,
            GENERAL_TEST_PDB.with_suffix(".annotated.pdb"),
            GENERAL_TEST_DF,
        ),
        (
            GENERAL_TEST_CIF,
            GENERAL_TEST_CIF.with_suffix(".annotated.pdb"),
            GENERAL_TEST_DF,
        ),
        (
            GENERAL_TEST_MMTF,
            GENERAL_TEST_MMTF.with_suffix(".annotated.pdb"),
            GENERAL_TEST_DF,
        ),
    ],
)
def test_volumize_pdb_and_save(
    in_pdb_file: Path, out_pdb_file: Path, out_annotation_df: Path
):
    temp_pdb = tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".pdb")
    temp_pdb_path = Path(temp_pdb.name)

    temp_df = tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8")
    temp_df_path = Path(temp_df.name)

    volumizer.volumize_pdb_and_save(
        in_pdb_file,
        temp_pdb_path,
        temp_df_path,
    )

    expected_annotated_pdb = pdb.load_structure(out_pdb_file)
    actual_annotated_pdb = pdb.load_structure(
        temp_pdb_path,
    )
    if not (expected_annotated_pdb == actual_annotated_pdb):
        raise AssertionError("Annotated PDB fixture mismatch")

    expected_annotation_df = pd.read_json(out_annotation_df)
    actual_annotation_df = pd.read_json(temp_df_path)
    assert (
        expected_annotation_df.iloc[0].to_dict()
        == actual_annotation_df.iloc[0].to_dict()
    )
