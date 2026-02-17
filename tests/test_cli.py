from pathlib import Path
from types import SimpleNamespace
import json

from volumizer import cli, rcsb
from volumizer.paths import TEST_DIR


TEST_PDB = TEST_DIR / "pdbs" / "cavity.pdb"


def test_resolve_input_structures_for_pdb_id(monkeypatch, tmp_path: Path):
    out_path = tmp_path / "1ABC.cif"
    out_path.write_text("dummy", encoding="utf-8")

    monkeypatch.setattr(
        rcsb,
        "download_structure_cif",
        lambda pdb_id, output_dir, overwrite=False, timeout=60.0: out_path,
    )

    args = SimpleNamespace(
        input=None,
        pdb_id="1abc",
        cluster_identity=None,
        overwrite=False,
        resume=False,
        timeout=60.0,
        max_structures=None,
        cluster_method=None,
        cluster_allow_all_methods=False,
        cluster_max_resolution=None,
        fail_fast=False,
    )
    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", out_path)]


def test_resolve_input_structures_for_cluster_identity(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False: [
            "1ABC",
            "2DEF",
        ],
    )

    metadata_by_id = {
        "1ABC": {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [1.5]},
        },
        "2DEF": {
            "exptl": [{"method": "ELECTRON MICROSCOPY"}],
            "rcsb_entry_info": {"resolution_combined": [2.8]},
        },
    }
    monkeypatch.setattr(
        rcsb,
        "fetch_entry_metadata",
        lambda pdb_id, timeout=60.0: metadata_by_id[pdb_id],
    )

    def _download(pdb_id, output_dir, overwrite=False, timeout=60.0):
        path = output_dir / f"{pdb_id}.cif"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("dummy", encoding="utf-8")
        return path

    monkeypatch.setattr(rcsb, "download_structure_cif", _download)

    args = SimpleNamespace(
        input=None,
        pdb_id=None,
        cluster_identity=30,
        overwrite=False,
        resume=False,
        timeout=60.0,
        max_structures=None,
        cluster_method=None,
        cluster_allow_all_methods=False,
        cluster_max_resolution=None,
        fail_fast=False,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert len(resolved) == 2
    assert resolved[0][0] == "1abc"
    assert resolved[1][0] == "2def"


def test_resolve_cluster_identity_default_method_filter_excludes_nmr(
    monkeypatch,
    tmp_path: Path,
):
    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False: [
            "1ABC",
            "2DEF",
        ],
    )

    metadata_by_id = {
        "1ABC": {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [1.9]},
        },
        "2DEF": {
            "exptl": [{"method": "SOLUTION NMR"}],
            "rcsb_entry_info": {},
        },
    }
    monkeypatch.setattr(
        rcsb,
        "fetch_entry_metadata",
        lambda pdb_id, timeout=60.0: metadata_by_id[pdb_id],
    )

    def _download(pdb_id, output_dir, overwrite=False, timeout=60.0):
        path = output_dir / f"{pdb_id}.cif"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("dummy", encoding="utf-8")
        return path

    monkeypatch.setattr(rcsb, "download_structure_cif", _download)

    args = SimpleNamespace(
        input=None,
        pdb_id=None,
        cluster_identity=30,
        overwrite=False,
        resume=False,
        timeout=60.0,
        max_structures=None,
        cluster_method=None,
        cluster_allow_all_methods=False,
        cluster_max_resolution=None,
        fail_fast=False,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]


def test_resolve_cluster_identity_max_resolution_filter(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False: [
            "1ABC",
            "2DEF",
        ],
    )

    metadata_by_id = {
        "1ABC": {
            "exptl": [{"method": "ELECTRON MICROSCOPY"}],
            "rcsb_entry_info": {"resolution_combined": [2.0]},
        },
        "2DEF": {
            "exptl": [{"method": "ELECTRON MICROSCOPY"}],
            "rcsb_entry_info": {"resolution_combined": [3.1]},
        },
    }
    monkeypatch.setattr(
        rcsb,
        "fetch_entry_metadata",
        lambda pdb_id, timeout=60.0: metadata_by_id[pdb_id],
    )

    def _download(pdb_id, output_dir, overwrite=False, timeout=60.0):
        path = output_dir / f"{pdb_id}.cif"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("dummy", encoding="utf-8")
        return path

    monkeypatch.setattr(rcsb, "download_structure_cif", _download)

    args = SimpleNamespace(
        input=None,
        pdb_id=None,
        cluster_identity=30,
        overwrite=False,
        resume=False,
        timeout=60.0,
        max_structures=None,
        cluster_method=["em"],
        cluster_allow_all_methods=False,
        cluster_max_resolution=2.5,
        fail_fast=False,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]


def test_analyze_structure_file_writes_cif_and_json(tmp_path: Path):
    result = cli.analyze_structure_file(
        source_label="cavity",
        input_path=TEST_PDB,
        output_dir=tmp_path,
        min_voxels=2,
        min_volume=None,
        overwrite=True,
    )

    structure_path = Path(result["structure_output"])
    annotation_path = Path(result["annotation_output"])

    assert structure_path.is_file()
    assert annotation_path.is_file()

    payload = json.loads(annotation_path.read_text(encoding="utf-8"))
    assert payload["source"] == "cavity"
    assert payload["num_volumes"] == len(payload["volumes"])
    assert payload["num_volumes"] >= 0
    assert isinstance(payload["volumes"], list)


def test_cli_main_single_input_writes_summary(monkeypatch, tmp_path: Path):
    monkeypatch.delenv("VOLUMIZER_BACKEND", raising=False)

    exit_code = cli.main(
        [
            "--input",
            str(TEST_PDB),
            "--output-dir",
            str(tmp_path),
            "--overwrite",
        ]
    )
    assert exit_code == 0

    summary_path = tmp_path / "run.summary.json"
    assert summary_path.is_file()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["num_processed"] == 1
    assert summary["num_failed"] == 0
    assert summary["num_skipped"] == 0


def test_cli_main_resume_skips_existing_output(monkeypatch, tmp_path: Path):
    monkeypatch.delenv("VOLUMIZER_BACKEND", raising=False)

    structure_path = tmp_path / "cavity.annotated.cif"
    annotation_path = tmp_path / "cavity.annotation.json"

    structure_path.write_text("data_dummy", encoding="utf-8")
    annotation_path.write_text(
        json.dumps(
            {
                "source": "cavity",
                "num_volumes": 3,
                "largest_type": "buried",
                "largest_volume": 123.0,
            }
        ),
        encoding="utf-8",
    )

    exit_code = cli.main(
        [
            "--input",
            str(TEST_PDB),
            "--output-dir",
            str(tmp_path),
            "--resume",
        ]
    )
    assert exit_code == 0

    summary_path = tmp_path / "run.summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["num_processed"] == 0
    assert summary["num_failed"] == 0
    assert summary["num_skipped"] == 1
    assert summary["skipped"][0]["source"] == "cavity"
    assert summary["skipped"][0]["num_volumes"] == 3
