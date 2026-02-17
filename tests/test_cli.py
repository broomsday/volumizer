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
        timeout=60.0,
        max_structures=None,
    )
    resolved = cli.resolve_input_structures(args, tmp_path)
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
        timeout=60.0,
        max_structures=None,
    )

    resolved = cli.resolve_input_structures(args, tmp_path)
    assert len(resolved) == 2
    assert resolved[0][0] == "1abc"
    assert resolved[1][0] == "2def"


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
