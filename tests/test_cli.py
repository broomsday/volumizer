from pathlib import Path
from types import SimpleNamespace
import json

from volumizer import cli, rcsb
from volumizer.paths import TEST_DIR


TEST_PDB = TEST_DIR / "pdbs" / "cavity.pdb"


def _make_args(tmp_path: Path, **overrides) -> SimpleNamespace:
    defaults = {
        "input": None,
        "pdb_id": None,
        "cluster_identity": None,
        "output_dir": tmp_path,
        "download_dir": None,
        "max_structures": None,
        "cluster_method": None,
        "cluster_allow_all_methods": False,
        "cluster_max_resolution": None,
        "metadata_cache": None,
        "no_metadata_cache": False,
        "resolution": 3.0,
        "min_voxels": 2,
        "min_volume": None,
        "backend": None,
        "keep_non_protein": False,
        "jobs": 1,
        "timeout": 60.0,
        "retries": 0,
        "retry_delay": 0.0,
        "overwrite": False,
        "resume": False,
        "dry_run": False,
        "fail_fast": False,
    }
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def test_resolve_input_structures_for_pdb_id(monkeypatch, tmp_path: Path):
    out_path = tmp_path / "1ABC.cif"
    out_path.write_text("dummy", encoding="utf-8")

    monkeypatch.setattr(
        rcsb,
        "download_structure_cif",
        lambda pdb_id, output_dir, overwrite=False, timeout=60.0, retries=0, retry_delay=1.0: out_path,
    )

    args = _make_args(tmp_path, pdb_id="1abc")
    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", out_path)]


def test_resolve_input_structures_for_pdb_id_dry_run_skips_download(
    monkeypatch,
    tmp_path: Path,
):
    monkeypatch.setattr(
        rcsb,
        "download_structure_cif",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("should not download")),
    )

    args = _make_args(tmp_path, pdb_id="1abc", dry_run=True)
    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]


def test_resolve_input_structures_for_cluster_identity(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False, retries=0, retry_delay=1.0: [
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
        lambda pdb_id, timeout=60.0, retries=0, retry_delay=1.0: metadata_by_id[pdb_id],
    )

    def _download(pdb_id, output_dir, overwrite=False, timeout=60.0, retries=0, retry_delay=1.0):
        path = output_dir / f"{pdb_id}.cif"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("dummy", encoding="utf-8")
        return path

    monkeypatch.setattr(rcsb, "download_structure_cif", _download)

    args = _make_args(tmp_path, cluster_identity=30)

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
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False, retries=0, retry_delay=1.0: [
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
        lambda pdb_id, timeout=60.0, retries=0, retry_delay=1.0: metadata_by_id[pdb_id],
    )

    def _download(pdb_id, output_dir, overwrite=False, timeout=60.0, retries=0, retry_delay=1.0):
        path = output_dir / f"{pdb_id}.cif"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("dummy", encoding="utf-8")
        return path

    monkeypatch.setattr(rcsb, "download_structure_cif", _download)

    args = _make_args(tmp_path, cluster_identity=30)

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]


def test_resolve_cluster_identity_max_resolution_filter(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False, retries=0, retry_delay=1.0: [
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
        lambda pdb_id, timeout=60.0, retries=0, retry_delay=1.0: metadata_by_id[pdb_id],
    )

    def _download(pdb_id, output_dir, overwrite=False, timeout=60.0, retries=0, retry_delay=1.0):
        path = output_dir / f"{pdb_id}.cif"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("dummy", encoding="utf-8")
        return path

    monkeypatch.setattr(rcsb, "download_structure_cif", _download)

    args = _make_args(
        tmp_path,
        cluster_identity=30,
        cluster_method=["em"],
        cluster_max_resolution=2.5,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]


def test_resolve_cluster_identity_uses_metadata_cache(monkeypatch, tmp_path: Path):
    cache_path = tmp_path / "metadata-cache.json"
    # Backward-compatible legacy format: positive entries at top level.
    cache_path.write_text(
        json.dumps(
            {
                "1ABC": {
                    "exptl": [{"method": "X-RAY DIFFRACTION"}],
                    "rcsb_entry_info": {"resolution_combined": [2.0]},
                }
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False, retries=0, retry_delay=1.0: [
            "1ABC",
            "2DEF",
        ],
    )

    fetch_calls = []

    def _fetch_entry_metadata(pdb_id, timeout=60.0, retries=0, retry_delay=1.0):
        fetch_calls.append(pdb_id)
        return {
            "exptl": [{"method": "ELECTRON MICROSCOPY"}],
            "rcsb_entry_info": {"resolution_combined": [2.6]},
        }

    monkeypatch.setattr(rcsb, "fetch_entry_metadata", _fetch_entry_metadata)

    def _download(pdb_id, output_dir, overwrite=False, timeout=60.0, retries=0, retry_delay=1.0):
        path = output_dir / f"{pdb_id}.cif"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("dummy", encoding="utf-8")
        return path

    monkeypatch.setattr(rcsb, "download_structure_cif", _download)

    args = _make_args(
        tmp_path,
        cluster_identity=30,
        metadata_cache=cache_path,
        jobs=2,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert len(resolved) == 2
    assert fetch_calls == ["2DEF"]

    updated_cache = json.loads(cache_path.read_text(encoding="utf-8"))
    assert updated_cache["cache_format"] == 2
    assert "1ABC" in updated_cache["entries"]
    assert "2DEF" in updated_cache["entries"]
    assert updated_cache["negative_entries"] == {}


def test_resolve_cluster_identity_uses_negative_metadata_cache(monkeypatch, tmp_path: Path):
    cache_path = tmp_path / "metadata-cache.json"
    cache_path.write_text(
        json.dumps(
            {
                "cache_format": 2,
                "entries": {
                    "1ABC": {
                        "exptl": [{"method": "X-RAY DIFFRACTION"}],
                        "rcsb_entry_info": {"resolution_combined": [2.0]},
                    }
                },
                "negative_entries": {
                    "2DEF": {
                        "reason": "permanent_metadata_error",
                        "status_code": 404,
                        "error": "not found",
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False, retries=0, retry_delay=1.0: [
            "1ABC",
            "2DEF",
        ],
    )

    fetch_calls = []
    monkeypatch.setattr(
        rcsb,
        "fetch_entry_metadata",
        lambda pdb_id, timeout=60.0, retries=0, retry_delay=1.0: fetch_calls.append(pdb_id),
    )

    def _download(pdb_id, output_dir, overwrite=False, timeout=60.0, retries=0, retry_delay=1.0):
        path = output_dir / f"{pdb_id}.cif"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("dummy", encoding="utf-8")
        return path

    monkeypatch.setattr(rcsb, "download_structure_cif", _download)

    args = _make_args(
        tmp_path,
        cluster_identity=30,
        metadata_cache=cache_path,
        max_structures=2,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]
    assert fetch_calls == []


def test_resolve_cluster_identity_updates_negative_cache_on_permanent_failure(
    monkeypatch,
    tmp_path: Path,
):
    cache_path = tmp_path / "metadata-cache.json"

    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False, retries=0, retry_delay=1.0: [
            "2DEF",
            "1ABC",
        ],
    )

    def _fetch_entry_metadata(pdb_id, timeout=60.0, retries=0, retry_delay=1.0):
        if pdb_id == "2DEF":
            raise rcsb.RCSBFetchError(
                "HTTP error while fetching https://example.org/2DEF: 404",
                url="https://example.org/2DEF",
                status_code=404,
                reason="not found",
                permanent=True,
            )
        return {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.2]},
        }

    monkeypatch.setattr(rcsb, "fetch_entry_metadata", _fetch_entry_metadata)

    def _download(pdb_id, output_dir, overwrite=False, timeout=60.0, retries=0, retry_delay=1.0):
        path = output_dir / f"{pdb_id}.cif"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("dummy", encoding="utf-8")
        return path

    monkeypatch.setattr(rcsb, "download_structure_cif", _download)

    args = _make_args(
        tmp_path,
        cluster_identity=30,
        metadata_cache=cache_path,
        max_structures=1,
        jobs=1,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]

    cache_payload = json.loads(cache_path.read_text(encoding="utf-8"))
    assert cache_payload["negative_entries"]["2DEF"]["status_code"] == 404


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


def test_run_cli_parallel_jobs_processes_multiple_structures(monkeypatch, tmp_path: Path):
    args = _make_args(tmp_path, jobs=2)

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [
            ("first", tmp_path / "first.cif"),
            ("second", tmp_path / "second.cif"),
        ],
    )

    def _analyze(source_label, input_path, output_dir, min_voxels, min_volume, overwrite):
        return {
            "source": source_label,
            "input_path": str(input_path),
            "structure_output": str(output_dir / f"{source_label}.annotated.cif"),
            "annotation_output": str(output_dir / f"{source_label}.annotation.json"),
            "num_volumes": 1,
            "largest_type": "pore",
            "largest_volume": 12.0,
        }

    monkeypatch.setattr(cli, "analyze_structure_file", _analyze)

    exit_code = cli.run_cli(args)
    assert exit_code == 0

    summary_path = tmp_path / "run.summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["num_processed"] == 2
    assert summary["num_failed"] == 0


def test_run_cli_dry_run_writes_plan_and_skips_analysis(monkeypatch, tmp_path: Path):
    args = _make_args(tmp_path, dry_run=True)

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [
            ("first", tmp_path / "first.cif"),
            ("second", tmp_path / "second.cif"),
        ],
    )
    monkeypatch.setattr(
        cli,
        "analyze_structure_file",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("should not analyze")),
    )

    exit_code = cli.run_cli(args)
    assert exit_code == 0

    summary_path = tmp_path / "run.summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["config"]["dry_run"] is True
    assert summary["num_planned"] == 2
    assert summary["num_processed"] == 0
    assert summary["num_failed"] == 0
    assert len(summary["planned"]) == 2


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
