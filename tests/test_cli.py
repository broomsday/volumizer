from pathlib import Path
from types import SimpleNamespace
import json

from volumizer import cli, rcsb
from volumizer.paths import TEST_DIR


TEST_PDB = TEST_DIR / "pdbs" / "cavity.pdb"


def _make_args(tmp_path: Path, **overrides) -> SimpleNamespace:
    defaults = {
        "command": "analyze",
        "cache_command": None,
        "input": None,
        "pdb_id": None,
        "manifest": None,
        "from_summary": None,
        "only": "failed",
        "cluster_identity": None,
        "output_dir": tmp_path,
        "download_dir": None,
        "max_structures": None,
        "num_shards": None,
        "shard_index": None,
        "write_manifest": None,
        "failures_manifest": None,
        "cluster_method": None,
        "cluster_allow_all_methods": False,
        "cluster_max_resolution": None,
        "metadata_cache": None,
        "no_metadata_cache": False,
        "checkpoint": None,
        "no_checkpoint": False,
        "progress_jsonl": None,
        "progress_interval": 30.0,
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


def test_normalize_argv_for_subcommands_infers_analyze():
    normalized = cli._normalize_argv_for_subcommands(
        ["--input", "input.cif", "--output-dir", "out"]
    )
    assert normalized[0] == "analyze"


def test_normalize_argv_for_subcommands_infers_cluster():
    normalized = cli._normalize_argv_for_subcommands(
        ["--cluster-identity", "30", "--output-dir", "out"]
    )
    assert normalized[0] == "cluster"


def test_resolve_input_structures_for_pdb_id(monkeypatch, tmp_path: Path):
    out_path = tmp_path / "1ABC.cif"
    out_path.write_text("dummy", encoding="utf-8")

    monkeypatch.setattr(
        rcsb,
        "download_structure_cif",
        lambda pdb_id, output_dir, overwrite=False, timeout=60.0, retries=0, retry_delay=1.0: out_path,
    )

    args = _make_args(tmp_path, command="analyze", pdb_id="1abc")
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

    args = _make_args(tmp_path, command="analyze", pdb_id="1abc", dry_run=True)
    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]


def test_resolve_input_structures_for_manifest_local_input(tmp_path: Path):
    manifest_path = tmp_path / "structures.manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "manifest_format": 1,
                "structures": [
                    {
                        "source": "local-cavity",
                        "input_path": str(TEST_PDB),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    args = _make_args(tmp_path, command="analyze", manifest=manifest_path)
    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("local-cavity", TEST_PDB)]


def test_resolve_input_structures_for_manifest_pdb_id_dry_run(
    monkeypatch,
    tmp_path: Path,
):
    manifest_path = tmp_path / "structures.manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "manifest_format": 1,
                "structures": [
                    {
                        "pdb_id": "1abc",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        rcsb,
        "download_structure_cif",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("should not download")),
    )

    args = _make_args(
        tmp_path,
        command="analyze",
        manifest=manifest_path,
        dry_run=True,
    )
    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]


def test_resolve_cluster_identity_writes_manifest(monkeypatch, tmp_path: Path):
    manifest_path = tmp_path / "cluster-selection.manifest.json"

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
            "rcsb_entry_info": {"resolution_combined": [1.8]},
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

    args = _make_args(
        tmp_path,
        command="cluster",
        cluster_identity=30,
        write_manifest=manifest_path,
        dry_run=True,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["manifest_format"] == 1
    assert payload["source_command"] == "cluster"
    assert payload["structures"] == [{"source": "1abc", "pdb_id": "1ABC"}]
    assert any(
        rejection.get("pdb_id") == "2DEF"
        and rejection.get("reason") == "experimental_method"
        for rejection in payload["rejections"]
    )


def test_resolve_input_structures_for_from_summary_failed(tmp_path: Path):
    previous_dir = tmp_path / "previous"
    downloads_dir = previous_dir / "downloads"
    downloads_dir.mkdir(parents=True)

    failed_path = downloads_dir / "1ABC.cif"
    failed_path.write_text("dummy", encoding="utf-8")
    ok_path = downloads_dir / "2DEF.cif"
    ok_path.write_text("dummy", encoding="utf-8")

    summary_path = previous_dir / "run.summary.json"
    summary_path.write_text(
        json.dumps(
            {
                "errors": [{"source": "1abc", "input_path": "downloads/1ABC.cif"}],
                "results": [{"source": "2def", "input_path": "downloads/2DEF.cif"}],
                "skipped": [],
                "planned": [],
            }
        ),
        encoding="utf-8",
    )

    args = _make_args(
        tmp_path,
        command="analyze",
        from_summary=summary_path,
        only="failed",
    )
    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", failed_path)]


def test_resolve_input_structures_for_from_summary_missing_input_raises(tmp_path: Path):
    summary_path = tmp_path / "run.summary.json"
    summary_path.write_text(
        json.dumps(
            {
                "errors": [{"source": "1abc", "input_path": "downloads/1ABC.cif"}],
                "results": [],
                "skipped": [],
                "planned": [],
            }
        ),
        encoding="utf-8",
    )

    args = _make_args(
        tmp_path,
        command="analyze",
        from_summary=summary_path,
        only="failed",
    )

    try:
        cli.resolve_input_structures(args, tmp_path, tmp_path)
        assert False, "expected FileNotFoundError"
    except FileNotFoundError as error:
        assert "Summary replay input structure does not exist" in str(error)


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

    args = _make_args(tmp_path, command="cluster", cluster_identity=30)

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert len(resolved) == 2
    assert resolved[0][0] == "1abc"
    assert resolved[1][0] == "2def"


def test_resolve_cluster_identity_applies_deterministic_shard(
    monkeypatch,
    tmp_path: Path,
):
    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False, retries=0, retry_delay=1.0: [
            "1ABC",
            "2DEF",
            "3GHI",
            "4JKL",
        ],
    )

    metadata_by_id = {
        "1ABC": {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.0]},
        },
        "2DEF": {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.0]},
        },
        "3GHI": {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.0]},
        },
        "4JKL": {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.0]},
        },
    }
    fetch_calls = []

    def _fetch(pdb_id, timeout=60.0, retries=0, retry_delay=1.0):
        fetch_calls.append(pdb_id)
        return metadata_by_id[pdb_id]

    monkeypatch.setattr(rcsb, "fetch_entry_metadata", _fetch)

    args = _make_args(
        tmp_path,
        command="cluster",
        cluster_identity=30,
        num_shards=2,
        shard_index=1,
        dry_run=True,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [
        ("2def", tmp_path / "2DEF.cif"),
        ("4jkl", tmp_path / "4JKL.cif"),
    ]
    assert fetch_calls == ["2DEF", "4JKL"]


def test_run_cli_cluster_shard_requires_both_flags(tmp_path: Path):
    args = _make_args(
        tmp_path,
        command="cluster",
        cluster_identity=30,
        num_shards=2,
        shard_index=None,
    )

    try:
        cli.run_cli(args)
        assert False, "expected ValueError"
    except ValueError as error:
        assert "Use --num-shards and --shard-index together." in str(error)


def test_run_cli_cluster_shard_index_bounds_checked(tmp_path: Path):
    args = _make_args(
        tmp_path,
        command="cluster",
        cluster_identity=30,
        num_shards=2,
        shard_index=2,
    )

    try:
        cli.run_cli(args)
        assert False, "expected ValueError"
    except ValueError as error:
        assert "--shard-index must be < --num-shards." in str(error)


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

    args = _make_args(tmp_path, command="cluster", cluster_identity=30)

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
        command="cluster",
        cluster_identity=30,
        cluster_method=["em"],
        cluster_max_resolution=2.5,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]


def test_resolve_cluster_identity_uses_metadata_cache(monkeypatch, tmp_path: Path):
    cache_path = tmp_path / "metadata-cache.json"
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
        command="cluster",
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
        command="cluster",
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
        command="cluster",
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
    args = _make_args(tmp_path, command="analyze", jobs=2)

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


def test_run_cli_parallel_jobs_emits_periodic_progress(
    monkeypatch,
    tmp_path: Path,
    capsys,
):
    args = _make_args(
        tmp_path,
        command="analyze",
        jobs=2,
        progress_interval=0.001,
    )

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

    stderr = capsys.readouterr().err
    assert "analysis progress:" in stderr
    assert "eta=" in stderr


def test_run_cli_progress_interval_zero_disables_periodic_progress(
    monkeypatch,
    tmp_path: Path,
    capsys,
):
    args = _make_args(
        tmp_path,
        command="analyze",
        jobs=1,
        progress_interval=0.0,
    )

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [
            ("single", tmp_path / "single.cif"),
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

    stderr = capsys.readouterr().err
    assert "analysis progress:" not in stderr


def test_run_cli_writes_failures_manifest(monkeypatch, tmp_path: Path):
    failures_manifest = tmp_path / "failed.manifest.json"
    args = _make_args(
        tmp_path,
        command="analyze",
        jobs=1,
        failures_manifest=failures_manifest,
    )

    first_input = tmp_path / "first.cif"
    second_input = tmp_path / "second.cif"

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [
            ("first", first_input),
            ("second", second_input),
        ],
    )

    def _analyze(source_label, input_path, output_dir, min_voxels, min_volume, overwrite):
        if source_label == "first":
            raise RuntimeError("boom")
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
    assert exit_code == 1
    assert failures_manifest.is_file()

    payload = json.loads(failures_manifest.read_text(encoding="utf-8"))
    assert payload["manifest_format"] == 1
    assert payload["source_command"] == "analyze"
    assert payload["num_failed"] == 1
    assert payload["structures"] == [
        {
            "source": "first",
            "input_path": str(first_input.resolve()),
        }
    ]

    summary = json.loads((tmp_path / "run.summary.json").read_text(encoding="utf-8"))
    assert summary["config"]["failures_manifest"] == str(failures_manifest)


def test_run_cli_dry_run_writes_plan_and_skips_analysis(monkeypatch, tmp_path: Path):
    args = _make_args(tmp_path, command="analyze", dry_run=True)

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


def test_run_cli_writes_progress_jsonl_events(monkeypatch, tmp_path: Path):
    progress_path = tmp_path / "run.progress.jsonl"
    args = _make_args(
        tmp_path,
        command="analyze",
        dry_run=True,
        progress_jsonl=progress_path,
    )

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [
            ("first", tmp_path / "first.cif"),
            ("second", tmp_path / "second.cif"),
        ],
    )

    exit_code = cli.run_cli(args)
    assert exit_code == 0
    assert progress_path.is_file()

    events = [
        json.loads(line)
        for line in progress_path.read_text(encoding="utf-8").splitlines()
        if len(line.strip()) > 0
    ]
    event_types = [event["event"] for event in events]

    assert "run_started" in event_types
    assert event_types.count("structure_planned") == 2
    assert event_types[-1] == "run_completed"


def test_run_cli_resume_uses_checkpoint_state(monkeypatch, tmp_path: Path):
    checkpoint_path = tmp_path / "state.checkpoint.json"

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [
            ("first", tmp_path / "first.cif"),
            ("second", tmp_path / "second.cif"),
        ],
    )

    def _analyze_first_pass(source_label, input_path, output_dir, min_voxels, min_volume, overwrite):
        if source_label == "first":
            structure_out = output_dir / "first.annotated.cif"
            annotation_out = output_dir / "first.annotation.json"
            structure_out.write_text("dummy", encoding="utf-8")
            annotation_out.write_text("{}", encoding="utf-8")
            return {
                "source": source_label,
                "input_path": str(input_path),
                "structure_output": str(structure_out),
                "annotation_output": str(annotation_out),
                "num_volumes": 1,
                "largest_type": "pore",
                "largest_volume": 1.0,
            }
        raise RuntimeError("boom")

    monkeypatch.setattr(cli, "analyze_structure_file", _analyze_first_pass)

    first_args = _make_args(
        tmp_path,
        command="analyze",
        checkpoint=checkpoint_path,
        jobs=1,
        fail_fast=True,
    )
    first_exit = cli.run_cli(first_args)
    assert first_exit == 1

    first_checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    assert len(first_checkpoint["results"]) == 1
    assert len(first_checkpoint["errors"]) == 1

    def _analyze_second_pass(source_label, input_path, output_dir, min_voxels, min_volume, overwrite):
        if source_label == "first":
            raise RuntimeError("first should be skipped on resume")
        structure_out = output_dir / "second.annotated.cif"
        annotation_out = output_dir / "second.annotation.json"
        structure_out.write_text("dummy", encoding="utf-8")
        annotation_out.write_text("{}", encoding="utf-8")
        return {
            "source": source_label,
            "input_path": str(input_path),
            "structure_output": str(structure_out),
            "annotation_output": str(annotation_out),
            "num_volumes": 2,
            "largest_type": "buried",
            "largest_volume": 2.0,
        }

    monkeypatch.setattr(cli, "analyze_structure_file", _analyze_second_pass)

    second_args = _make_args(
        tmp_path,
        command="analyze",
        checkpoint=checkpoint_path,
        jobs=1,
        resume=True,
    )
    second_exit = cli.run_cli(second_args)
    assert second_exit == 0

    summary = json.loads((tmp_path / "run.summary.json").read_text(encoding="utf-8"))
    assert summary["num_processed"] == 2
    assert summary["num_failed"] == 0


def test_cache_subcommand_inspect_and_clear_negative(tmp_path: Path):
    cache_path = tmp_path / "metadata-cache.json"
    cache_path.write_text(
        json.dumps(
            {
                "cache_format": 2,
                "entries": {"1ABC": {"dummy": 1}},
                "negative_entries": {"2DEF": {"status_code": 404}},
            }
        ),
        encoding="utf-8",
    )

    inspect_exit = cli.main(
        [
            "cache",
            "inspect",
            "--metadata-cache",
            str(cache_path),
        ]
    )
    assert inspect_exit == 0

    clear_exit = cli.main(
        [
            "cache",
            "clear-negative",
            "--metadata-cache",
            str(cache_path),
        ]
    )
    assert clear_exit == 0

    payload = json.loads(cache_path.read_text(encoding="utf-8"))
    assert payload["negative_entries"] == {}


def test_cli_main_analyze_subcommand_writes_summary(monkeypatch, tmp_path: Path):
    monkeypatch.delenv("VOLUMIZER_BACKEND", raising=False)

    exit_code = cli.main(
        [
            "analyze",
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


def test_cli_main_analyze_manifest_subcommand_writes_summary(monkeypatch, tmp_path: Path):
    monkeypatch.delenv("VOLUMIZER_BACKEND", raising=False)

    manifest_path = tmp_path / "structures.manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "manifest_format": 1,
                "structures": [
                    {
                        "source": "cavity-from-manifest",
                        "input_path": str(TEST_PDB),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    exit_code = cli.main(
        [
            "analyze",
            "--manifest",
            str(manifest_path),
            "--output-dir",
            str(tmp_path),
            "--overwrite",
        ]
    )
    assert exit_code == 0

    summary_path = tmp_path / "run.summary.json"
    assert summary_path.is_file()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["config"]["manifest"] == str(manifest_path)
    assert summary["num_processed"] == 1
    assert summary["num_failed"] == 0


def test_cli_main_analyze_from_summary_failed_writes_summary(monkeypatch, tmp_path: Path):
    monkeypatch.delenv("VOLUMIZER_BACKEND", raising=False)

    summary_path = tmp_path / "previous.summary.json"
    summary_path.write_text(
        json.dumps(
            {
                "errors": [{"source": "failed-cavity", "input_path": str(TEST_PDB)}],
                "results": [{"source": "ok-cavity", "input_path": str(TEST_PDB)}],
                "skipped": [],
                "planned": [],
            }
        ),
        encoding="utf-8",
    )

    exit_code = cli.main(
        [
            "analyze",
            "--from-summary",
            str(summary_path),
            "--only",
            "failed",
            "--output-dir",
            str(tmp_path),
            "--overwrite",
        ]
    )
    assert exit_code == 0

    summary = json.loads((tmp_path / "run.summary.json").read_text(encoding="utf-8"))
    assert summary["config"]["from_summary"] == str(summary_path)
    assert summary["config"]["from_summary_only"] == "failed"
    assert summary["num_processed"] == 1
    assert summary["num_failed"] == 0


def test_cli_main_cluster_subcommand_dry_run_writes_summary(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False, retries=0, retry_delay=1.0: [
            "1ABC"
        ],
    )
    monkeypatch.setattr(
        rcsb,
        "fetch_entry_metadata",
        lambda pdb_id, timeout=60.0, retries=0, retry_delay=1.0: {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.0]},
        },
    )

    exit_code = cli.main(
        [
            "cluster",
            "--cluster-identity",
            "30",
            "--output-dir",
            str(tmp_path),
            "--dry-run",
        ]
    )
    assert exit_code == 0

    summary_path = tmp_path / "run.summary.json"
    assert summary_path.is_file()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["config"]["command"] == "cluster"
    assert summary["num_processed"] == 0
    assert summary["num_planned"] == 1


def test_cli_main_legacy_single_input_writes_summary(monkeypatch, tmp_path: Path):
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
            "analyze",
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
