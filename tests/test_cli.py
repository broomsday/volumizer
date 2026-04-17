from pathlib import Path
from types import SimpleNamespace
import json
import time

import pandas as pd

from volumizer import cli, rcsb
from volumizer.paths import TEST_DIR


TEST_PDB = TEST_DIR / "pdbs" / "cavity.pdb"


class _DummyOutputStructure:
    def __init__(self, res_name: list[str]):
        self.res_name = list(res_name)

    def __len__(self) -> int:
        return len(self.res_name)

    def __getitem__(self, indices):
        return _DummyOutputStructure([self.res_name[index] for index in indices])

    def __add__(self, other):
        return _DummyOutputStructure([*self.res_name, *other.res_name])


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
        "cluster_max_residues": 20000,
        "metadata_cache": None,
        "no_metadata_cache": False,
        "checkpoint": None,
        "no_checkpoint": False,
        "progress_jsonl": None,
        "progress_interval": 30.0,
        "resolution": 3.0,
        "min_voxels": 4,
        "min_volume": None,
        "include_hubs": False,
        "backend": None,
        "surface_connectivity": "custom18",
        "merge_mouth_gap_voxels": 1,
        "assembly_policy": "biological",
        "keep_non_protein": False,
        "jobs": 1,
        "timeout": 60.0,
        "worker_timeout_seconds": 1800.0,
        "retries": 0,
        "retry_delay": 0.0,
        "overwrite": False,
        "reannotate": False,
        "resume": False,
        "dry_run": False,
        "fail_fast": False,
    }
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def _patch_cluster_fetch(
    monkeypatch,
    representative_to_members: dict[str, list[str]],
) -> None:
    representative_ids = list(representative_to_members.keys())

    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, include_non_pdb=False, retries=0, retry_delay=1.0: representative_ids,
    )
    monkeypatch.setattr(
        rcsb,
        "fetch_cluster_representative_member_entry_ids",
        lambda identity, max_structures=None, timeout=60.0, retries=0, retry_delay=1.0: {
            representative_id: list(member_ids)
            for representative_id, member_ids in representative_to_members.items()
        },
    )


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



def test_normalize_argv_for_subcommands_preserves_version_flag():
    assert cli._normalize_argv_for_subcommands(["--version"]) == ["--version"]
    assert cli._normalize_argv_for_subcommands(["-V"]) == ["-V"]


def test_main_version_flag_outputs_version(capsys):
    try:
        cli.main(["--version"])
        assert False, "expected SystemExit for --version"
    except SystemExit as error:
        assert error.code == 0

    stdout = capsys.readouterr().out.strip()
    assert stdout.startswith("volumizer ")


def test_build_parser_defaults_min_voxels_to_four():
    args = cli.build_parser().parse_args(
        ["analyze", "--input", "input.cif", "--output-dir", "out"]
    )
    assert args.min_voxels == 4


def test_build_parser_defaults_surface_connectivity_to_custom18():
    args = cli.build_parser().parse_args(
        ["analyze", "--input", "input.cif", "--output-dir", "out"]
    )
    assert args.surface_connectivity == "custom18"


def test_build_parser_defaults_merge_mouth_gap_voxels_to_one():
    args = cli.build_parser().parse_args(
        ["analyze", "--input", "input.cif", "--output-dir", "out"]
    )
    assert args.merge_mouth_gap_voxels == 1


def test_build_parser_defaults_cluster_max_residues_to_twenty_thousand():
    args = cli.build_parser().parse_args(
        ["cluster", "--cluster-identity", "30", "--output-dir", "out"]
    )
    assert args.cluster_max_residues == 20000


def test_build_analysis_worker_command_include_hubs_flag():
    command = cli._build_analysis_worker_command(
        source_label="sample",
        input_path=Path("input.cif"),
        output_dir=Path("out"),
        min_voxels=2,
        min_volume=None,
        overwrite=False,
        assembly_policy="biological",
        resolution=3.0,
        keep_non_protein=False,
        backend=None,
        surface_connectivity="custom18",
        merge_mouth_gap_voxels=1,
        max_residues=None,
        include_hubs=False,
    )
    assert "--include-hubs" not in command
    assert command[command.index("--surface-connectivity") + 1] == "custom18"
    assert command[command.index("--merge-mouth-gap-voxels") + 1] == "1"

    include_hubs_command = cli._build_analysis_worker_command(
        source_label="sample",
        input_path=Path("input.cif"),
        output_dir=Path("out"),
        min_voxels=2,
        min_volume=None,
        overwrite=False,
        assembly_policy="biological",
        resolution=3.0,
        keep_non_protein=False,
        backend=None,
        surface_connectivity="18",
        merge_mouth_gap_voxels=-1,
        max_residues=None,
        include_hubs=True,
    )
    assert "--include-hubs" in include_hubs_command
    assert include_hubs_command[
        include_hubs_command.index("--surface-connectivity") + 1
    ] == "18"
    assert include_hubs_command[
        include_hubs_command.index("--merge-mouth-gap-voxels") + 1
    ] == "-1"


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

    _patch_cluster_fetch(
        monkeypatch,
        {
            "1ABC": ["1ABC", "1ABD"],
            "2DEF": ["2DEF", "2DEG"],
        },
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
    assert payload["structures"] == [
        {
            "source": "1abc",
            "pdb_id": "1ABC",
            "cluster_member_pdb_ids": ["1ABC", "1ABD"],
        }
    ]
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
    _patch_cluster_fetch(
        monkeypatch,
        {
            "1ABC": ["1ABC", "1ABD"],
            "2DEF": ["2DEF", "2DEG"],
        },
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
    _patch_cluster_fetch(
        monkeypatch,
        {
            "1ABC": ["1ABC"],
            "2DEF": ["2DEF"],
            "3GHI": ["3GHI"],
            "4JKL": ["4JKL"],
        },
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
    _patch_cluster_fetch(
        monkeypatch,
        {
            "1ABC": ["1ABC"],
            "2DEF": ["2DEF"],
        },
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
    _patch_cluster_fetch(
        monkeypatch,
        {
            "1ABC": ["1ABC"],
            "2DEF": ["2DEF"],
        },
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


def test_resolve_cluster_identity_uses_metadata_store(monkeypatch, tmp_path: Path):
    metadata_store = tmp_path / "entry_metadata"
    cli._write_metadata_record(
        metadata_store,
        "1ABC",
        entry_metadata={
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.0]},
        },
    )

    _patch_cluster_fetch(
        monkeypatch,
        {
            "1ABC": ["1ABC"],
            "2DEF": ["2DEF"],
        },
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
        metadata_cache=metadata_store,
        jobs=2,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert len(resolved) == 2
    assert fetch_calls == ["2DEF"]

    cached_entry, cached_negative = cli._load_metadata_record(metadata_store, "2DEF")
    assert cached_negative is None
    assert cached_entry is not None
    assert cached_entry["rcsb_entry_info"]["resolution_combined"] == [2.6]


def test_resolve_cluster_identity_uses_negative_metadata_store(monkeypatch, tmp_path: Path):
    metadata_store = tmp_path / "entry_metadata"
    cli._write_metadata_record(
        metadata_store,
        "1ABC",
        entry_metadata={
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.0]},
        },
    )
    cli._write_metadata_record(
        metadata_store,
        "2DEF",
        negative_entry={
            "reason": "permanent_metadata_error",
            "status_code": 404,
            "error": "not found",
        },
    )

    _patch_cluster_fetch(
        monkeypatch,
        {
            "1ABC": ["1ABC"],
            "2DEF": ["2DEF"],
        },
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
        metadata_cache=metadata_store,
        max_structures=2,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]
    assert fetch_calls == []


def test_resolve_cluster_identity_updates_negative_metadata_store_on_permanent_failure(
    monkeypatch,
    tmp_path: Path,
):
    metadata_store = tmp_path / "entry_metadata"

    _patch_cluster_fetch(
        monkeypatch,
        {
            "2DEF": ["2DEF"],
            "1ABC": ["1ABC"],
        },
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
        metadata_cache=metadata_store,
        max_structures=1,
        jobs=1,
    )

    resolved = cli.resolve_input_structures(args, tmp_path, tmp_path)
    assert resolved == [("1abc", tmp_path / "1ABC.cif")]

    cached_entry, cached_negative = cli._load_metadata_record(metadata_store, "2DEF")
    assert cached_entry is None
    assert cached_negative is not None
    assert cached_negative["status_code"] == 404


def test_load_metadata_cache_scans_metadata_store_directory(tmp_path: Path):
    metadata_store = tmp_path / "entry_metadata"
    cli._write_metadata_record(
        metadata_store,
        "1ABC",
        entry_metadata={
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.0]},
        },
    )
    cli._write_metadata_record(
        metadata_store,
        "2DEF",
        entry_metadata={
            "exptl": [{"method": "ELECTRON MICROSCOPY"}],
            "rcsb_entry_info": {"resolution_combined": [3.1]},
        },
    )
    cli._write_metadata_record(
        metadata_store,
        "3GHI",
        negative_entry={
            "reason": "permanent_metadata_error",
            "status_code": 404,
            "error": "not found",
        },
    )

    entries, negative_entries = cli._load_metadata_cache(metadata_store)
    assert "1ABC" in entries
    assert "2DEF" in entries
    assert negative_entries["3GHI"]["status_code"] == 404


def test_load_metadata_record_recovers_non_utf8_bytes(tmp_path: Path):
    metadata_store = tmp_path / "entry_metadata"
    metadata_store.mkdir()
    record_path = metadata_store / "2DEF.json"
    good_json = json.dumps(
        {
            "metadata_record_format": 1,
            "pdb_id": "2DEF",
            "kind": "entry",
            "entry_metadata": {
                "exptl": [{"method": "ELECTRON MICROSCOPY"}],
                "rcsb_entry_info": {"resolution_combined": [3.1]},
            },
        }
    )
    raw_bytes = good_json.encode("utf-8")
    idx = raw_bytes.index(b"ELECTRON")
    corrupted = raw_bytes[:idx] + b"\x91" + raw_bytes[idx + 1 :]
    record_path.write_bytes(corrupted)

    entry_metadata, negative_entry = cli._load_metadata_record(metadata_store, "2DEF")
    assert negative_entry is None
    assert entry_metadata is not None
    assert entry_metadata["rcsb_entry_info"]["resolution_combined"] == [3.1]


def test_load_metadata_cache_returns_empty_for_non_directory(tmp_path: Path):
    cache_path = tmp_path / "metadata-cache.json"
    cache_path.write_text("{}", encoding="utf-8")
    entries, negative_entries = cli._load_metadata_cache(cache_path)
    assert entries == {}
    assert negative_entries == {}


def test_download_cluster_structures_reuses_existing_downloads_without_resume(
    monkeypatch,
    tmp_path: Path,
):
    existing_path = tmp_path / "1ABC.cif"
    existing_path.write_text("dummy", encoding="utf-8")

    monkeypatch.setattr(
        rcsb,
        "download_structure_cif",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("existing download should not be fetched again")
        ),
    )

    downloaded = cli._download_cluster_structures(
        ["1ABC"],
        _make_args(tmp_path, command="cluster", cluster_identity=30),
        tmp_path,
        tmp_path,
    )

    assert downloaded == {"1ABC": existing_path}


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


def test_analyze_structure_file_excludes_hubs_from_written_outputs_by_default(
    monkeypatch,
    tmp_path: Path,
):
    input_path = tmp_path / "sample.cif"
    input_path.write_text("dummy", encoding="utf-8")
    saved_structure = {}

    def _save_structure(structure, output_path):
        saved_structure["res_name"] = list(structure.res_name)
        Path(output_path).write_text("filtered-structure", encoding="utf-8")

    dummy_pdb = SimpleNamespace(
        load_structure=lambda input_path, assembly_policy="biological": "input",
        save_structure=_save_structure,
        ensure_b_factor_annotation=lambda structure: None,
        compute_sse_fractions=lambda prepared_structure: {
            "frac_alpha": 0.1,
            "frac_beta": 0.2,
            "frac_coil": 0.7,
        },
    )
    dummy_volumizer = SimpleNamespace(
        prepare_pdb_structure=lambda structure: _DummyOutputStructure(["PRO"]),
        annotate_structure_volumes=lambda prepared_structure, min_voxels=2, min_volume=None: (
            pd.DataFrame(
                [
                    {"id": 0, "type": "hub", "volume": 50.0},
                    {"id": 1, "type": "pore", "volume": 12.0},
                ]
            ),
            _DummyOutputStructure(["HUB", "POR"]),
        ),
    )
    dummy_utils = SimpleNamespace(
        get_active_backend=lambda: "native",
        VOXEL_SIZE=3.0,
    )

    monkeypatch.setattr(cli, "_pdb_module", lambda: dummy_pdb)
    monkeypatch.setattr(cli, "_volumizer_module", lambda: dummy_volumizer)
    monkeypatch.setattr(cli, "_utils_module", lambda: dummy_utils)

    result = cli.analyze_structure_file(
        source_label="sample",
        input_path=input_path,
        output_dir=tmp_path,
        min_voxels=2,
        min_volume=None,
        overwrite=False,
        assembly_policy="biological",
    )

    payload = json.loads((tmp_path / "sample.annotation.json").read_text(encoding="utf-8"))
    assert result["num_volumes"] == 1
    assert result["largest_type"] == "pore"
    assert payload["num_volumes"] == 1
    assert payload["largest_type"] == "pore"
    assert [volume["type"] for volume in payload["volumes"]] == ["pore"]
    assert saved_structure["res_name"] == ["PRO", "POR"]


def test_analyze_structure_file_include_hubs_preserves_written_outputs(
    monkeypatch,
    tmp_path: Path,
):
    input_path = tmp_path / "sample.cif"
    input_path.write_text("dummy", encoding="utf-8")
    saved_structure = {}

    def _save_structure(structure, output_path):
        saved_structure["res_name"] = list(structure.res_name)
        Path(output_path).write_text("filtered-structure", encoding="utf-8")

    dummy_pdb = SimpleNamespace(
        load_structure=lambda input_path, assembly_policy="biological": "input",
        save_structure=_save_structure,
        ensure_b_factor_annotation=lambda structure: None,
        compute_sse_fractions=lambda prepared_structure: {
            "frac_alpha": 0.1,
            "frac_beta": 0.2,
            "frac_coil": 0.7,
        },
    )
    dummy_volumizer = SimpleNamespace(
        prepare_pdb_structure=lambda structure: _DummyOutputStructure(["PRO"]),
        annotate_structure_volumes=lambda prepared_structure, min_voxels=2, min_volume=None: (
            pd.DataFrame(
                [
                    {"id": 0, "type": "hub", "volume": 50.0},
                    {"id": 1, "type": "pore", "volume": 12.0},
                ]
            ),
            _DummyOutputStructure(["HUB", "POR"]),
        ),
    )
    dummy_utils = SimpleNamespace(
        get_active_backend=lambda: "native",
        VOXEL_SIZE=3.0,
    )

    monkeypatch.setattr(cli, "_pdb_module", lambda: dummy_pdb)
    monkeypatch.setattr(cli, "_volumizer_module", lambda: dummy_volumizer)
    monkeypatch.setattr(cli, "_utils_module", lambda: dummy_utils)

    result = cli.analyze_structure_file(
        source_label="sample",
        input_path=input_path,
        output_dir=tmp_path,
        min_voxels=2,
        min_volume=None,
        overwrite=False,
        assembly_policy="biological",
        include_hubs=True,
    )

    payload = json.loads((tmp_path / "sample.annotation.json").read_text(encoding="utf-8"))
    assert result["num_volumes"] == 2
    assert result["largest_type"] == "hub"
    assert payload["num_volumes"] == 2
    assert payload["largest_type"] == "hub"
    assert [volume["type"] for volume in payload["volumes"]] == ["hub", "pore"]
    assert saved_structure["res_name"] == ["PRO", "HUB", "POR"]


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

    def _analyze(source_label, input_path, output_dir, min_voxels, min_volume, overwrite, assembly_policy="biological"):
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


def test_run_cli_parallel_jobs_batches_checkpoint_persistence_after_completion(
    monkeypatch,
    tmp_path: Path,
):
    args = _make_args(tmp_path, command="analyze", jobs=2)

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [
            ("first", tmp_path / "first.cif"),
            ("second", tmp_path / "second.cif"),
            ("third", tmp_path / "third.cif"),
        ],
    )

    persist_snapshots = []

    def _persist_checkpoint(self):
        persist_snapshots.append(
            {
                "results": len(self.results),
                "errors": len(self.errors),
                "skipped": len(self.skipped),
                "planned": len(self.planned),
            }
        )

    monkeypatch.setattr(cli._RunTracker, "persist_checkpoint", _persist_checkpoint)

    def _analyze(source_label, input_path, output_dir, min_voxels, min_volume, overwrite, assembly_policy="biological"):
        if source_label == "first":
            return {
                "source": source_label,
                "input_path": str(input_path),
                "structure_output": str(output_dir / f"{source_label}.annotated.cif"),
                "annotation_output": str(output_dir / f"{source_label}.annotation.json"),
                "num_volumes": 1,
                "largest_type": "pore",
                "largest_volume": 12.0,
            }
        if source_label == "second":
            raise RuntimeError("boom")
        raise cli.PostAssemblyResidueLimitExceeded(
            actual_residues=12000,
            max_residues=10000,
            assembly_policy=assembly_policy,
        )

    monkeypatch.setattr(cli, "analyze_structure_file", _analyze)

    exit_code = cli.run_cli(args)
    assert exit_code == 1

    summary = json.loads((tmp_path / "run.summary.json").read_text(encoding="utf-8"))
    assert summary["num_processed"] == 1
    assert summary["num_failed"] == 1
    assert summary["num_skipped"] == 1

    assert persist_snapshots == [
        {"results": 0, "errors": 0, "skipped": 0, "planned": 0},
        {"results": 1, "errors": 1, "skipped": 1, "planned": 0},
        {"results": 1, "errors": 1, "skipped": 1, "planned": 0},
    ]


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

    def _analyze(source_label, input_path, output_dir, min_voxels, min_volume, overwrite, assembly_policy="biological"):
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
    assert "active=" in stderr
    assert "queued=" in stderr


def test_run_cli_parallel_jobs_emits_heartbeat_before_first_completion(
    monkeypatch,
    tmp_path: Path,
    capsys,
):
    args = _make_args(
        tmp_path,
        command="analyze",
        jobs=2,
        progress_interval=0.01,
    )

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [
            ("first", tmp_path / "first.cif"),
            ("second", tmp_path / "second.cif"),
        ],
    )

    def _analyze(source_label, input_path, output_dir, min_voxels, min_volume, overwrite, assembly_policy="biological"):
        time.sleep(0.2)
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
    assert stderr.count("analysis progress:") >= 2


def test_run_cli_cluster_native_uses_isolated_workers(
    monkeypatch,
    tmp_path: Path,
    capsys,
):
    args = _make_args(
        tmp_path,
        command="cluster",
        cluster_identity=30,
        jobs=2,
    )

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir, metadata_cache_path=None: [
            ("first", tmp_path / "first.cif"),
            ("second", tmp_path / "second.cif"),
        ],
    )

    monkeypatch.setattr(
        cli,
        "_native_backend_module",
        lambda: SimpleNamespace(
            BACKEND_ENV=cli.BACKEND_ENV,
            clear_backend_cache=lambda: None,
            active_backend=lambda: "native",
        ),
    )
    monkeypatch.setattr(
        cli,
        "_utils_module",
        lambda: SimpleNamespace(
            set_resolution=lambda resolution: None,
            set_non_protein=lambda keep_non_protein: None,
            set_surface_component_connectivity_mode=lambda mode: None,
            set_surface_mouth_merge_gap_voxels=lambda gap: None,
        ),
    )
    monkeypatch.setattr(cli, "analyze_structure_file", lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("cluster native path should use isolated workers")))

    def _worker(**kwargs):
        source_label = kwargs["source_label"]
        return {
            "source": source_label,
            "input_path": str(kwargs["input_path"]),
            "structure_output": str(kwargs["output_dir"] / f"{source_label}.annotated.cif"),
            "annotation_output": str(kwargs["output_dir"] / f"{source_label}.annotation.json"),
            "num_volumes": 1,
            "largest_type": "pore",
            "largest_volume": 12.0,
        }

    monkeypatch.setattr(cli, "_run_isolated_analysis_worker", _worker)

    exit_code = cli.run_cli(args)
    assert exit_code == 0

    summary = json.loads((tmp_path / "run.summary.json").read_text(encoding="utf-8"))
    assert summary["num_processed"] == 2
    assert summary["num_failed"] == 0

    stderr = capsys.readouterr().err
    assert "analysis worker mode: isolated subprocesses" in stderr


def test_run_cli_cluster_native_records_isolated_worker_failure(
    monkeypatch,
    tmp_path: Path,
):
    args = _make_args(
        tmp_path,
        command="cluster",
        cluster_identity=30,
        jobs=2,
    )

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir, metadata_cache_path=None: [
            ("first", tmp_path / "first.cif"),
            ("second", tmp_path / "second.cif"),
        ],
    )
    monkeypatch.setattr(
        cli,
        "_native_backend_module",
        lambda: SimpleNamespace(
            BACKEND_ENV=cli.BACKEND_ENV,
            clear_backend_cache=lambda: None,
            active_backend=lambda: "native",
        ),
    )
    monkeypatch.setattr(
        cli,
        "_utils_module",
        lambda: SimpleNamespace(
            set_resolution=lambda resolution: None,
            set_non_protein=lambda keep_non_protein: None,
            set_surface_component_connectivity_mode=lambda mode: None,
            set_surface_mouth_merge_gap_voxels=lambda gap: None,
        ),
    )

    def _worker(**kwargs):
        if kwargs["source_label"] == "first":
            raise RuntimeError("analysis worker terminated by signal 11 (SIGSEGV)")
        return {
            "source": kwargs["source_label"],
            "input_path": str(kwargs["input_path"]),
            "structure_output": str(kwargs["output_dir"] / "second.annotated.cif"),
            "annotation_output": str(kwargs["output_dir"] / "second.annotation.json"),
            "num_volumes": 1,
            "largest_type": "pore",
            "largest_volume": 12.0,
        }

    monkeypatch.setattr(cli, "_run_isolated_analysis_worker", _worker)

    exit_code = cli.run_cli(args)
    assert exit_code == 1

    summary = json.loads((tmp_path / "run.summary.json").read_text(encoding="utf-8"))
    assert summary["num_processed"] == 1
    assert summary["num_failed"] == 1
    assert "SIGSEGV" in summary["errors"][0]["error"]


def test_run_cli_cluster_native_records_post_assembly_residue_limit_skip(
    monkeypatch,
    tmp_path: Path,
    capsys,
):
    args = _make_args(
        tmp_path,
        command="cluster",
        cluster_identity=30,
        jobs=2,
    )

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir, metadata_cache_path=None: [
            ("first", tmp_path / "first.cif"),
            ("second", tmp_path / "second.cif"),
        ],
    )
    monkeypatch.setattr(
        cli,
        "_native_backend_module",
        lambda: SimpleNamespace(
            BACKEND_ENV=cli.BACKEND_ENV,
            clear_backend_cache=lambda: None,
            active_backend=lambda: "native",
        ),
    )
    monkeypatch.setattr(
        cli,
        "_utils_module",
        lambda: SimpleNamespace(
            set_resolution=lambda resolution: None,
            set_non_protein=lambda keep_non_protein: None,
            set_surface_component_connectivity_mode=lambda mode: None,
            set_surface_mouth_merge_gap_voxels=lambda gap: None,
        ),
    )

    def _worker(**kwargs):
        if kwargs["source_label"] == "first":
            raise cli.PostAssemblyResidueLimitExceeded(
                actual_residues=12000,
                max_residues=10000,
                assembly_policy="biological",
            )
        return {
            "source": kwargs["source_label"],
            "input_path": str(kwargs["input_path"]),
            "structure_output": str(kwargs["output_dir"] / "second.annotated.cif"),
            "annotation_output": str(kwargs["output_dir"] / "second.annotation.json"),
            "num_volumes": 1,
            "largest_type": "pore",
            "largest_volume": 12.0,
        }

    monkeypatch.setattr(cli, "_run_isolated_analysis_worker", _worker)

    exit_code = cli.run_cli(args)
    assert exit_code == 0

    summary = json.loads((tmp_path / "run.summary.json").read_text(encoding="utf-8"))
    assert summary["num_processed"] == 1
    assert summary["num_failed"] == 0
    assert summary["num_skipped"] == 1
    assert summary["skipped"][0]["reason"] == "post_assembly_residue_limit"
    assert summary["skipped"][0]["actual_residues"] == 12000
    assert summary["skipped"][0]["max_residues"] == 10000
    status_payload = json.loads(
        (tmp_path / "first.status.json").read_text(encoding="utf-8")
    )
    assert status_payload["kind"] == "terminal_skip"
    assert status_payload["reason"] == "post_assembly_residue_limit"
    assert status_payload["actual_residues"] == 12000
    assert status_payload["max_residues"] == 10000

    stderr = capsys.readouterr().err
    assert "skipping first:" in stderr


def test_run_isolated_analysis_worker_reports_signal(monkeypatch, tmp_path: Path):
    calls = []

    monkeypatch.setattr(
        cli.subprocess,
        "run",
        lambda *args, **kwargs: (
            calls.append((args, kwargs))
            or SimpleNamespace(
                returncode=-11,
                stdout="",
                stderr="",
            )
        ),
    )

    try:
        cli._run_isolated_analysis_worker(
            source_label="boom",
            input_path=tmp_path / "boom.cif",
            output_dir=tmp_path,
            min_voxels=2,
            min_volume=None,
            overwrite=False,
            assembly_policy="biological",
            resolution=3.0,
            keep_non_protein=False,
            backend="native",
            surface_connectivity="custom18",
            merge_mouth_gap_voxels=1,
        )
        assert False, "expected RuntimeError"
    except RuntimeError as error:
        assert "signal 11 (SIGSEGV)" in str(error)
        assert "attempt 1 backend=native:" in str(error)
        assert "attempt 3 backend=python:" in str(error)

    assert len(calls) == 3
    env = calls[0][1]["env"]
    assert env["PYTHONFAULTHANDLER"] == "1"
    for variable in cli._ANALYSIS_WORKER_THREAD_ENV_VARS:
        assert env[variable] == "1"


def test_run_isolated_analysis_worker_retries_native_signal_then_succeeds(
    monkeypatch,
    tmp_path: Path,
):
    responses = [
        SimpleNamespace(returncode=-11, stdout="", stderr=""),
        SimpleNamespace(
            returncode=0,
            stdout=json.dumps(
                {
                    "source": "boom",
                    "pdb_id": "1ABC",
                    "input_path": str(tmp_path / "boom.cif"),
                    "structure_output": str(tmp_path / "boom.annotated.cif"),
                    "annotation_output": str(tmp_path / "boom.annotation.json"),
                    "num_volumes": 1,
                    "largest_type": "pore",
                    "largest_volume": 12.0,
                }
            ),
            stderr="",
        ),
    ]
    calls = []

    def _run(*args, **kwargs):
        calls.append((args, kwargs))
        return responses.pop(0)

    monkeypatch.setattr(cli.subprocess, "run", _run)

    result = cli._run_isolated_analysis_worker(
        source_label="boom",
        input_path=tmp_path / "boom.cif",
        output_dir=tmp_path,
        min_voxels=2,
        min_volume=None,
        overwrite=False,
        assembly_policy="biological",
        resolution=3.0,
        keep_non_protein=False,
        backend="native",
        surface_connectivity="custom18",
        merge_mouth_gap_voxels=1,
        max_residues=10000,
    )

    assert result["source"] == "boom"
    assert len(calls) == 2
    first_command = calls[0][0][0]
    second_command = calls[1][0][0]
    assert first_command[first_command.index("--backend") + 1] == "native"
    assert second_command[second_command.index("--backend") + 1] == "native"
    assert "--overwrite" not in first_command
    assert "--overwrite" in second_command


def test_run_isolated_analysis_worker_retries_native_timeout_then_succeeds(
    monkeypatch,
    tmp_path: Path,
):
    responses = [
        cli.subprocess.TimeoutExpired(
            cmd=["python", "-m", "volumizer._analysis_worker"],
            timeout=12.5,
        ),
        SimpleNamespace(
            returncode=0,
            stdout=json.dumps(
                {
                    "source": "boom",
                    "pdb_id": "1ABC",
                    "input_path": str(tmp_path / "boom.cif"),
                    "structure_output": str(tmp_path / "boom.annotated.cif"),
                    "annotation_output": str(tmp_path / "boom.annotation.json"),
                    "num_volumes": 1,
                    "largest_type": "pore",
                    "largest_volume": 12.0,
                }
            ),
            stderr="",
        ),
    ]
    calls = []

    def _run(*args, **kwargs):
        calls.append((args, kwargs))
        response = responses.pop(0)
        if isinstance(response, BaseException):
            raise response
        return response

    monkeypatch.setattr(cli.subprocess, "run", _run)

    result = cli._run_isolated_analysis_worker(
        source_label="boom",
        input_path=tmp_path / "boom.cif",
        output_dir=tmp_path,
        min_voxels=2,
        min_volume=None,
        overwrite=False,
        assembly_policy="biological",
        resolution=3.0,
        keep_non_protein=False,
        backend="native",
        surface_connectivity="18",
        merge_mouth_gap_voxels=1,
        max_residues=10000,
        worker_timeout_seconds=12.5,
    )

    assert result["source"] == "boom"
    assert len(calls) == 2
    assert calls[0][1]["timeout"] == 12.5
    assert calls[1][1]["timeout"] == 12.5
    assert "--overwrite" not in calls[0][0][0]
    assert "--overwrite" in calls[1][0][0]


def test_run_isolated_analysis_worker_falls_back_to_python_backend(
    monkeypatch,
    tmp_path: Path,
):
    responses = [
        SimpleNamespace(returncode=-11, stdout="", stderr=""),
        SimpleNamespace(returncode=-11, stdout="", stderr=""),
        SimpleNamespace(
            returncode=0,
            stdout=json.dumps(
                {
                    "source": "boom",
                    "pdb_id": "1ABC",
                    "input_path": str(tmp_path / "boom.cif"),
                    "structure_output": str(tmp_path / "boom.annotated.cif"),
                    "annotation_output": str(tmp_path / "boom.annotation.json"),
                    "num_volumes": 1,
                    "largest_type": "pore",
                    "largest_volume": 12.0,
                }
            ),
            stderr="",
        ),
    ]
    calls = []

    def _run(*args, **kwargs):
        calls.append((args, kwargs))
        return responses.pop(0)

    monkeypatch.setattr(cli.subprocess, "run", _run)

    result = cli._run_isolated_analysis_worker(
        source_label="boom",
        input_path=tmp_path / "boom.cif",
        output_dir=tmp_path,
        min_voxels=2,
        min_volume=None,
        overwrite=False,
        assembly_policy="biological",
        resolution=3.0,
        keep_non_protein=False,
        backend="native",
        surface_connectivity="custom18",
        merge_mouth_gap_voxels=1,
        max_residues=10000,
    )

    assert result["source"] == "boom"
    assert len(calls) == 3
    backends = [
        call[0][0][call[0][0].index("--backend") + 1]
        for call in calls
    ]
    assert backends == ["native", "native", "python"]
    assert "--overwrite" not in calls[0][0][0]
    assert "--overwrite" in calls[1][0][0]
    assert "--overwrite" in calls[2][0][0]


def test_run_isolated_analysis_worker_reports_post_assembly_residue_limit(
    monkeypatch,
    tmp_path: Path,
):
    monkeypatch.setattr(
        cli.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(
            returncode=0,
            stdout=json.dumps(
                {
                    "status": "skipped_post_assembly_residue_limit",
                    "actual_residues": 12000,
                    "max_residues": 10000,
                    "assembly_policy": "biological",
                }
            ),
            stderr="",
        ),
    )

    try:
        cli._run_isolated_analysis_worker(
            source_label="boom",
            input_path=tmp_path / "boom.cif",
            output_dir=tmp_path,
            min_voxels=2,
            min_volume=None,
            overwrite=False,
            assembly_policy="biological",
            resolution=3.0,
            keep_non_protein=False,
            backend="native",
            surface_connectivity="custom18",
            merge_mouth_gap_voxels=1,
            max_residues=10000,
        )
        assert False, "expected PostAssemblyResidueLimitExceeded"
    except cli.PostAssemblyResidueLimitExceeded as error:
        assert error.actual_residues == 12000
        assert error.max_residues == 10000


def test_run_analysis_command_rejects_negative_worker_timeout(tmp_path: Path):
    try:
        cli._run_analysis_command(
            _make_args(
                tmp_path,
                worker_timeout_seconds=-1.0,
            )
        )
        assert False, "expected ValueError"
    except ValueError as error:
        assert "--worker-timeout-seconds must be >= 0." == str(error)


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

    def _analyze(source_label, input_path, output_dir, min_voxels, min_volume, overwrite, assembly_policy="biological"):
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

    def _analyze(source_label, input_path, output_dir, min_voxels, min_volume, overwrite, assembly_policy="biological"):
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


def test_run_cli_resume_ignores_checkpoint_state_without_outputs(monkeypatch, tmp_path: Path):
    checkpoint_path = tmp_path / "state.checkpoint.json"

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [
            ("first", tmp_path / "first.cif"),
            ("second", tmp_path / "second.cif"),
        ],
    )

    def _analyze_first_pass(source_label, input_path, output_dir, min_voxels, min_volume, overwrite, assembly_policy="biological"):
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
    (tmp_path / "first.annotated.cif").unlink()
    (tmp_path / "first.annotation.json").unlink()

    analyzed_sources = []

    def _analyze_second_pass(source_label, input_path, output_dir, min_voxels, min_volume, overwrite, assembly_policy="biological"):
        analyzed_sources.append(source_label)
        structure_out = output_dir / f"{source_label}.annotated.cif"
        annotation_out = output_dir / f"{source_label}.annotation.json"
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
    assert analyzed_sources == ["first", "second"]


def test_analyze_structure_file_removes_partial_outputs_before_reanalysis(
    monkeypatch,
    tmp_path: Path,
):
    input_path = tmp_path / "sample.cif"
    input_path.write_text("dummy", encoding="utf-8")

    stale_structure_output = tmp_path / "sample.annotated.cif"
    stale_structure_output.write_text("stale", encoding="utf-8")
    annotation_output = tmp_path / "sample.annotation.json"

    dummy_pdb = SimpleNamespace(
        load_structure=lambda input_path, assembly_policy="biological": "input",
        save_structure=lambda structure, output_path: Path(output_path).write_text(
            "fresh-structure",
            encoding="utf-8",
        ),
        ensure_b_factor_annotation=lambda structure: None,
        compute_sse_fractions=lambda prepared_structure: {
            "frac_alpha": 0.1,
            "frac_beta": 0.2,
            "frac_coil": 0.7,
        },
    )
    dummy_volumizer = SimpleNamespace(
        prepare_pdb_structure=lambda structure: "prepared",
        annotate_structure_volumes=lambda prepared_structure, min_voxels=2, min_volume=None: (
            pd.DataFrame([{"id": 1, "type": "pore", "volume": 12.0}]),
            "annotation",
        ),
    )
    dummy_utils = SimpleNamespace(
        get_active_backend=lambda: "native",
        VOXEL_SIZE=3.0,
    )

    monkeypatch.setattr(cli, "_pdb_module", lambda: dummy_pdb)
    monkeypatch.setattr(cli, "_volumizer_module", lambda: dummy_volumizer)
    monkeypatch.setattr(cli, "_utils_module", lambda: dummy_utils)

    result = cli.analyze_structure_file(
        source_label="sample",
        input_path=input_path,
        output_dir=tmp_path,
        min_voxels=2,
        min_volume=None,
        overwrite=False,
        assembly_policy="biological",
    )

    assert result["num_volumes"] == 1
    assert stale_structure_output.read_text(encoding="utf-8") == "fresh-structure"
    assert annotation_output.is_file()


def test_analyze_structure_file_removes_stale_status_record_before_analysis(
    monkeypatch,
    tmp_path: Path,
):
    input_path = tmp_path / "sample.cif"
    input_path.write_text("dummy", encoding="utf-8")
    status_path = tmp_path / "sample.status.json"
    status_path.write_text(
        json.dumps(
            {
                "status_record_format": 1,
                "kind": "terminal_skip",
                "source": "sample",
                "reason": "post_assembly_residue_limit",
            }
        ),
        encoding="utf-8",
    )

    dummy_pdb = SimpleNamespace(
        load_structure=lambda input_path, assembly_policy="biological": "input",
        save_structure=lambda structure, output_path: Path(output_path).write_text(
            "fresh-structure",
            encoding="utf-8",
        ),
        ensure_b_factor_annotation=lambda structure: None,
        compute_sse_fractions=lambda prepared_structure: {
            "frac_alpha": 0.1,
            "frac_beta": 0.2,
            "frac_coil": 0.7,
        },
    )
    dummy_volumizer = SimpleNamespace(
        prepare_pdb_structure=lambda structure: "prepared",
        annotate_structure_volumes=lambda prepared_structure, min_voxels=2, min_volume=None: (
            pd.DataFrame([{"id": 1, "type": "pore", "volume": 12.0}]),
            "annotation",
        ),
    )
    dummy_utils = SimpleNamespace(
        get_active_backend=lambda: "native",
        VOXEL_SIZE=3.0,
    )

    monkeypatch.setattr(cli, "_pdb_module", lambda: dummy_pdb)
    monkeypatch.setattr(cli, "_volumizer_module", lambda: dummy_volumizer)
    monkeypatch.setattr(cli, "_utils_module", lambda: dummy_utils)

    result = cli.analyze_structure_file(
        source_label="sample",
        input_path=input_path,
        output_dir=tmp_path,
        min_voxels=2,
        min_volume=None,
        overwrite=False,
        assembly_policy="biological",
    )

    assert result["num_volumes"] == 1
    assert not status_path.exists()


def test_analyze_structure_file_enforces_post_assembly_residue_limit(
    monkeypatch,
    tmp_path: Path,
):
    input_path = tmp_path / "sample.cif"
    input_path.write_text("dummy", encoding="utf-8")

    dummy_pdb = SimpleNamespace(
        load_structure=lambda input_path, assembly_policy="biological": "input",
        get_structure_residue_count=lambda structure: 10001,
    )
    dummy_volumizer = SimpleNamespace(
        prepare_pdb_structure=lambda structure: "prepared",
        annotate_structure_volumes=lambda prepared_structure, min_voxels=2, min_volume=None: (_ for _ in ()).throw(
            AssertionError("annotation should not run for oversize assemblies")
        ),
    )

    monkeypatch.setattr(cli, "_pdb_module", lambda: dummy_pdb)
    monkeypatch.setattr(cli, "_volumizer_module", lambda: dummy_volumizer)

    try:
        cli.analyze_structure_file(
            source_label="sample",
            input_path=input_path,
            output_dir=tmp_path,
            min_voxels=2,
            min_volume=None,
            overwrite=False,
            assembly_policy="biological",
            max_residues=10000,
        )
        assert False, "expected PostAssemblyResidueLimitExceeded"
    except cli.PostAssemblyResidueLimitExceeded as error:
        assert error.actual_residues == 10001
        assert error.max_residues == 10000


def test_run_cli_reports_exception_type_when_message_is_empty(
    monkeypatch,
    tmp_path: Path,
    capsys,
):
    args = _make_args(tmp_path, command="analyze", jobs=1)

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [
            ("blank", tmp_path / "blank.cif"),
        ],
    )

    def _analyze(
        source_label,
        input_path,
        output_dir,
        min_voxels,
        min_volume,
        overwrite,
        assembly_policy="biological",
    ):
        raise StopIteration()

    monkeypatch.setattr(cli, "analyze_structure_file", _analyze)

    exit_code = cli.run_cli(args)
    assert exit_code == 1

    stderr = capsys.readouterr().err
    assert "error for blank: StopIteration()" in stderr

    summary = json.loads((tmp_path / "run.summary.json").read_text(encoding="utf-8"))
    assert summary["num_failed"] == 1
    assert summary["errors"][0]["error"] == "StopIteration()"
    assert summary["errors"][0]["error_type"] == "StopIteration"
    assert "StopIteration" in summary["errors"][0]["traceback"]


def test_cache_subcommand_inspect_and_clear_negative(tmp_path: Path):
    cache_path = tmp_path / "entry_metadata"
    cli._write_metadata_record(
        cache_path,
        "1ABC",
        entry_metadata={"dummy": 1},
    )
    cli._write_metadata_record(
        cache_path,
        "2DEF",
        negative_entry={"status_code": 404, "reason": "permanent_metadata_error"},
    )
    cli._write_metadata_record(
        cache_path,
        "3GHI",
        negative_entry={"status_code": 410, "reason": "permanent_metadata_error"},
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

    entries, negative_entries = cli._load_metadata_cache(cache_path)
    assert "1ABC" in entries
    assert negative_entries == {}


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
                        "cluster_member_pdb_ids": ["1ABC", "1ABD"],
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
    assert summary["results"][0]["cluster_member_pdb_ids"] == ["1ABC", "1ABD"]


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
    monkeypatch.delenv("VOLUMIZER_BACKEND", raising=False)
    _patch_cluster_fetch(
        monkeypatch,
        {
            "1ABC": ["1ABC", "1ABD"],
        },
    )
    monkeypatch.setattr(
        rcsb,
        "fetch_entry_metadata",
        lambda pdb_id, timeout=60.0, retries=0, retry_delay=1.0: {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.0]},
        },
    )
    monkeypatch.setattr(
        cli,
        "_native_backend_module",
        lambda: (_ for _ in ()).throw(
            AssertionError("dry-run cluster should not initialize backend")
        ),
    )
    monkeypatch.setattr(
        cli,
        "_utils_module",
        lambda: (_ for _ in ()).throw(
            AssertionError("dry-run cluster should not import analysis utils")
        ),
    )
    monkeypatch.setattr(
        cli,
        "_pdb_module",
        lambda: (_ for _ in ()).throw(
            AssertionError("dry-run cluster should not import pdb helpers")
        ),
    )
    monkeypatch.setattr(
        cli,
        "_volumizer_module",
        lambda: (_ for _ in ()).throw(
            AssertionError("dry-run cluster should not import volumizer pipeline")
        ),
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
    assert summary["config"]["backend"] == cli.BACKEND_AUTO
    assert summary["num_processed"] == 0
    assert summary["num_planned"] == 1
    assert summary["planned"][0]["cluster_member_pdb_ids"] == ["1ABC", "1ABD"]


def test_run_analysis_command_passes_cluster_metadata_store_dir(monkeypatch, tmp_path: Path):
    cache_path = tmp_path / "entry_metadata"

    seen: dict[str, Path | None] = {}

    def _resolve_inputs(args, download_dir, output_dir, metadata_cache_path=None):
        seen["metadata_cache_path"] = metadata_cache_path
        return []

    monkeypatch.setattr(cli, "resolve_input_structures", _resolve_inputs)

    exit_code = cli._run_analysis_command(
        _make_args(
            tmp_path,
            command="cluster",
            cluster_identity=30,
            metadata_cache=cache_path,
            dry_run=True,
        )
    )

    assert exit_code == 0
    assert seen["metadata_cache_path"] == cache_path

    summary = json.loads((tmp_path / "run.summary.json").read_text(encoding="utf-8"))
    assert summary["config"]["metadata_cache"] == str(cache_path)


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
                "input_path": str(TEST_PDB),
                "output_structure_cif": str(structure_path),
                "num_volumes": 3,
                "volumes": [{}, {}, {}],
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


def test_run_cli_resume_emits_compact_resume_summary(
    monkeypatch,
    tmp_path: Path,
    capsys,
):
    args = _make_args(tmp_path, command="analyze", resume=True, jobs=1)

    first_structure = tmp_path / "first.annotated.cif"
    first_annotation = tmp_path / "first.annotation.json"
    first_structure.write_text("data_dummy", encoding="utf-8")
    first_annotation.write_text(
        json.dumps(
            {
                "source": "first",
                "input_path": str(tmp_path / "first.cif"),
                "output_structure_cif": str(first_structure),
                "num_volumes": 1,
                "volumes": [{}],
            }
        ),
        encoding="utf-8",
    )

    second_structure = tmp_path / "second.annotated.cif"
    second_annotation = tmp_path / "second.annotation.json"
    second_structure.write_text("data_dummy", encoding="utf-8")
    second_annotation.write_text(
        json.dumps(
            {
                "source": "second",
                "input_path": str(tmp_path / "second.cif"),
                "output_structure_cif": str(second_structure),
                "num_volumes": 0,
                "volumes": [],
            }
        ),
        encoding="utf-8",
    )

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
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("valid resumed outputs should be skipped")
        ),
    )

    exit_code = cli.run_cli(args)
    assert exit_code == 0

    stderr = capsys.readouterr().err
    assert (
        "resume scan: skipped_valid=2, skipped_terminal=0, "
        "reanalyze_partial=0, reanalyze_invalid=0"
    ) in stderr
    assert "skipping first: outputs already exist (--resume)" not in stderr
    assert "skipping second: outputs already exist (--resume)" not in stderr


def test_run_cli_resume_emits_progress_during_resume_scan(
    monkeypatch,
    tmp_path: Path,
    capsys,
):
    args = _make_args(
        tmp_path,
        command="analyze",
        resume=True,
        jobs=1,
        progress_interval=0.01,
    )

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [
            ("first", tmp_path / "first.cif"),
            ("second", tmp_path / "second.cif"),
        ],
    )

    def _resume_action(*, source_label, input_path, output_dir, assembly_policy, max_residues):
        time.sleep(0.05)
        return {
            "action": "skip",
            "reason_kind": "valid",
            "skip_entry": {
                "source": source_label,
                "input_path": str(input_path),
                "reason": "resume_existing_outputs",
            },
            "overwrite_existing_outputs": False,
            "log_message": None,
        }

    monkeypatch.setattr(cli, "_get_resume_action", _resume_action)
    monkeypatch.setattr(
        cli,
        "analyze_structure_file",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("resume scan should skip all structures")
        ),
    )

    exit_code = cli.run_cli(args)
    assert exit_code == 0

    stderr = capsys.readouterr().err
    assert "resume scan: validating existing outputs for 2 structures..." in stderr
    assert "resume scan progress:" in stderr
    assert (
        "resume scan: skipped_valid=2, skipped_terminal=0, "
        "reanalyze_partial=0, reanalyze_invalid=0"
    ) in stderr


def test_run_cli_resume_reanalyzes_invalid_existing_outputs(
    monkeypatch,
    tmp_path: Path,
    capsys,
):
    args = _make_args(tmp_path, command="analyze", resume=True, jobs=1)

    structure_path = tmp_path / "cavity.annotated.cif"
    annotation_path = tmp_path / "cavity.annotation.json"
    structure_path.write_text("data_dummy", encoding="utf-8")
    annotation_path.write_text(
        json.dumps(
            {
                "source": "cavity",
                "input_path": str(TEST_PDB),
                "output_structure_cif": str(structure_path),
                "num_volumes": 3,
                "volumes": [{}],
                "largest_type": "buried",
                "largest_volume": 123.0,
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir: [("cavity", TEST_PDB)],
    )

    seen_overwrite_flags: list[bool] = []

    def _analyze(
        source_label,
        input_path,
        output_dir,
        min_voxels,
        min_volume,
        overwrite,
        assembly_policy="biological",
    ):
        seen_overwrite_flags.append(bool(overwrite))
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
    assert seen_overwrite_flags == [True]

    stderr = capsys.readouterr().err
    assert (
        "resume scan: skipped_valid=0, skipped_terminal=0, "
        "reanalyze_partial=0, reanalyze_invalid=1"
    ) in stderr
    assert "resume invalid examples: cavity (" in stderr
    assert "annotation num_volumes does not match volumes length" in stderr

    summary = json.loads((tmp_path / "run.summary.json").read_text(encoding="utf-8"))
    assert summary["num_processed"] == 1
    assert summary["num_failed"] == 0
    assert summary["num_skipped"] == 0


def test_run_cli_resume_skips_cached_post_assembly_limit_status(
    monkeypatch,
    tmp_path: Path,
    capsys,
):
    args = _make_args(
        tmp_path,
        command="cluster",
        cluster_identity=30,
        resume=True,
        jobs=1,
        cluster_max_residues=10000,
    )

    status_path = tmp_path / "oversize.status.json"
    status_path.write_text(
        json.dumps(
            {
                "status_record_format": 1,
                "kind": "terminal_skip",
                "source": "oversize",
                "pdb_id": "1ABC",
                "input_path": str(tmp_path / "oversize.cif"),
                "structure_output": str(tmp_path / "oversize.annotated.cif"),
                "annotation_output": str(tmp_path / "oversize.annotation.json"),
                "reason": "post_assembly_residue_limit",
                "actual_residues": 12000,
                "max_residues": 10000,
                "assembly_policy": "biological",
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir, metadata_cache_path=None: [
            ("oversize", tmp_path / "oversize.cif"),
        ],
    )
    monkeypatch.setattr(
        cli,
        "_native_backend_module",
        lambda: SimpleNamespace(
            BACKEND_ENV=cli.BACKEND_ENV,
            clear_backend_cache=lambda: None,
            active_backend=lambda: "native",
        ),
    )
    monkeypatch.setattr(
        cli,
        "_utils_module",
        lambda: SimpleNamespace(
            set_resolution=lambda resolution: None,
            set_non_protein=lambda keep_non_protein: None,
            set_surface_component_connectivity_mode=lambda mode: None,
            set_surface_mouth_merge_gap_voxels=lambda gap: None,
        ),
    )
    monkeypatch.setattr(
        cli,
        "_run_isolated_analysis_worker",
        lambda **kwargs: (_ for _ in ()).throw(
            AssertionError("cached terminal skip should avoid analysis")
        ),
    )

    exit_code = cli.run_cli(args)
    assert exit_code == 0

    summary = json.loads((tmp_path / "run.summary.json").read_text(encoding="utf-8"))
    assert summary["num_processed"] == 0
    assert summary["num_failed"] == 0
    assert summary["num_skipped"] == 1
    assert summary["skipped"][0]["reason"] == "post_assembly_residue_limit"
    assert summary["skipped"][0]["actual_residues"] == 12000
    assert summary["skipped"][0]["max_residues"] == 10000
    assert summary["skipped"][0]["status_kind"] == "terminal_skip"

    stderr = capsys.readouterr().err
    assert (
        "resume scan: skipped_valid=0, skipped_terminal=1, "
        "reanalyze_partial=0, reanalyze_invalid=0"
    ) in stderr
    assert "resume terminal examples: oversize (" in stderr


def test_run_cli_resume_ignores_cached_post_assembly_limit_when_threshold_relaxed(
    monkeypatch,
    tmp_path: Path,
):
    args = _make_args(
        tmp_path,
        command="cluster",
        cluster_identity=30,
        resume=True,
        jobs=1,
        cluster_max_residues=20000,
    )

    status_path = tmp_path / "oversize.status.json"
    status_path.write_text(
        json.dumps(
            {
                "status_record_format": 1,
                "kind": "terminal_skip",
                "source": "oversize",
                "pdb_id": "1ABC",
                "input_path": str(tmp_path / "oversize.cif"),
                "structure_output": str(tmp_path / "oversize.annotated.cif"),
                "annotation_output": str(tmp_path / "oversize.annotation.json"),
                "reason": "post_assembly_residue_limit",
                "actual_residues": 12000,
                "max_residues": 10000,
                "assembly_policy": "biological",
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        cli,
        "resolve_input_structures",
        lambda args, download_dir, output_dir, metadata_cache_path=None: [
            ("oversize", tmp_path / "oversize.cif"),
        ],
    )
    monkeypatch.setattr(
        cli,
        "_native_backend_module",
        lambda: SimpleNamespace(
            BACKEND_ENV=cli.BACKEND_ENV,
            clear_backend_cache=lambda: None,
            active_backend=lambda: "python",
        ),
    )
    monkeypatch.setattr(
        cli,
        "_utils_module",
        lambda: SimpleNamespace(
            set_resolution=lambda resolution: None,
            set_non_protein=lambda keep_non_protein: None,
            set_surface_component_connectivity_mode=lambda mode: None,
            set_surface_mouth_merge_gap_voxels=lambda gap: None,
        ),
    )

    analyze_calls: list[str] = []

    def _analyze(
        source_label,
        input_path,
        output_dir,
        min_voxels,
        min_volume,
        overwrite,
        assembly_policy="biological",
        max_residues=None,
    ):
        analyze_calls.append(source_label)
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
    assert analyze_calls == ["oversize"]
