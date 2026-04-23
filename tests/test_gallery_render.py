import json
from pathlib import Path
import sqlite3
import threading
import time

import pytest

from volumizer import gallery_index, gallery_render
from volumizer.paths import TEST_DIR


TEST_INPUT_PDB = TEST_DIR / "pdbs" / "cavity.pdb"


def _write_annotation(path: Path, volume: float) -> None:
    path.write_text(
        json.dumps(
            {
                "source": path.stem,
                "num_volumes": 1,
                "volumes": [
                    {
                        "id": 0,
                        "type": "pore",
                        "volume": volume,
                        "x": 12.0,
                        "y": 5.0,
                        "z": 3.0,
                    }
                ],
            }
        ),
        encoding="utf-8",
    )


def _build_render_fixture_db(
    tmp_path: Path,
    *,
    run_id: str = "render-fixture",
    missing_structures: set[str] | None = None,
    source_labels: list[str] | None = None,
) -> tuple[Path, str, Path]:
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)

    missing = missing_structures or set()
    labels = source_labels or ["hit-a", "hit-b"]

    results = []
    for index, source_label in enumerate(labels):
        volume = 200.0 - (index * 10.0)
        annotation_path = run_dir / f"{source_label}.annotation.json"
        structure_output_path = run_dir / f"{source_label}.annotated.cif"
        _write_annotation(annotation_path, volume)

        if source_label not in missing:
            structure_output_path.write_text("data_test\n#\n", encoding="utf-8")

        results.append(
            {
                "source": source_label,
                "input_path": str(TEST_INPUT_PDB),
                "structure_output": str(structure_output_path),
                "annotation_output": str(annotation_path),
            }
        )

    summary_path = run_dir / "run.summary.json"
    summary_path.write_text(
        json.dumps(
            {
                "config": {
                    "assembly_policy": "biological",
                    "resolution": 3.0,
                    "keep_non_protein": False,
                    "output_dir": str(run_dir),
                },
                "results": results,
                "errors": [],
                "skipped": [],
                "planned": [],
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    db_path = tmp_path / "gallery.db"
    gallery_index.build_gallery_index(
        summary_path=summary_path,
        db_path=db_path,
        run_id=run_id,
        replace_run=False,
        strict=True,
    )

    return db_path, run_id, run_dir


def test_render_gallery_thumbnails_updates_rows_and_paths(tmp_path: Path):
    db_path, run_id, _ = _build_render_fixture_db(tmp_path)
    render_root = tmp_path / "renders"

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        assert structure_path.is_file()
        assert width == 160
        assert height == 120
        assert style["width"] == 160
        assert style["height"] == 120
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(f"png-{axis}".encode("utf-8"))

    result = gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        width=160,
        height=120,
        render_fn=fake_render,
    )

    assert result["total_jobs"] == 2
    assert result["rendered_jobs"] == 2
    assert result["failed_jobs"] == 0
    assert result["skipped_jobs"] == 0

    with sqlite3.connect(db_path) as connection:
        rows = connection.execute(
            """
            SELECT s.source_label, r.render_status, r.x_png_path, r.y_png_path, r.z_png_path,
                   r.render_style_hash, r.render_error
            FROM renders r
            INNER JOIN structures s ON s.structure_id = r.structure_id
            ORDER BY s.source_label ASC
            """
        ).fetchall()

    assert len(rows) == 2
    for source_label, status, x_path, y_path, z_path, style_hash, render_error in rows:
        assert source_label in {"hit-a", "hit-b"}
        assert status == "done"
        assert render_error is None
        assert style_hash == result["style_hash"]
        for raw_path in (x_path, y_path, z_path):
            assert Path(raw_path).is_file()


def test_render_gallery_thumbnails_skips_completed_rows_when_cache_is_fresh(tmp_path: Path):
    db_path, run_id, _ = _build_render_fixture_db(tmp_path)
    render_root = tmp_path / "renders"

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(b"ok")

    gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        render_fn=fake_render,
    )

    second = gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        render_fn=fake_render,
    )

    assert second["total_jobs"] == 0
    assert second["rendered_jobs"] == 0
    assert second["failed_jobs"] == 0
    assert second["skipped_jobs"] == 2


def test_render_gallery_thumbnails_rerenders_missing_cached_file(tmp_path: Path):
    db_path, run_id, _ = _build_render_fixture_db(tmp_path)
    render_root = tmp_path / "renders"

    calls: list[str] = []

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        calls.append(output_dir.name)
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(b"ok")

    gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        render_fn=fake_render,
    )

    missing_path = render_root / run_id / "hit-a" / "y.png"
    missing_path.unlink()

    calls.clear()
    second = gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        render_fn=fake_render,
    )

    assert second["total_jobs"] == 1
    assert second["rendered_jobs"] == 1
    assert second["failed_jobs"] == 0
    assert second["skipped_jobs"] == 1
    assert calls == ["hit-a"]


def test_render_gallery_thumbnails_rerenders_when_style_hash_changes(tmp_path: Path):
    db_path, run_id, _ = _build_render_fixture_db(tmp_path)
    render_root = tmp_path / "renders"

    calls: list[str] = []

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        calls.append(output_dir.name)
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(b"ok")

    first = gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        render_fn=fake_render,
    )

    calls.clear()
    second = gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        style={"background_hex": "#f6f6f6"},
        render_fn=fake_render,
    )

    assert second["style_hash"] != first["style_hash"]
    assert second["total_jobs"] == 2
    assert second["rendered_jobs"] == 2
    assert second["failed_jobs"] == 0
    assert second["skipped_jobs"] == 0
    assert calls == ["hit-a", "hit-b"]


def test_render_gallery_thumbnails_reuses_existing_outputs_without_rendering(tmp_path: Path):
    db_path, run_id, _ = _build_render_fixture_db(tmp_path)
    render_root = tmp_path / "renders"

    for source_label in ("hit-a", "hit-b"):
        output_dir = render_root / run_id / source_label
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(f"cached-{axis}".encode("utf-8"))

    result = gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        reuse_existing_only=True,
        render_fn=lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("render_fn should not be called in reuse_existing_only mode")
        ),
    )

    assert result["total_jobs"] == 0
    assert result["rendered_jobs"] == 0
    assert result["failed_jobs"] == 0
    assert result["reused_jobs"] == 2
    assert result["skipped_jobs"] == 0

    with sqlite3.connect(db_path) as connection:
        rows = connection.execute(
            """
            SELECT s.source_label, r.render_status, r.x_png_path, r.y_png_path, r.z_png_path,
                   r.render_style_hash, r.render_error
            FROM renders r
            INNER JOIN structures s ON s.structure_id = r.structure_id
            ORDER BY s.source_label ASC
            """
        ).fetchall()

    assert len(rows) == 2
    for source_label, status, x_path, y_path, z_path, style_hash, render_error in rows:
        assert source_label in {"hit-a", "hit-b"}
        assert status == "done"
        assert render_error is None
        assert style_hash == result["style_hash"]
        for raw_path in (x_path, y_path, z_path):
            assert Path(raw_path).is_file()


def test_render_gallery_thumbnails_handles_failed_and_missing_structures(tmp_path: Path):
    db_path, run_id, _ = _build_render_fixture_db(
        tmp_path,
        missing_structures={"hit-b"},
    )
    render_root = tmp_path / "renders"

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(b"ok")

    first = gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        render_fn=fake_render,
    )

    assert first["total_jobs"] == 2
    assert first["rendered_jobs"] == 1
    assert first["failed_jobs"] == 1
    assert first["skipped_jobs"] == 0

    with sqlite3.connect(db_path) as connection:
        failed_row = connection.execute(
            """
            SELECT s.source_label, r.render_status, r.render_error
            FROM renders r
            INNER JOIN structures s ON s.structure_id = r.structure_id
            WHERE s.source_label = 'hit-b'
            """
        ).fetchone()

        assert failed_row is not None
        assert failed_row[0] == "hit-b"
        assert failed_row[1] == "failed"
        assert "annotated structure not found" in failed_row[2]

        connection.execute(
            """
            UPDATE renders
            SET render_status = 'failed', render_error = 'synthetic failure'
            WHERE structure_id IN (
                SELECT structure_id FROM structures WHERE source_label = 'hit-a'
            )
            """
        )
        connection.commit()

    skipped_failed = gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        include_failed=False,
        render_fn=fake_render,
    )
    assert skipped_failed["total_jobs"] == 0
    assert skipped_failed["rendered_jobs"] == 0
    assert skipped_failed["failed_jobs"] == 0
    assert skipped_failed["skipped_jobs"] == 2

    include_failed = gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        include_failed=True,
        render_fn=fake_render,
    )
    assert include_failed["total_jobs"] == 2
    assert include_failed["rendered_jobs"] == 1
    assert include_failed["failed_jobs"] == 1
    assert include_failed["skipped_jobs"] == 0


def test_render_gallery_thumbnails_emits_progress_with_rate_and_eta(
    tmp_path: Path,
    capsys,
):
    db_path, run_id, _ = _build_render_fixture_db(tmp_path)
    render_root = tmp_path / "renders"

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(b"ok")

    gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        render_fn=fake_render,
    )

    stderr = capsys.readouterr().err
    assert "[gallery-render] progress:" in stderr
    assert "rate=" in stderr
    assert "eta=" in stderr
    assert "skipped=" in stderr


def test_render_gallery_thumbnails_progress_interval_zero_disables_progress(
    tmp_path: Path,
    capsys,
):
    db_path, run_id, _ = _build_render_fixture_db(tmp_path)
    render_root = tmp_path / "renders"

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(b"ok")

    gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        render_fn=fake_render,
        progress_interval=0.0,
    )

    stderr = capsys.readouterr().err
    assert "[gallery-render] progress:" not in stderr


def test_render_gallery_thumbnails_emits_heartbeat_before_first_completion(
    tmp_path: Path,
    capsys,
):
    db_path, run_id, _ = _build_render_fixture_db(tmp_path)
    render_root = tmp_path / "renders"

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        time.sleep(0.05)
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(b"ok")

    gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        limit=1,
        render_fn=fake_render,
        progress_interval=0.01,
    )

    stderr = capsys.readouterr().err
    assert stderr.count("[gallery-render] progress: 0/1") >= 2


def test_render_gallery_thumbnails_parallel_mixed_results_keep_db_updates_safe(
    tmp_path: Path,
):
    db_path, run_id, _ = _build_render_fixture_db(tmp_path)
    render_root = tmp_path / "renders"
    started: list[str] = []
    barrier = threading.Barrier(2, timeout=2.0)

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        started.append(output_dir.name)
        barrier.wait()
        if output_dir.name == "hit-a":
            time.sleep(0.05)
            output_dir.mkdir(parents=True, exist_ok=True)
            for axis in ("x", "y", "z"):
                (output_dir / f"{axis}.png").write_bytes(b"ok")
            return
        raise RuntimeError("synthetic parallel failure")

    result = gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        jobs=2,
        render_fn=fake_render,
    )

    assert sorted(started) == ["hit-a", "hit-b"]
    assert result["jobs"] == 2
    assert result["total_jobs"] == 2
    assert result["rendered_jobs"] == 1
    assert result["failed_jobs"] == 1

    with sqlite3.connect(db_path) as connection:
        rows = connection.execute(
            """
            SELECT s.source_label, r.render_status, r.render_error
            FROM renders r
            INNER JOIN structures s ON s.structure_id = r.structure_id
            ORDER BY s.source_label ASC
            """
        ).fetchall()

    assert rows[0][0] == "hit-a"
    assert rows[0][1] == "done"
    assert rows[0][2] is None
    assert rows[1][0] == "hit-b"
    assert rows[1][1] == "failed"
    assert rows[1][2] == "synthetic parallel failure"


def test_render_gallery_thumbnails_parallel_only_creates_inflight_output_dirs(
    tmp_path: Path,
):
    source_labels = [f"hit-{index}" for index in range(6)]
    db_path, run_id, _ = _build_render_fixture_db(
        tmp_path,
        source_labels=source_labels,
    )
    render_root = tmp_path / "renders"
    render_run_dir = render_root / run_id
    entered_count = 0
    entered_event = threading.Event()
    release_event = threading.Event()
    result_holder: dict[str, dict] = {}
    error_holder: dict[str, BaseException] = {}
    lock = threading.Lock()

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        nonlocal entered_count
        with lock:
            entered_count += 1
            if entered_count >= 2:
                entered_event.set()
        assert release_event.wait(timeout=2.0)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(b"ok")

    def run_render() -> None:
        try:
            result_holder["result"] = gallery_render.render_gallery_thumbnails(
                db_path=db_path,
                render_root=render_root,
                run_id=run_id,
                jobs=2,
                render_fn=fake_render,
                progress_interval=0.0,
            )
        except BaseException as error:  # pragma: no cover - test coordination
            error_holder["error"] = error

    worker = threading.Thread(target=run_render)
    worker.start()

    assert entered_event.wait(timeout=2.0)
    time.sleep(0.05)

    existing_dirs = sorted(
        path.name for path in render_run_dir.iterdir() if path.is_dir()
    )
    assert len(existing_dirs) == 2

    release_event.set()
    worker.join(timeout=5.0)
    assert not worker.is_alive()
    assert "error" not in error_holder
    assert result_holder["result"]["rendered_jobs"] == len(source_labels)


def test_render_gallery_thumbnails_writes_timing_jsonl(tmp_path: Path):
    db_path, run_id, _ = _build_render_fixture_db(tmp_path)
    render_root = tmp_path / "renders"
    timing_path = tmp_path / "timings.jsonl"

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(b"ok")

    gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        timing_jsonl=timing_path,
        render_fn=fake_render,
    )

    events = [
        json.loads(line)
        for line in timing_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]

    assert events[0]["event"] == "planning"
    assert events[-1]["event"] == "summary"
    assert sum(1 for item in events if item["event"] == "job") == 2

    planning = next(item for item in events if item["event"] == "planning")
    assert planning["queue_scan_seconds"] >= 0
    assert planning["jobs"] == 1

    job_events = [item for item in events if item["event"] == "job"]
    assert all(item["renderer"] == "custom" for item in job_events)
    assert all(item["wall_seconds"] >= 0 for item in job_events)
    assert all(item["dispatch_wait_seconds"] >= 0 for item in job_events)

    summary = events[-1]
    assert summary["rendered_jobs"] == 2
    assert summary["failed_jobs"] == 0
    assert summary["dispatch_seconds"] >= 0
    assert summary["wall_seconds"] >= 0


def test_render_gallery_thumbnails_rejects_invalid_jobs_and_timeout(tmp_path: Path):
    db_path, run_id, _ = _build_render_fixture_db(tmp_path)
    render_root = tmp_path / "renders"

    with pytest.raises(ValueError, match="jobs must be >= 1"):
        gallery_render.render_gallery_thumbnails(
            db_path=db_path,
            render_root=render_root,
            run_id=run_id,
            jobs=0,
        )

    with pytest.raises(ValueError, match="worker_timeout_seconds must be >= 0\\."):
        gallery_render.render_gallery_thumbnails(
            db_path=db_path,
            render_root=render_root,
            run_id=run_id,
            worker_timeout_seconds=-1.0,
        )


def test_coerce_absolute_path_keeps_absolute_and_expands_relative():
    relative = Path("data") / "sample.cif"
    absolute = gallery_render._coerce_absolute_path(str(relative))
    assert absolute == Path.cwd() / relative

    already_absolute = Path("/tmp/sample.cif")
    assert gallery_render._coerce_absolute_path(str(already_absolute)) == already_absolute


def test_render_single_structure_with_node_passes_backend_flags_and_parses_timing(
    tmp_path: Path,
    monkeypatch,
):
    structure_path = tmp_path / "sample.cif"
    structure_path.write_text("data_test\n#\n", encoding="utf-8")
    output_dir = tmp_path / "renders"
    captured: dict[str, object] = {}

    def fake_run(command: list[str], **kwargs):
        captured["command"] = command
        captured["timeout"] = kwargs.get("timeout")

        class _Completed:
            returncode = 0
            stdout = (
                '__VOLUMIZER_TIMING__ '
                '{"requested_backend":"auto","backend_used":"software","fallback_used":true}\n'
            )
            stderr = ""

        return _Completed()

    monkeypatch.setattr(gallery_render.subprocess, "run", fake_run)

    result = gallery_render._render_single_structure_with_node(
        renderer_script=tmp_path / "renderer.mjs",
        node_executable="node",
        structure_path=structure_path,
        output_dir=output_dir,
        width=320,
        height=240,
        style={"background_hex": "#ffffff"},
        render_backend="auto",
        axis_render_mode="fast",
        worker_timeout_seconds=12.5,
    )

    command = captured["command"]
    assert isinstance(command, list)
    assert "--render-backend" in command
    assert command[command.index("--render-backend") + 1] == "auto"
    assert "--axis-render-mode" in command
    assert command[command.index("--axis-render-mode") + 1] == "fast"
    assert "--annotation" not in command
    assert captured["timeout"] == 12.5
    assert result["timing"]["backend_used"] == "software"
    assert result["timing"]["fallback_used"] is True
