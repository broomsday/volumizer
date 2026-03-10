import json
from pathlib import Path
import sqlite3

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
) -> tuple[Path, str, Path]:
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)

    missing = missing_structures or set()

    results = []
    for source_label, volume in (("hit-a", 200.0), ("hit-b", 80.0)):
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
