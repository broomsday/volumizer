from pathlib import Path
import sqlite3

import pytest

from volumizer import gallery_rerender
from tests.test_gallery_render import _build_render_fixture_db


def test_rerender_selected_gallery_thumbnails_replaces_output_dir(tmp_path: Path):
    db_path, run_id, _ = _build_render_fixture_db(
        tmp_path,
        source_labels=["9bq2", "hit-b"],
    )
    render_root = tmp_path / "renders"
    output_dir = render_root / run_id / "9bq2"
    output_dir.mkdir(parents=True, exist_ok=True)
    for axis in ("x", "y", "z"):
        (output_dir / f"{axis}.png").write_bytes(b"old")
    (output_dir / "stale.txt").write_text("stale", encoding="utf-8")

    with sqlite3.connect(db_path) as connection:
        connection.execute(
            """
            UPDATE renders
            SET render_status = 'failed', render_error = 'old failure'
            WHERE structure_id IN (
                SELECT structure_id FROM structures WHERE source_label = '9bq2'
            )
            """
        )
        connection.commit()

    calls: list[str] = []

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        calls.append(output_dir.name)
        assert structure_path.is_file()
        assert not (output_dir / "stale.txt").exists()
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(f"new-{axis}".encode("utf-8"))

    result = gallery_rerender.rerender_selected_gallery_thumbnails(
        db_path=db_path,
        selectors=["9bq2"],
        render_root=render_root,
        render_fn=fake_render,
        progress_interval=0.0,
    )

    assert result["total_jobs"] == 1
    assert result["rendered_jobs"] == 1
    assert result["failed_jobs"] == 0
    assert calls == ["9bq2"]
    assert not (output_dir / "stale.txt").exists()
    for axis in ("x", "y", "z"):
        assert (output_dir / f"{axis}.png").read_bytes() == f"new-{axis}".encode("utf-8")

    with sqlite3.connect(db_path) as connection:
        row = connection.execute(
            """
            SELECT r.render_status, r.render_error, r.x_png_path, r.y_png_path, r.z_png_path
            FROM renders r
            INNER JOIN structures s ON s.structure_id = r.structure_id
            WHERE s.source_label = '9bq2'
            """
        ).fetchone()

    assert row == (
        "done",
        None,
        str(output_dir / "x.png"),
        str(output_dir / "y.png"),
        str(output_dir / "z.png"),
    )


def test_rerender_selected_gallery_thumbnails_matches_pdb_id_and_alias(tmp_path: Path):
    db_path, run_id, _ = _build_render_fixture_db(
        tmp_path,
        source_labels=["cluster-target"],
    )
    render_root = tmp_path / "renders"

    with sqlite3.connect(db_path) as connection:
        structure_id = connection.execute(
            "SELECT structure_id FROM structures WHERE source_label = 'cluster-target'"
        ).fetchone()[0]
        connection.execute(
            "UPDATE structures SET pdb_id = '9BQ2' WHERE structure_id = ?",
            (structure_id,),
        )
        connection.execute(
            """
            INSERT INTO structure_pdb_aliases (structure_id, alias_pdb_id, alias_source)
            VALUES (?, ?, ?)
            """,
            (structure_id, "ALT9", "test"),
        )
        connection.commit()

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

    by_pdb_id = gallery_rerender.rerender_selected_gallery_thumbnails(
        db_path=db_path,
        selectors=["9bq2"],
        render_root=render_root,
        render_fn=fake_render,
        progress_interval=0.0,
    )
    assert by_pdb_id["rendered_jobs"] == 1
    assert calls == ["cluster-target"]

    calls.clear()
    by_alias = gallery_rerender.rerender_selected_gallery_thumbnails(
        db_path=db_path,
        selectors=["alt9"],
        render_root=render_root,
        render_fn=fake_render,
        progress_interval=0.0,
    )
    assert by_alias["rendered_jobs"] == 1
    assert calls == ["cluster-target"]


def test_rerender_selected_gallery_thumbnails_rejects_unmatched_selector(tmp_path: Path):
    db_path, _, _ = _build_render_fixture_db(tmp_path)

    with pytest.raises(ValueError, match="No gallery rows matched selector"):
        gallery_rerender.rerender_selected_gallery_thumbnails(
            db_path=db_path,
            selectors=["missing-id"],
            render_root=tmp_path / "renders",
            progress_interval=0.0,
        )
