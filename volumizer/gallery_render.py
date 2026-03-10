"""
Thumbnail rendering queue for gallery hits.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
import subprocess
from typing import Any, Callable


AXIS_FILENAMES = {
    "x": "x.png",
    "y": "y.png",
    "z": "z.png",
}

DEFAULT_RENDER_STYLE: dict[str, Any] = {
    "background_hex": "#ffffff",
    "width": 320,
    "height": 240,
    "camera_axes": ["x", "y", "z"],
    "representation": "default",
}


def _sanitize_path_token(token: str) -> str:
    sanitized = "".join(
        char if (char.isalnum() or char in {"-", "_"}) else "_" for char in token
    )
    sanitized = sanitized.strip("_")
    return sanitized if len(sanitized) > 0 else "item"


def _compute_style_hash(style: dict[str, Any]) -> str:
    canonical = json.dumps(style, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    return digest[:16]


def _safe_non_negative_int(value: int, name: str) -> int:
    normalized = int(value)
    if normalized < 0:
        raise ValueError(f"{name} must be >= 0")
    return normalized


def _resolve_format_from_suffix(structure_path: Path) -> str:
    suffix = structure_path.suffix.lower()
    if suffix == ".pdb":
        return "pdb"
    if suffix in {".cif", ".mmcif"}:
        return "mmcif"
    if suffix == ".bcif":
        return "bcif"
    return "mmcif"


def _should_render_row(
    *,
    status: str,
    include_failed: bool,
    force: bool,
    style_mismatch: bool,
    outputs_exist: bool,
) -> bool:
    if force:
        return True

    normalized_status = str(status).strip().lower()
    if normalized_status == "pending":
        return True
    if normalized_status == "failed":
        return include_failed
    if normalized_status == "done":
        return (not outputs_exist) or style_mismatch

    return True


def _render_single_structure_with_node(
    *,
    renderer_script: Path,
    node_executable: str,
    structure_path: Path,
    output_dir: Path,
    width: int,
    height: int,
    style: dict[str, Any],
) -> None:
    format_name = _resolve_format_from_suffix(structure_path)

    command = [
        node_executable,
        str(renderer_script),
        "--structure",
        str(structure_path),
        "--format",
        format_name,
        "--out-dir",
        str(output_dir),
        "--width",
        str(width),
        "--height",
        str(height),
        "--style-json",
        json.dumps(style, sort_keys=True),
    ]

    completed = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        message = (completed.stderr or completed.stdout or "renderer failed").strip()
        raise RuntimeError(message)


RenderFunction = Callable[[Path, Path, int, int, dict[str, Any]], None]


def render_gallery_thumbnails(
    *,
    db_path: Path,
    render_root: Path,
    run_id: str | None = None,
    limit: int | None = None,
    include_failed: bool = False,
    force: bool = False,
    width: int = 320,
    height: int = 240,
    node_executable: str = "node",
    renderer_script: Path | None = None,
    style: dict[str, Any] | None = None,
    render_fn: RenderFunction | None = None,
) -> dict[str, Any]:
    db_path = Path(db_path).resolve()
    if not db_path.is_file():
        raise FileNotFoundError(f"Gallery DB does not exist: {db_path}")

    render_root = Path(render_root).resolve()
    render_root.mkdir(parents=True, exist_ok=True)

    width = _safe_non_negative_int(width, "width")
    height = _safe_non_negative_int(height, "height")

    normalized_limit = None
    if limit is not None:
        normalized_limit = _safe_non_negative_int(limit, "limit")

    effective_style = dict(DEFAULT_RENDER_STYLE)
    if style is not None:
        effective_style.update(style)
    effective_style["width"] = width
    effective_style["height"] = height

    style_hash = _compute_style_hash(effective_style)
    now_utc = datetime.now(tz=timezone.utc).isoformat()

    if renderer_script is None:
        renderer_script = Path(__file__).resolve().parents[1] / "scripts" / "molstar_render_single.mjs"

    query = (
        "SELECT "
        "s.structure_id, s.run_id, s.source_label, s.annotated_cif_path, "
        "r.render_status, r.render_style_hash "
        "FROM structures s "
        "INNER JOIN renders r ON r.structure_id = s.structure_id "
        "WHERE s.annotated_cif_path IS NOT NULL AND length(s.annotated_cif_path) > 0"
    )
    params: list[Any] = []

    if run_id is not None:
        query += " AND s.run_id = ?"
        params.append(str(run_id))

    query += " ORDER BY s.structure_id ASC"

    total_jobs = 0
    rendered_jobs = 0
    failed_jobs = 0
    skipped_jobs = 0

    with sqlite3.connect(db_path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(query, params).fetchall()

        for row in rows:
            if normalized_limit is not None and total_jobs >= normalized_limit:
                break

            structure_id = int(row["structure_id"])
            source_label = str(row["source_label"])
            run_token = _sanitize_path_token(str(row["run_id"]))
            source_token = _sanitize_path_token(source_label)

            output_dir = render_root / run_token / source_token
            output_dir.mkdir(parents=True, exist_ok=True)
            x_path = output_dir / AXIS_FILENAMES["x"]
            y_path = output_dir / AXIS_FILENAMES["y"]
            z_path = output_dir / AXIS_FILENAMES["z"]

            status = str(row["render_status"] or "pending")
            existing_style_hash = row["render_style_hash"]
            style_mismatch = str(existing_style_hash or "") != style_hash
            outputs_exist = x_path.is_file() and y_path.is_file() and z_path.is_file()

            should_render = _should_render_row(
                status=status,
                include_failed=include_failed,
                force=force,
                style_mismatch=style_mismatch,
                outputs_exist=outputs_exist,
            )
            if not should_render:
                skipped_jobs += 1
                continue

            total_jobs += 1

            structure_path = Path(str(row["annotated_cif_path"])).resolve()
            if not structure_path.is_file():
                failed_jobs += 1
                connection.execute(
                    """
                    UPDATE renders
                    SET render_status = ?, render_error = ?, updated_at = ?
                    WHERE structure_id = ?
                    """,
                    (
                        "failed",
                        f"annotated structure not found: {structure_path}",
                        now_utc,
                        structure_id,
                    ),
                )
                continue

            try:
                if render_fn is not None:
                    render_fn(structure_path, output_dir, width, height, effective_style)
                else:
                    _render_single_structure_with_node(
                        renderer_script=Path(renderer_script),
                        node_executable=node_executable,
                        structure_path=structure_path,
                        output_dir=output_dir,
                        width=width,
                        height=height,
                        style=effective_style,
                    )

                for output_path in (x_path, y_path, z_path):
                    if not output_path.is_file():
                        raise RuntimeError(f"missing output image: {output_path}")

                connection.execute(
                    """
                    UPDATE renders
                    SET x_png_path = ?,
                        y_png_path = ?,
                        z_png_path = ?,
                        render_style_hash = ?,
                        render_status = ?,
                        render_error = NULL,
                        updated_at = ?
                    WHERE structure_id = ?
                    """,
                    (
                        str(x_path),
                        str(y_path),
                        str(z_path),
                        style_hash,
                        "done",
                        now_utc,
                        structure_id,
                    ),
                )
                rendered_jobs += 1
            except Exception as error:  # pragma: no cover - defensive path
                failed_jobs += 1
                connection.execute(
                    """
                    UPDATE renders
                    SET render_style_hash = ?,
                        render_status = ?,
                        render_error = ?,
                        updated_at = ?
                    WHERE structure_id = ?
                    """,
                    (
                        style_hash,
                        "failed",
                        str(error),
                        now_utc,
                        structure_id,
                    ),
                )

        connection.commit()

    return {
        "db": str(db_path),
        "render_root": str(render_root),
        "run_id": run_id,
        "style_hash": style_hash,
        "total_jobs": total_jobs,
        "rendered_jobs": rendered_jobs,
        "failed_jobs": failed_jobs,
        "skipped_jobs": skipped_jobs,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render gallery thumbnails (x/y/z) for indexed structures."
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=Path("data") / "gallery.db",
        help="SQLite gallery DB path (default: data/gallery.db).",
    )
    parser.add_argument(
        "--render-root",
        type=Path,
        default=Path("data") / "renders",
        help="Root directory for rendered thumbnails (default: data/renders).",
    )
    parser.add_argument("--run-id", type=str, default=None, help="Optional run-id filter.")
    parser.add_argument("--limit", type=int, default=None, help="Optional max jobs to process.")
    parser.add_argument(
        "--include-failed",
        action="store_true",
        help="Include failed rows alongside pending rows.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Ignore render status and process all matching structures.",
    )
    parser.add_argument("--width", type=int, default=320, help="Thumbnail width in pixels.")
    parser.add_argument("--height", type=int, default=240, help="Thumbnail height in pixels.")
    parser.add_argument(
        "--node-executable",
        type=str,
        default="node",
        help="Node.js executable used for Mol* renderer script.",
    )
    parser.add_argument(
        "--renderer-script",
        type=Path,
        default=Path("scripts") / "molstar_render_single.mjs",
        help="Path to Node Mol* single-structure renderer script.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    result = render_gallery_thumbnails(
        db_path=Path(args.db),
        render_root=Path(args.render_root),
        run_id=args.run_id,
        limit=args.limit,
        include_failed=bool(args.include_failed),
        force=bool(args.force),
        width=args.width,
        height=args.height,
        node_executable=args.node_executable,
        renderer_script=Path(args.renderer_script),
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
