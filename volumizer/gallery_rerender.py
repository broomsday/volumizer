"""
Targeted thumbnail rerendering for selected gallery structures.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sqlite3
from typing import Any

from volumizer import gallery_render


def _normalize_selector(value: str) -> str:
    normalized = str(value).strip().lower()
    if len(normalized) == 0:
        raise ValueError("selectors must not be empty")
    return normalized


def _resolve_selected_structures(
    *,
    db_path: Path,
    selectors: list[str],
    run_id: str | None = None,
) -> list[dict[str, Any]]:
    db_path = Path(db_path).resolve()
    if not db_path.is_file():
        raise FileNotFoundError(f"Gallery DB does not exist: {db_path}")

    normalized_to_raw: dict[str, str] = {}
    for selector in selectors:
        normalized = _normalize_selector(selector)
        normalized_to_raw.setdefault(normalized, str(selector))

    query = (
        "SELECT "
        "s.structure_id, s.run_id, s.source_label, s.pdb_id, s.annotated_cif_path, "
        "spa.alias_pdb_id "
        "FROM structures s "
        "INNER JOIN renders r ON r.structure_id = s.structure_id "
        "LEFT JOIN structure_pdb_aliases spa ON spa.structure_id = s.structure_id "
        "WHERE s.annotated_cif_path IS NOT NULL AND length(s.annotated_cif_path) > 0"
    )
    params: list[Any] = []
    if run_id is not None:
        query += " AND s.run_id = ?"
        params.append(str(run_id))
    query += " ORDER BY s.structure_id ASC"

    matched_structure_rows: dict[int, sqlite3.Row] = {}
    matched_selectors_by_structure: dict[int, set[str]] = {}
    selector_hits: dict[str, set[int]] = {key: set() for key in normalized_to_raw}

    with sqlite3.connect(db_path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(query, params).fetchall()

    normalized_selectors = set(normalized_to_raw)
    for row in rows:
        candidate_terms = {
            _normalize_selector(str(row["source_label"])),
        }
        pdb_id = row["pdb_id"]
        alias_pdb_id = row["alias_pdb_id"]
        if pdb_id:
            candidate_terms.add(_normalize_selector(str(pdb_id)))
        if alias_pdb_id:
            candidate_terms.add(_normalize_selector(str(alias_pdb_id)))

        matched_terms = candidate_terms.intersection(normalized_selectors)
        if len(matched_terms) == 0:
            continue

        structure_id = int(row["structure_id"])
        matched_structure_rows.setdefault(structure_id, row)
        matched_selectors_by_structure.setdefault(structure_id, set()).update(matched_terms)
        for matched_term in matched_terms:
            selector_hits[matched_term].add(structure_id)

    unmatched = [
        normalized_to_raw[normalized]
        for normalized, hits in selector_hits.items()
        if len(hits) == 0
    ]
    if len(unmatched) > 0:
        unmatched_text = ", ".join(sorted(unmatched))
        raise ValueError(f"No gallery rows matched selector(s): {unmatched_text}")

    resolved: list[dict[str, Any]] = []
    for structure_id in sorted(matched_structure_rows):
        row = matched_structure_rows[structure_id]
        resolved.append(
            {
                "structure_id": structure_id,
                "run_id": str(row["run_id"]),
                "source_label": str(row["source_label"]),
                "pdb_id": (str(row["pdb_id"]) if row["pdb_id"] is not None else None),
                "matched_selectors": sorted(matched_selectors_by_structure[structure_id]),
                "row": row,
            }
        )

    return resolved


def _reset_selected_render_rows(
    *,
    connection: sqlite3.Connection,
    structure_ids: list[int],
) -> None:
    placeholders = ", ".join("?" for _ in structure_ids)
    connection.execute(
        f"""
        UPDATE renders
        SET x_png_path = NULL,
            y_png_path = NULL,
            z_png_path = NULL,
            render_style_hash = NULL,
            render_status = ?,
            render_error = NULL,
            updated_at = ?
        WHERE structure_id IN ({placeholders})
        """,
        (
            "pending",
            datetime.now(tz=timezone.utc).isoformat(),
            *structure_ids,
        ),
    )


def rerender_selected_gallery_thumbnails(
    *,
    db_path: Path,
    selectors: list[str],
    render_root: Path = Path("data") / "renders",
    run_id: str | None = None,
    width: int = 320,
    height: int = 240,
    node_executable: str = "node",
    renderer_script: Path | None = None,
    style: dict[str, Any] | None = None,
    render_fn: gallery_render.RenderFunction | None = None,
    jobs: int = 1,
    render_backend: str = gallery_render.DEFAULT_RENDER_BACKEND,
    axis_render_mode: str = gallery_render.DEFAULT_AXIS_RENDER_MODE,
    worker_timeout_seconds: float | None = None,
    timing_jsonl: Path | None = None,
    progress_interval: float = gallery_render.DEFAULT_PROGRESS_INTERVAL_SECONDS,
) -> dict[str, Any]:
    db_path = Path(db_path).resolve()
    render_root = Path(render_root).resolve()
    render_root.mkdir(parents=True, exist_ok=True)

    resolved_structures = _resolve_selected_structures(
        db_path=db_path,
        selectors=selectors,
        run_id=run_id,
    )
    structure_ids = [int(item["structure_id"]) for item in resolved_structures]

    with sqlite3.connect(db_path) as connection:
        for item in resolved_structures:
            job = gallery_render._build_planned_job(
                row=item["row"],
                render_root=render_root,
            )
            if job.output_dir.exists():
                shutil.rmtree(job.output_dir)
        _reset_selected_render_rows(
            connection=connection,
            structure_ids=structure_ids,
        )
        connection.commit()

    result = gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        structure_ids=structure_ids,
        force=True,
        width=width,
        height=height,
        node_executable=node_executable,
        renderer_script=renderer_script,
        style=style,
        render_fn=render_fn,
        jobs=jobs,
        render_backend=render_backend,
        axis_render_mode=axis_render_mode,
        worker_timeout_seconds=worker_timeout_seconds,
        timing_jsonl=timing_jsonl,
        progress_interval=progress_interval,
    )
    result["selectors"] = [str(selector) for selector in selectors]
    result["matched_structures"] = [
        {
            "structure_id": int(item["structure_id"]),
            "run_id": item["run_id"],
            "source_label": item["source_label"],
            "pdb_id": item["pdb_id"],
            "matched_selectors": [
                selector
                for selector in selectors
                if _normalize_selector(selector) in set(item["matched_selectors"])
            ],
        }
        for item in resolved_structures
    ]
    return result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rerender gallery thumbnails for selected source labels or PDB IDs, "
            "replacing the existing render directories for those matches."
        )
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
    parser.add_argument(
        "--run-id",
        type=str,
        default=None,
        help="Optional run-id filter when the DB contains multiple runs.",
    )
    parser.add_argument("--width", type=int, default=320, help="Thumbnail width in pixels.")
    parser.add_argument("--height", type=int, default=240, help="Thumbnail height in pixels.")
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Parallel structure render worker count (default: 1).",
    )
    parser.add_argument(
        "--node-executable",
        type=str,
        default="node",
        help="Node.js executable used for the Mol* renderer script.",
    )
    parser.add_argument(
        "--renderer-script",
        type=Path,
        default=Path("scripts") / "molstar_render_single.mjs",
        help="Path to the Node Mol* single-structure renderer script.",
    )
    parser.add_argument(
        "--render-backend",
        choices=sorted(gallery_render.VALID_RENDER_BACKENDS),
        default=gallery_render.DEFAULT_RENDER_BACKEND,
        help="Browser rendering backend: software, hardware, or auto (default: software).",
    )
    parser.add_argument(
        "--axis-render-mode",
        choices=sorted(gallery_render.VALID_AXIS_RENDER_MODES),
        default=gallery_render.DEFAULT_AXIS_RENDER_MODE,
        help=(
            "Axis rendering mode: compatibility loads axis-specific atom-filtered "
            "structures into one reused Mol* viewer, fast captures full-structure "
            "views without clipping."
        ),
    )
    parser.add_argument(
        "--worker-timeout-seconds",
        type=float,
        default=0.0,
        help=(
            "Max seconds per external renderer job before it is terminated "
            "(default: 0, disables timeout)."
        ),
    )
    parser.add_argument(
        "--timing-jsonl",
        type=Path,
        default=None,
        help="Optional JSONL path for structured timing output.",
    )
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=gallery_render.DEFAULT_PROGRESS_INTERVAL_SECONDS,
        help=(
            "Seconds between human-readable progress/ETA updates "
            f"(default: {int(gallery_render.DEFAULT_PROGRESS_INTERVAL_SECONDS)}, set <= 0 to disable)."
        ),
    )
    parser.add_argument(
        "selectors",
        nargs="+",
        help="One or more source labels, PDB IDs, or alias PDB IDs to rerender.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    result = rerender_selected_gallery_thumbnails(
        db_path=Path(args.db),
        selectors=[str(selector) for selector in args.selectors],
        render_root=Path(args.render_root),
        run_id=args.run_id,
        width=args.width,
        height=args.height,
        node_executable=args.node_executable,
        renderer_script=Path(args.renderer_script),
        jobs=args.jobs,
        render_backend=args.render_backend,
        axis_render_mode=args.axis_render_mode,
        worker_timeout_seconds=float(args.worker_timeout_seconds),
        timing_jsonl=args.timing_jsonl,
        progress_interval=float(args.progress_interval),
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
