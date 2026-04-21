"""
Thumbnail rendering queue for gallery hits.
"""

from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import time
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


DEFAULT_PROGRESS_INTERVAL_SECONDS = 30.0
DEFAULT_RENDER_BACKEND = "software"
VALID_RENDER_BACKENDS = {"software", "hardware", "auto"}
DEFAULT_AXIS_RENDER_MODE = "compatibility"
VALID_AXIS_RENDER_MODES = {"compatibility", "fast"}
NODE_TIMING_PREFIX = "__VOLUMIZER_TIMING__ "


@dataclass(frozen=True)
class PlannedRenderJob:
    structure_id: int
    run_id: str
    source_label: str
    structure_path: Path
    annotation_path: Path | None
    output_dir: Path
    x_path: Path
    y_path: Path
    z_path: Path


def _adopt_existing_render_outputs(
    *,
    connection: sqlite3.Connection,
    structure_id: int,
    x_path: Path,
    y_path: Path,
    z_path: Path,
    style_hash: str,
) -> None:
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
            datetime.now(tz=timezone.utc).isoformat(),
            int(structure_id),
        ),
    )


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


def _coerce_absolute_path(raw_path: str) -> Path:
    candidate = Path(str(raw_path)).expanduser()
    if candidate.is_absolute():
        return candidate
    return Path.cwd() / candidate


def _warn(message: str) -> None:
    print(f"[gallery-render] {message}", file=sys.stderr)


def _safe_non_negative_int(value: int, name: str) -> int:
    normalized = int(value)
    if normalized < 0:
        raise ValueError(f"{name} must be >= 0")
    return normalized


def _safe_positive_int(value: int, name: str) -> int:
    normalized = int(value)
    if normalized < 1:
        raise ValueError(f"{name} must be >= 1")
    return normalized


def _safe_non_negative_float(value: float, name: str) -> float:
    normalized = float(value)
    if normalized < 0:
        raise ValueError(f"{name} must be >= 0.")
    return normalized


def _validate_choice(value: str, valid_values: set[str], name: str) -> str:
    normalized = str(value).strip().lower()
    if normalized not in valid_values:
        choices = ", ".join(sorted(valid_values))
        raise ValueError(f"{name} must be one of: {choices}.")
    return normalized


def _format_eta_seconds(eta_seconds: float | None) -> str:
    if eta_seconds is None or eta_seconds < 0:
        return "unknown"

    total_seconds = int(round(eta_seconds))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours > 0:
        return f"{hours}:{minutes:02d}:{seconds:02d}"
    return f"{minutes:02d}:{seconds:02d}"


def _emit_render_progress(
    *,
    completed: int,
    total: int,
    failed: int,
    skipped: int,
    started_at: float,
    last_emitted_at: float,
    interval: float,
    force: bool = False,
) -> float:
    if interval <= 0:
        return last_emitted_at

    now = time.monotonic()
    if not force and (now - last_emitted_at) < interval:
        return last_emitted_at

    elapsed = max(0.0, now - started_at)
    rate = completed / elapsed if elapsed > 0 else 0.0
    remaining = max(0, total - completed)
    eta_seconds = 0.0 if remaining == 0 else ((remaining / rate) if rate > 0 else None)
    percent = (100.0 * completed / total) if total > 0 else 100.0
    _warn(
        "progress: "
        f"{completed}/{total} ({percent:.1f}%), "
        f"failed={failed}, skipped={skipped}, rate={rate:.2f}/s, "
        f"eta={_format_eta_seconds(eta_seconds)}"
    )
    return now


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
    annotation_path: Path | None = None,
    output_dir: Path,
    width: int,
    height: int,
    style: dict[str, Any],
    render_backend: str,
    axis_render_mode: str,
    worker_timeout_seconds: float | None = None,
) -> dict[str, Any]:
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
        "--render-backend",
        render_backend,
        "--axis-render-mode",
        axis_render_mode,
    ]
    if annotation_path is not None:
        command.extend(["--annotation", str(annotation_path)])

    timeout_seconds = None
    if worker_timeout_seconds is not None and worker_timeout_seconds > 0:
        timeout_seconds = float(worker_timeout_seconds)

    try:
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        detail = (error.stderr or error.stdout or "").strip()
        message = (
            f"renderer timed out after {float(error.timeout):.1f}s"
            if error.timeout is not None
            else "renderer timed out"
        )
        if detail:
            message = f"{message}: {detail}"
        raise RuntimeError(message) from error

    stdout_text, timing_payload = _extract_node_timing_payload(completed.stdout or "")
    if completed.returncode != 0:
        message = (completed.stderr or stdout_text or "renderer failed").strip()
        raise RuntimeError(message)
    return {
        "stderr": (completed.stderr or "").strip() or None,
        "timing": timing_payload,
    }


def _extract_node_timing_payload(output: str) -> tuple[str, dict[str, Any] | None]:
    if not output:
        return "", None

    remaining_lines: list[str] = []
    timing_payload: dict[str, Any] | None = None
    for raw_line in output.splitlines():
        line = raw_line.strip()
        if line.startswith(NODE_TIMING_PREFIX):
            payload_text = line[len(NODE_TIMING_PREFIX) :]
            try:
                parsed = json.loads(payload_text)
            except json.JSONDecodeError:
                remaining_lines.append(raw_line)
                continue
            if isinstance(parsed, dict):
                timing_payload = parsed
            continue
        remaining_lines.append(raw_line)

    return "\n".join(remaining_lines).strip(), timing_payload


def _append_jsonl(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True))
        handle.write("\n")


def _build_render_signature(
    *,
    style: dict[str, Any],
    render_backend: str,
    axis_render_mode: str,
) -> dict[str, Any]:
    signature = dict(style)
    if render_backend != DEFAULT_RENDER_BACKEND:
        signature["render_backend"] = render_backend
    if axis_render_mode != DEFAULT_AXIS_RENDER_MODE:
        signature["axis_render_mode"] = axis_render_mode
    return signature


def _build_planned_job(
    *,
    row: sqlite3.Row,
    render_root: Path,
) -> PlannedRenderJob:
    source_label = str(row["source_label"])
    run_token = _sanitize_path_token(str(row["run_id"]))
    source_token = _sanitize_path_token(source_label)

    output_dir = render_root / run_token / source_token
    return PlannedRenderJob(
        structure_id=int(row["structure_id"]),
        run_id=str(row["run_id"]),
        source_label=source_label,
        structure_path=_coerce_absolute_path(str(row["annotated_cif_path"])),
        annotation_path=(
            _coerce_absolute_path(str(row["annotation_json_path"]))
            if row["annotation_json_path"]
            else None
        ),
        output_dir=output_dir,
        x_path=output_dir / AXIS_FILENAMES["x"],
        y_path=output_dir / AXIS_FILENAMES["y"],
        z_path=output_dir / AXIS_FILENAMES["z"],
    )


RenderFunction = Callable[[Path, Path, int, int, dict[str, Any]], Any]


def _execute_render_job(
    *,
    job: PlannedRenderJob,
    submitted_at: float,
    width: int,
    height: int,
    style: dict[str, Any],
    render_fn: RenderFunction | None,
    renderer_script: Path,
    node_executable: str,
    render_backend: str,
    axis_render_mode: str,
    worker_timeout_seconds: float | None,
) -> dict[str, Any]:
    started_at = time.monotonic()
    result: dict[str, Any] = {
        "structure_id": job.structure_id,
        "run_id": job.run_id,
        "source_label": job.source_label,
        "x_path": str(job.x_path),
        "y_path": str(job.y_path),
        "z_path": str(job.z_path),
        "dispatch_wait_seconds": max(0.0, started_at - submitted_at),
        "renderer": "custom" if render_fn is not None else "node",
        "status": "failed",
        "error": None,
        "node_result": None,
    }

    try:
        if not job.structure_path.is_file():
            raise FileNotFoundError(f"annotated structure not found: {job.structure_path}")

        job.output_dir.mkdir(parents=True, exist_ok=True)

        if render_fn is not None:
            maybe_result = render_fn(job.structure_path, job.output_dir, width, height, style)
        else:
            maybe_result = _render_single_structure_with_node(
                renderer_script=renderer_script,
                node_executable=node_executable,
                structure_path=job.structure_path,
                annotation_path=job.annotation_path,
                output_dir=job.output_dir,
                width=width,
                height=height,
                style=style,
                render_backend=render_backend,
                axis_render_mode=axis_render_mode,
                worker_timeout_seconds=worker_timeout_seconds,
            )

        for output_path in (job.x_path, job.y_path, job.z_path):
            if not output_path.is_file():
                raise RuntimeError(f"missing output image: {output_path}")

        result["status"] = "done"
        if isinstance(maybe_result, dict):
            result["node_result"] = maybe_result
    except Exception as error:  # pragma: no cover - defensive path
        result["error"] = str(error)
    finally:
        result["wall_seconds"] = max(0.0, time.monotonic() - started_at)

    return result


def render_gallery_thumbnails(
    *,
    db_path: Path,
    render_root: Path,
    run_id: str | None = None,
    limit: int | None = None,
    include_failed: bool = False,
    force: bool = False,
    reuse_existing_only: bool = False,
    width: int = 320,
    height: int = 240,
    node_executable: str = "node",
    renderer_script: Path | None = None,
    style: dict[str, Any] | None = None,
    render_fn: RenderFunction | None = None,
    jobs: int = 1,
    render_backend: str = DEFAULT_RENDER_BACKEND,
    axis_render_mode: str = DEFAULT_AXIS_RENDER_MODE,
    worker_timeout_seconds: float | None = None,
    timing_jsonl: Path | None = None,
    progress_interval: float = DEFAULT_PROGRESS_INTERVAL_SECONDS,
) -> dict[str, Any]:
    db_path = Path(db_path).resolve()
    if not db_path.is_file():
        raise FileNotFoundError(f"Gallery DB does not exist: {db_path}")

    render_root = Path(render_root).resolve()
    render_root.mkdir(parents=True, exist_ok=True)

    width = _safe_non_negative_int(width, "width")
    height = _safe_non_negative_int(height, "height")
    progress_interval = float(progress_interval)
    if progress_interval < 0:
        raise ValueError("progress_interval must be >= 0")
    jobs = _safe_positive_int(jobs, "jobs")
    render_backend = _validate_choice(
        render_backend,
        VALID_RENDER_BACKENDS,
        "render_backend",
    )
    axis_render_mode = _validate_choice(
        axis_render_mode,
        VALID_AXIS_RENDER_MODES,
        "axis_render_mode",
    )
    normalized_worker_timeout: float | None = None
    if worker_timeout_seconds is not None:
        normalized_worker_timeout = _safe_non_negative_float(
            worker_timeout_seconds,
            "worker_timeout_seconds",
        )
        if normalized_worker_timeout <= 0:
            normalized_worker_timeout = None

    normalized_limit = None
    if limit is not None:
        normalized_limit = _safe_non_negative_int(limit, "limit")

    effective_style = dict(DEFAULT_RENDER_STYLE)
    if style is not None:
        effective_style.update(style)
    effective_style["width"] = width
    effective_style["height"] = height

    render_signature = _build_render_signature(
        style=effective_style,
        render_backend=render_backend,
        axis_render_mode=axis_render_mode,
    )
    style_hash = _compute_style_hash(render_signature)
    timing_path = Path(timing_jsonl).resolve() if timing_jsonl is not None else None

    if renderer_script is None:
        renderer_script = Path(__file__).resolve().parents[1] / "scripts" / "molstar_render_single.mjs"

    query = (
        "SELECT "
        "s.structure_id, s.run_id, s.source_label, "
        "s.annotated_cif_path, s.annotation_json_path, "
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
    planning_started_at = time.monotonic()

    total_jobs = 0
    rendered_jobs = 0
    failed_jobs = 0
    skipped_jobs = 0
    reused_jobs = 0
    queue_scan_seconds = 0.0
    dispatch_seconds = 0.0

    with sqlite3.connect(db_path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(query, params).fetchall()
        planned_jobs: list[PlannedRenderJob] = []

        for row in rows:
            if normalized_limit is not None and (len(planned_jobs) + reused_jobs) >= normalized_limit:
                break

            job = _build_planned_job(row=row, render_root=render_root)

            status = str(row["render_status"] or "pending")
            existing_style_hash = row["render_style_hash"]
            style_mismatch = str(existing_style_hash or "") != style_hash
            outputs_exist = job.x_path.is_file() and job.y_path.is_file() and job.z_path.is_file()

            if reuse_existing_only:
                if outputs_exist:
                    _adopt_existing_render_outputs(
                        connection=connection,
                        structure_id=job.structure_id,
                        x_path=job.x_path,
                        y_path=job.y_path,
                        z_path=job.z_path,
                        style_hash=style_hash,
                    )
                    reused_jobs += 1
                else:
                    skipped_jobs += 1
                continue

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

            planned_jobs.append(job)

        total_jobs = len(planned_jobs)
        queue_scan_seconds = max(0.0, time.monotonic() - planning_started_at)
        if timing_path is not None:
            _append_jsonl(
                timing_path,
                {
                    "event": "planning",
                    "db": str(db_path),
                    "render_root": str(render_root),
                    "run_id": run_id,
                    "jobs": jobs,
                    "queue_scan_seconds": queue_scan_seconds,
                    "planned_jobs": total_jobs,
                    "reused_jobs": reused_jobs,
                    "skipped_jobs": skipped_jobs,
                    "render_backend": render_backend,
                    "axis_render_mode": axis_render_mode,
                    "timestamp_utc": datetime.now(tz=timezone.utc).isoformat(),
                },
            )
        started_at = time.monotonic()
        last_progress_emit = started_at
        last_progress_emit = _emit_render_progress(
            completed=0,
            total=total_jobs,
            failed=failed_jobs,
            skipped=skipped_jobs,
            started_at=started_at,
            last_emitted_at=last_progress_emit,
            interval=progress_interval,
            force=True,
        )

        if total_jobs > 0:
            dispatch_started_at = time.monotonic()
            max_workers = min(jobs, total_jobs)
            renderer_path = Path(renderer_script)
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_job: dict[Future[dict[str, Any]], PlannedRenderJob] = {}
                pending_jobs = iter(planned_jobs)

                for _ in range(max_workers):
                    try:
                        job = next(pending_jobs)
                    except StopIteration:
                        break
                    submitted_at = time.monotonic()
                    future = executor.submit(
                        _execute_render_job,
                        job=job,
                        submitted_at=submitted_at,
                        width=width,
                        height=height,
                        style=effective_style,
                        render_fn=render_fn,
                        renderer_script=renderer_path,
                        node_executable=node_executable,
                        render_backend=render_backend,
                        axis_render_mode=axis_render_mode,
                        worker_timeout_seconds=normalized_worker_timeout,
                    )
                    future_to_job[future] = job
                dispatch_seconds = max(0.0, time.monotonic() - dispatch_started_at)

                pending_futures = set(future_to_job)
                while pending_futures:
                    wait_timeout = progress_interval if progress_interval > 0 else None
                    completed_futures, pending_futures = wait(
                        pending_futures,
                        timeout=wait_timeout,
                        return_when=FIRST_COMPLETED,
                    )
                    if not completed_futures:
                        last_progress_emit = _emit_render_progress(
                            completed=rendered_jobs + failed_jobs,
                            total=total_jobs,
                            failed=failed_jobs,
                            skipped=skipped_jobs,
                            started_at=started_at,
                            last_emitted_at=last_progress_emit,
                            interval=progress_interval,
                        )
                        continue

                    for future in completed_futures:
                        job = future_to_job[future]
                        result = future.result()
                        updated_at = datetime.now(tz=timezone.utc).isoformat()
                        node_result = result.get("node_result")
                        node_timing = None
                        if isinstance(node_result, dict):
                            node_timing = node_result.get("timing")

                        if result["status"] == "done":
                            rendered_jobs += 1
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
                                    result["x_path"],
                                    result["y_path"],
                                    result["z_path"],
                                    style_hash,
                                    "done",
                                    updated_at,
                                    job.structure_id,
                                ),
                            )
                        else:
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
                                    result["error"],
                                    updated_at,
                                    job.structure_id,
                                ),
                            )

                        connection.commit()

                        if (
                            isinstance(node_timing, dict)
                            and bool(node_timing.get("fallback_used"))
                        ):
                            _warn(
                                "render backend auto-fallback used for "
                                f"{job.source_label}: "
                                f"{node_timing.get('requested_backend', render_backend)} -> "
                                f"{node_timing.get('backend_used', DEFAULT_RENDER_BACKEND)}"
                            )

                        if timing_path is not None:
                            _append_jsonl(
                                timing_path,
                                {
                                    "event": "job",
                                    "structure_id": job.structure_id,
                                    "run_id": job.run_id,
                                    "source_label": job.source_label,
                                    "status": result["status"],
                                    "error": result["error"],
                                    "renderer": result["renderer"],
                                    "dispatch_wait_seconds": result["dispatch_wait_seconds"],
                                    "wall_seconds": result["wall_seconds"],
                                    "render_backend": render_backend,
                                    "axis_render_mode": axis_render_mode,
                                    "node_timing": node_timing,
                                    "timestamp_utc": updated_at,
                                },
                            )

                        last_progress_emit = _emit_render_progress(
                            completed=rendered_jobs + failed_jobs,
                            total=total_jobs,
                            failed=failed_jobs,
                            skipped=skipped_jobs,
                            started_at=started_at,
                            last_emitted_at=last_progress_emit,
                            interval=progress_interval,
                        )

                        try:
                            next_job = next(pending_jobs)
                        except StopIteration:
                            continue

                        next_submitted_at = time.monotonic()
                        next_future = executor.submit(
                            _execute_render_job,
                            job=next_job,
                            submitted_at=next_submitted_at,
                            width=width,
                            height=height,
                            style=effective_style,
                            render_fn=render_fn,
                            renderer_script=renderer_path,
                            node_executable=node_executable,
                            render_backend=render_backend,
                            axis_render_mode=axis_render_mode,
                            worker_timeout_seconds=normalized_worker_timeout,
                        )
                        future_to_job[next_future] = next_job
                        pending_futures.add(next_future)

    _emit_render_progress(
        completed=rendered_jobs + failed_jobs,
        total=total_jobs,
        failed=failed_jobs,
        skipped=skipped_jobs,
        started_at=started_at,
        last_emitted_at=last_progress_emit,
        interval=progress_interval,
        force=True,
    )

    total_wall_seconds = max(0.0, time.monotonic() - started_at)
    if timing_path is not None:
        _append_jsonl(
            timing_path,
            {
                "event": "summary",
                "db": str(db_path),
                "render_root": str(render_root),
                "run_id": run_id,
                "jobs": jobs,
                "total_jobs": total_jobs,
                "rendered_jobs": rendered_jobs,
                "failed_jobs": failed_jobs,
                "skipped_jobs": skipped_jobs,
                "queue_scan_seconds": queue_scan_seconds,
                "dispatch_seconds": dispatch_seconds,
                "wall_seconds": total_wall_seconds,
                "render_backend": render_backend,
                "axis_render_mode": axis_render_mode,
                "timestamp_utc": datetime.now(tz=timezone.utc).isoformat(),
            },
        )

    return {
        "db": str(db_path),
        "render_root": str(render_root),
        "run_id": run_id,
        "style_hash": style_hash,
        "jobs": jobs,
        "total_jobs": total_jobs,
        "rendered_jobs": rendered_jobs,
        "failed_jobs": failed_jobs,
        "reused_jobs": reused_jobs,
        "skipped_jobs": skipped_jobs,
        "queue_scan_seconds": queue_scan_seconds,
        "dispatch_seconds": dispatch_seconds,
        "wall_seconds": total_wall_seconds,
        "render_backend": render_backend,
        "axis_render_mode": axis_render_mode,
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
    parser.add_argument(
        "--reuse-existing-only",
        action="store_true",
        help=(
            "Do not render new thumbnails. Instead, relink any complete existing "
            "x/y/z PNG triplets already present under --render-root into the DB."
        ),
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
        help="Node.js executable used for Mol* renderer script.",
    )
    parser.add_argument(
        "--renderer-script",
        type=Path,
        default=Path("scripts") / "molstar_render_single.mjs",
        help="Path to Node Mol* single-structure renderer script.",
    )
    parser.add_argument(
        "--render-backend",
        choices=sorted(VALID_RENDER_BACKENDS),
        default=DEFAULT_RENDER_BACKEND,
        help="Browser rendering backend: software, hardware, or auto (default: software).",
    )
    parser.add_argument(
        "--axis-render-mode",
        choices=sorted(VALID_AXIS_RENDER_MODES),
        default=DEFAULT_AXIS_RENDER_MODE,
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
        help="Optional JSONL path for planning, per-job, and renderer timing records.",
    )
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=DEFAULT_PROGRESS_INTERVAL_SECONDS,
        help=(
            "Seconds between human-readable progress/ETA updates "
            f"(default: {int(DEFAULT_PROGRESS_INTERVAL_SECONDS)}, set <= 0 to disable)."
        ),
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
        reuse_existing_only=bool(args.reuse_existing_only),
        width=args.width,
        height=args.height,
        jobs=args.jobs,
        node_executable=args.node_executable,
        renderer_script=Path(args.renderer_script),
        render_backend=args.render_backend,
        axis_render_mode=args.axis_render_mode,
        worker_timeout_seconds=args.worker_timeout_seconds,
        timing_jsonl=args.timing_jsonl,
        progress_interval=float(args.progress_interval),
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
