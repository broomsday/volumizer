"""
In-process upload analysis jobs for the local gallery web app.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
import queue
import re
import threading
import traceback
from typing import Any, Callable
from uuid import uuid4

from volumizer import gallery_index, gallery_render, utils
from volumizer.cli import (
    PostAssemblyResidueLimitExceeded,
    _run_isolated_analysis_worker,
)
from volumizer.constants import MIN_NUM_VOXELS
from volumizer.pdb import DEFAULT_ASSEMBLY_POLICY
from volumizer.web import db as web_db


UPLOAD_RUN_ID = "uploads"
ALLOWED_UPLOAD_EXTENSIONS = frozenset({".pdb", ".cif", ".mmcif"})
DEFAULT_MAX_UPLOAD_BYTES = 25 * 1024 * 1024
DEFAULT_MAX_RESIDUES = 100_000

JOB_QUEUED = "queued"
JOB_RUNNING = "running"
JOB_INDEXING = "indexing"
JOB_DONE = "done"
JOB_ERROR = "error"

THUMBNAIL_PENDING = "pending"
THUMBNAIL_RENDERING = "rendering"
THUMBNAIL_DONE = "done"
THUMBNAIL_FAILED = "failed"
THUMBNAIL_SKIPPED = "skipped"


AnalysisFunction = Callable[..., dict[str, Any]]
ThumbnailFunction = Callable[..., dict[str, Any] | None]


@dataclass
class Job:
    job_id: str
    status: str
    message: str
    structure_id: int | None
    error: str | None
    created_at: str
    updated_at: str
    thumbnail_status: str
    source_label: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class AnalysisRequest:
    filename: str
    content: bytes
    resolution: float
    assembly_policy: str


def _now_utc() -> str:
    return datetime.now(tz=timezone.utc).isoformat()


def _sanitize_source_label(filename: str) -> str:
    stem = Path(filename).stem.strip()
    candidate = re.sub(r"[^A-Za-z0-9_.-]+", "_", stem).strip("._-")
    if len(candidate) == 0:
        return "upload"
    return candidate[:80]


def _dedupe_source_label(base_label: str, existing_labels: set[str]) -> str:
    if base_label not in existing_labels:
        return base_label

    for index in range(2, 10_000):
        candidate = f"{base_label}-{index}"
        if candidate not in existing_labels:
            return candidate

    return f"{base_label}-{uuid4().hex[:8]}"


def _paths_from_analysis_result(
    *,
    result: dict[str, Any],
    source_label: str,
    output_dir: Path,
) -> tuple[Path, Path]:
    annotated_path = result.get("structure_output")
    annotation_path = result.get("annotation_output")
    if annotated_path is None:
        annotated_path = output_dir / f"{source_label}.annotated.cif"
    if annotation_path is None:
        annotation_path = output_dir / f"{source_label}.annotation.json"
    return Path(annotated_path), Path(annotation_path)


def default_analysis_function(
    *,
    source_label: str,
    input_path: Path,
    output_dir: Path,
    resolution: float,
    assembly_policy: str = DEFAULT_ASSEMBLY_POLICY,
    max_residues: int | None = DEFAULT_MAX_RESIDUES,
) -> dict[str, Any]:
    return _run_isolated_analysis_worker(
        source_label=source_label,
        input_path=input_path,
        output_dir=output_dir,
        min_voxels=MIN_NUM_VOXELS,
        min_volume=None,
        overwrite=True,
        assembly_policy=assembly_policy,
        resolution=resolution,
        keep_non_protein=False,
        backend=None,
        surface_connectivity=utils.get_surface_component_connectivity_mode(),
        merge_mouth_gap_voxels=utils.get_surface_mouth_merge_gap_voxels(),
        max_residues=max_residues,
        include_hubs=False,
        enable_necked_pocket_cavity=True,
    )


def default_thumbnail_function(
    *,
    db_path: Path,
    render_root: Path,
    structure_id: int,
) -> dict[str, Any]:
    return gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        structure_ids=[int(structure_id)],
        limit=1,
    )


class AnalysisJobRegistry:
    def __init__(
        self,
        *,
        db_path: Path,
        analysis_fn: AnalysisFunction | None = None,
        thumbnail_fn: ThumbnailFunction | None = None,
        max_upload_bytes: int = DEFAULT_MAX_UPLOAD_BYTES,
        max_residues: int | None = DEFAULT_MAX_RESIDUES,
    ) -> None:
        self.db_path = Path(db_path).expanduser().resolve()
        self.upload_root = self.db_path.parent / "runs" / UPLOAD_RUN_ID
        self.incoming_root = self.upload_root / "incoming"
        self.render_root = self.db_path.parent / "renders"
        self.analysis_fn = analysis_fn or default_analysis_function
        self.thumbnail_fn = thumbnail_fn or default_thumbnail_function
        self.max_upload_bytes = int(max_upload_bytes)
        self.max_residues = max_residues
        self._jobs: dict[str, Job] = {}
        self._requests: dict[str, AnalysisRequest] = {}
        self._lock = threading.Lock()
        self._queue: queue.Queue[str] = queue.Queue()
        self._worker: threading.Thread | None = None

    def submit_analysis(
        self,
        *,
        filename: str,
        content: bytes,
        resolution: float,
        assembly_policy: str,
    ) -> Job:
        job_id = uuid4().hex
        now = _now_utc()
        job = Job(
            job_id=job_id,
            status=JOB_QUEUED,
            message="Queued for analysis.",
            structure_id=None,
            error=None,
            created_at=now,
            updated_at=now,
            thumbnail_status=THUMBNAIL_PENDING,
        )
        request = AnalysisRequest(
            filename=filename,
            content=bytes(content),
            resolution=float(resolution),
            assembly_policy=str(assembly_policy),
        )
        with self._lock:
            self._jobs[job_id] = job
            self._requests[job_id] = request
        self._ensure_worker()
        self._queue.put(job_id)
        return job

    def get_job(self, job_id: str) -> Job | None:
        with self._lock:
            return self._jobs.get(str(job_id))

    def list_jobs(self) -> list[Job]:
        with self._lock:
            return sorted(
                self._jobs.values(),
                key=lambda job: job.created_at,
                reverse=True,
            )

    def _set_job(self, job_id: str, **updates: Any) -> Job | None:
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                return None
            for key, value in updates.items():
                setattr(job, key, value)
            job.updated_at = _now_utc()
            return job

    def _ensure_worker(self) -> None:
        with self._lock:
            if self._worker is not None and self._worker.is_alive():
                return
            self._worker = threading.Thread(
                target=self._worker_loop,
                name="volumizer-upload-analysis",
                daemon=True,
            )
            self._worker.start()

    def _worker_loop(self) -> None:
        while True:
            job_id = self._queue.get()
            try:
                self._run_job(job_id)
            finally:
                self._queue.task_done()

    def _run_job(self, job_id: str) -> None:
        with self._lock:
            request = self._requests.get(job_id)
        if request is None:
            return

        try:
            self.incoming_root.mkdir(parents=True, exist_ok=True)
            base_label = _sanitize_source_label(request.filename)
            existing_labels = set(web_db.list_source_labels(self.db_path, UPLOAD_RUN_ID))
            source_label = _dedupe_source_label(base_label, existing_labels)
            suffix = Path(request.filename).suffix.lower()
            input_path = self.incoming_root / f"{job_id}{suffix}"
            input_path.write_bytes(request.content)
            output_dir = self.upload_root / source_label

            self._set_job(
                job_id,
                status=JOB_RUNNING,
                message="Running volumizer analysis.",
                source_label=source_label,
            )
            result = self.analysis_fn(
                source_label=source_label,
                input_path=input_path,
                output_dir=output_dir,
                resolution=request.resolution,
                assembly_policy=request.assembly_policy,
                max_residues=self.max_residues,
            )

            annotated_cif_path, annotation_json_path = _paths_from_analysis_result(
                result=result,
                source_label=source_label,
                output_dir=output_dir,
            )
            self._set_job(
                job_id,
                status=JOB_INDEXING,
                message="Indexing analyzed structure.",
            )
            structure_id = gallery_index.index_single_structure(
                db_path=self.db_path,
                run_id=UPLOAD_RUN_ID,
                source_label=source_label,
                input_path=input_path,
                annotated_cif_path=annotated_cif_path,
                annotation_json_path=annotation_json_path,
                resolution=request.resolution,
                assembly_policy=request.assembly_policy,
            )

            self._set_job(
                job_id,
                status=JOB_DONE,
                message="Analysis complete.",
                structure_id=structure_id,
            )
            self._render_thumbnail(job_id=job_id, structure_id=structure_id)
        except PostAssemblyResidueLimitExceeded as error:
            self._set_job(
                job_id,
                status=JOB_ERROR,
                message="Uploaded structure is too large after assembly.",
                error=(
                    f"Structure has {error.actual_residues} residues after "
                    f"{error.assembly_policy} assembly; limit is {error.max_residues}."
                ),
                thumbnail_status=THUMBNAIL_SKIPPED,
            )
        except Exception as error:
            traceback.print_exc()
            self._set_job(
                job_id,
                status=JOB_ERROR,
                message="Analysis failed.",
                error=str(error),
                thumbnail_status=THUMBNAIL_SKIPPED,
            )
        finally:
            with self._lock:
                self._requests.pop(job_id, None)

    def _render_thumbnail(self, *, job_id: str, structure_id: int) -> None:
        if self.thumbnail_fn is None:
            self._set_job(job_id, thumbnail_status=THUMBNAIL_SKIPPED)
            return

        self._set_job(
            job_id,
            thumbnail_status=THUMBNAIL_RENDERING,
            message="Analysis complete; rendering thumbnail.",
        )
        try:
            self.thumbnail_fn(
                db_path=self.db_path,
                render_root=self.render_root,
                structure_id=structure_id,
            )
        except Exception:
            traceback.print_exc()
            self._set_job(
                job_id,
                thumbnail_status=THUMBNAIL_FAILED,
                message="Analysis complete; thumbnail rendering failed.",
            )
            return

        self._set_job(
            job_id,
            thumbnail_status=THUMBNAIL_DONE,
            message="Analysis complete.",
        )
