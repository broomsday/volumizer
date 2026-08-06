"""
FastAPI app for local browsing of indexed gallery results.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any

from fastapi import FastAPI, File, Form, HTTPException, Query, UploadFile
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles

from volumizer import gallery_query
from volumizer.molstar import build_volume_surface_style
from volumizer import pdb
from volumizer.web import db as web_db
from volumizer.web.jobs import (
    ALLOWED_UPLOAD_EXTENSIONS,
    DEFAULT_MAX_UPLOAD_BYTES,
    AnalysisJobRegistry,
)


DEFAULT_DB_PATH = Path("data") / "gallery.db"
GALLERY_DB_ENV = "VOLUMIZER_GALLERY_DB"
MOLSTAR_ASSET_ROOT_ENV = "MOLSTAR_ASSET_ROOT"
ENABLE_UPLOAD_ENV = "VOLUMIZER_ENABLE_UPLOAD"
MAX_UPLOAD_BYTES_ENV = "VOLUMIZER_UPLOAD_MAX_BYTES"
DEFAULT_MOLSTAR_ASSET_ROOT = Path("node_modules") / "molstar" / "build" / "viewer"
MOLSTAR_ASSET_FILENAMES = frozenset({"molstar.js", "molstar.css"})


def _resolve_db_path(explicit_db_path: Path | None = None) -> Path:
    if explicit_db_path is not None:
        return Path(explicit_db_path).expanduser().resolve()

    raw_env_path = os.getenv(GALLERY_DB_ENV)
    if raw_env_path:
        return Path(raw_env_path).expanduser().resolve()

    return (Path.cwd() / DEFAULT_DB_PATH).resolve()


def _resolve_molstar_asset_root(explicit_asset_root: Path | None = None) -> Path | None:
    candidates: list[Path] = []
    if explicit_asset_root is not None:
        candidates.append(Path(explicit_asset_root))

    raw_env_path = os.getenv(MOLSTAR_ASSET_ROOT_ENV)
    if raw_env_path:
        candidates.append(Path(raw_env_path))

    package_root = Path(__file__).resolve().parents[2]
    candidates.append(package_root / DEFAULT_MOLSTAR_ASSET_ROOT)
    candidates.append(Path(__file__).with_name("static") / "vendor" / "molstar")

    for candidate in candidates:
        resolved = candidate.expanduser().resolve()
        if all((resolved / filename).is_file() for filename in MOLSTAR_ASSET_FILENAMES):
            return resolved
    return None


def _resolve_upload_enabled(explicit_upload_enabled: bool | None = None) -> bool:
    if explicit_upload_enabled is not None:
        return bool(explicit_upload_enabled)

    raw_value = os.getenv(ENABLE_UPLOAD_ENV)
    if raw_value is None:
        return True

    normalized = raw_value.strip().lower()
    return normalized not in {"0", "false", "no", "off"}


def _resolve_max_upload_bytes(explicit_max_upload_bytes: int | None = None) -> int:
    if explicit_max_upload_bytes is not None:
        return int(explicit_max_upload_bytes)

    raw_value = os.getenv(MAX_UPLOAD_BYTES_ENV)
    if raw_value:
        try:
            return int(raw_value)
        except ValueError:
            return DEFAULT_MAX_UPLOAD_BYTES

    return DEFAULT_MAX_UPLOAD_BYTES


def _load_index_html(static_dir: Path) -> str:
    return (static_dir / "index.html").read_text(encoding="utf-8")


def _get_detail_or_404(db_path: Path, structure_id: int) -> dict[str, Any]:
    try:
        detail = web_db.get_hit_detail(db_path, structure_id)
    except FileNotFoundError as error:
        raise HTTPException(status_code=404, detail=str(error)) from error

    if detail is None:
        raise HTTPException(status_code=404, detail=f"Unknown structure_id: {structure_id}")
    return detail


def _get_artifacts_or_404(db_path: Path, structure_id: int) -> dict[str, Any]:
    try:
        artifacts = web_db.get_hit_artifacts(db_path, structure_id)
    except FileNotFoundError as error:
        raise HTTPException(status_code=404, detail=str(error)) from error

    if artifacts is None:
        raise HTTPException(status_code=404, detail=f"Unknown structure_id: {structure_id}")
    return artifacts


def _decorate_record(record: dict[str, Any]) -> dict[str, Any]:
    structure_id = int(record["structure_id"])
    record_with_urls = dict(record)
    record_with_urls["urls"] = {
        "detail": f"/api/hits/{structure_id}",
        "viewer_data": f"/api/hits/{structure_id}/viewer-data",
        "structure": (
            f"/files/structure/{structure_id}"
            if record.get("annotated_cif_path")
            else None
        ),
        "annotation": (
            f"/files/annotation/{structure_id}"
            if record.get("annotation_json_path")
            else None
        ),
        "thumbnail_x": (
            f"/files/render/{structure_id}/x.png"
            if record.get("x_png_path")
            else None
        ),
        "thumbnail_y": (
            f"/files/render/{structure_id}/y.png"
            if record.get("y_png_path")
            else None
        ),
        "thumbnail_z": (
            f"/files/render/{structure_id}/z.png"
            if record.get("z_png_path")
            else None
        ),
    }
    return record_with_urls


def _file_response_or_404(
    *,
    raw_path: str | None,
    description: str,
    media_type: str,
) -> FileResponse:
    path = web_db.resolve_artifact_path(raw_path)
    if path is None or not path.is_file():
        raise HTTPException(status_code=404, detail=f"{description} not found")

    return FileResponse(path, media_type=media_type, filename=path.name)


def _molstar_asset_response_or_404(
    *,
    asset_root: Path | None,
    filename: str,
    media_type: str,
) -> FileResponse:
    if filename not in MOLSTAR_ASSET_FILENAMES:
        raise HTTPException(status_code=404, detail=f"Unsupported Mol* asset: {filename}")

    if asset_root is None:
        raise HTTPException(
            status_code=404,
            detail=(
                "Mol* assets are unavailable. Install gallery dependencies with "
                "`npm ci` and `npm run gallery:install-browser`."
            ),
        )

    path = asset_root / filename
    if not path.is_file():
        raise HTTPException(status_code=404, detail=f"Mol* asset not found: {filename}")

    return FileResponse(path, media_type=media_type, filename=path.name)


def create_app(
    db_path: Path | None = None,
    molstar_asset_root: Path | None = None,
    upload_enabled: bool | None = None,
    job_registry: AnalysisJobRegistry | None = None,
    max_upload_bytes: int | None = None,
) -> FastAPI:
    resolved_db_path = _resolve_db_path(db_path)
    resolved_upload_enabled = _resolve_upload_enabled(upload_enabled)
    resolved_max_upload_bytes = _resolve_max_upload_bytes(max_upload_bytes)
    static_dir = Path(__file__).with_name("static")
    index_html = _load_index_html(static_dir)

    app = FastAPI(
        title="Volumizer Gallery",
        docs_url="/api/docs",
        redoc_url=None,
    )
    app.state.db_path = resolved_db_path
    app.state.molstar_asset_root = _resolve_molstar_asset_root(molstar_asset_root)
    app.state.upload_enabled = resolved_upload_enabled
    app.state.max_upload_bytes = resolved_max_upload_bytes
    app.state.analysis_jobs = job_registry or (
        AnalysisJobRegistry(
            db_path=resolved_db_path,
            max_upload_bytes=resolved_max_upload_bytes,
        )
        if resolved_upload_enabled
        else None
    )

    app.mount("/static", StaticFiles(directory=static_dir), name="static")

    @app.get("/", response_class=HTMLResponse)
    def serve_index() -> HTMLResponse:
        return HTMLResponse(index_html)

    @app.get("/api/health")
    def health() -> dict[str, Any]:
        return {
            "status": "ok",
            "db": str(app.state.db_path),
            "db_exists": Path(app.state.db_path).is_file(),
            "molstar_asset_root": (
                str(app.state.molstar_asset_root)
                if app.state.molstar_asset_root is not None
                else None
            ),
            "molstar_assets_available": app.state.molstar_asset_root is not None,
            "upload_enabled": app.state.upload_enabled,
            "max_upload_bytes": app.state.max_upload_bytes,
            "assembly_policies": list(pdb.VALID_ASSEMBLY_POLICIES),
            "renderer_available": app.state.analysis_jobs is not None,
        }

    @app.post("/api/analyze")
    async def submit_analysis(
        file: UploadFile = File(...),
        resolution: float = Form(...),
        assembly_policy: str = Form(pdb.DEFAULT_ASSEMBLY_POLICY),
    ) -> dict[str, Any]:
        if not app.state.upload_enabled or app.state.analysis_jobs is None:
            raise HTTPException(status_code=403, detail="Upload analysis is disabled")

        filename = Path(file.filename or "").name
        if len(filename.strip()) == 0:
            raise HTTPException(status_code=422, detail="Upload filename is required")

        suffix = Path(filename).suffix.lower()
        if suffix not in ALLOWED_UPLOAD_EXTENSIONS:
            raise HTTPException(
                status_code=422,
                detail=(
                    "Unsupported upload type. Expected one of: "
                    + ", ".join(sorted(ALLOWED_UPLOAD_EXTENSIONS))
                ),
            )

        if float(resolution) <= 0:
            raise HTTPException(status_code=422, detail="resolution must be greater than 0")

        if assembly_policy not in pdb.VALID_ASSEMBLY_POLICIES:
            raise HTTPException(
                status_code=422,
                detail=(
                    "Unsupported assembly_policy. Expected one of: "
                    + ", ".join(pdb.VALID_ASSEMBLY_POLICIES)
                ),
            )

        content = await file.read(int(app.state.max_upload_bytes) + 1)
        if len(content) == 0:
            raise HTTPException(status_code=422, detail="Uploaded file is empty")
        if len(content) > int(app.state.max_upload_bytes):
            raise HTTPException(status_code=413, detail="Uploaded file is too large")

        job = app.state.analysis_jobs.submit_analysis(
            filename=filename,
            content=content,
            resolution=float(resolution),
            assembly_policy=assembly_policy,
        )
        return {
            "job_id": job.job_id,
            "status": job.status,
        }

    @app.get("/api/analyze")
    def list_analysis_jobs() -> dict[str, Any]:
        registry = app.state.analysis_jobs
        if registry is None:
            return {"jobs": []}
        return {"jobs": [job.to_dict() for job in registry.list_jobs()]}

    @app.get("/api/analyze/{job_id}")
    def get_analysis_job(job_id: str) -> dict[str, Any]:
        registry = app.state.analysis_jobs
        if registry is None:
            raise HTTPException(status_code=404, detail=f"Unknown job_id: {job_id}")
        job = registry.get_job(job_id)
        if job is None:
            raise HTTPException(status_code=404, detail=f"Unknown job_id: {job_id}")
        return job.to_dict()

    @app.get("/assets/molstar.js")
    def get_molstar_js() -> FileResponse:
        return _molstar_asset_response_or_404(
            asset_root=app.state.molstar_asset_root,
            filename="molstar.js",
            media_type="application/javascript",
        )

    @app.get("/assets/molstar.css")
    def get_molstar_css() -> FileResponse:
        return _molstar_asset_response_or_404(
            asset_root=app.state.molstar_asset_root,
            filename="molstar.css",
            media_type="text/css",
        )

    @app.get("/api/runs")
    def list_runs() -> dict[str, Any]:
        try:
            runs = web_db.list_runs(app.state.db_path)
        except FileNotFoundError as error:
            raise HTTPException(status_code=404, detail=str(error)) from error

        return {
            "db": str(app.state.db_path),
            "runs": runs,
        }

    @app.get("/api/hits")
    def list_hits(
        run_id: str | None = None,
        pdb_id_query: str | None = None,
        pore_volume_min: float | None = None,
        pore_volume_max: float | None = None,
        pore_length_min: float | None = None,
        pore_length_max: float | None = None,
        pore_dmin_min: float | None = None,
        pore_dmin_max: float | None = None,
        pore_dmax_min: float | None = None,
        pore_dmax_max: float | None = None,
        num_pores_min: int | None = None,
        num_pores_max: int | None = None,
        pocket_volume_min: float | None = None,
        pocket_volume_max: float | None = None,
        pocket_length_min: float | None = None,
        pocket_length_max: float | None = None,
        pocket_dmin_min: float | None = None,
        pocket_dmin_max: float | None = None,
        pocket_dmax_min: float | None = None,
        pocket_dmax_max: float | None = None,
        num_pockets_min: int | None = None,
        num_pockets_max: int | None = None,
        cavity_volume_min: float | None = None,
        cavity_volume_max: float | None = None,
        cavity_length_min: float | None = None,
        cavity_length_max: float | None = None,
        cavity_dmin_min: float | None = None,
        cavity_dmin_max: float | None = None,
        cavity_dmax_min: float | None = None,
        cavity_dmax_max: float | None = None,
        num_cavities_min: int | None = None,
        num_cavities_max: int | None = None,
        chains_min: int | None = None,
        chains_max: int | None = None,
        residues_min: int | None = None,
        residues_max: int | None = None,
        seq_unique_chains_min: int | None = None,
        seq_unique_chains_max: int | None = None,
        frac_alpha_min: float | None = None,
        frac_alpha_max: float | None = None,
        frac_beta_min: float | None = None,
        frac_beta_max: float | None = None,
        frac_coil_min: float | None = None,
        frac_coil_max: float | None = None,
        pore_circularity_min: float | None = None,
        pore_circularity_max: float | None = None,
        pore_uniformity_min: float | None = None,
        pore_uniformity_max: float | None = None,
        pocket_circularity_min: float | None = None,
        pocket_circularity_max: float | None = None,
        pocket_uniformity_min: float | None = None,
        pocket_uniformity_max: float | None = None,
        cavity_circularity_min: float | None = None,
        cavity_circularity_max: float | None = None,
        cavity_uniformity_min: float | None = None,
        cavity_uniformity_max: float | None = None,
        limit: int = Query(default=24, ge=0),
        offset: int = Query(default=0, ge=0),
        sort_by: str = "largest_pore_volume",
        sort_dir: str = "desc",
    ) -> dict[str, Any]:
        try:
            result = gallery_query.query_gallery_index(
                db_path=app.state.db_path,
                run_id=run_id,
                pdb_id_query=pdb_id_query,
                pore_volume_min=pore_volume_min,
                pore_volume_max=pore_volume_max,
                pore_length_min=pore_length_min,
                pore_length_max=pore_length_max,
                pore_dmin_min=pore_dmin_min,
                pore_dmin_max=pore_dmin_max,
                pore_dmax_min=pore_dmax_min,
                pore_dmax_max=pore_dmax_max,
                num_pores_min=num_pores_min,
                num_pores_max=num_pores_max,
                pocket_volume_min=pocket_volume_min,
                pocket_volume_max=pocket_volume_max,
                pocket_length_min=pocket_length_min,
                pocket_length_max=pocket_length_max,
                pocket_dmin_min=pocket_dmin_min,
                pocket_dmin_max=pocket_dmin_max,
                pocket_dmax_min=pocket_dmax_min,
                pocket_dmax_max=pocket_dmax_max,
                num_pockets_min=num_pockets_min,
                num_pockets_max=num_pockets_max,
                cavity_volume_min=cavity_volume_min,
                cavity_volume_max=cavity_volume_max,
                cavity_length_min=cavity_length_min,
                cavity_length_max=cavity_length_max,
                cavity_dmin_min=cavity_dmin_min,
                cavity_dmin_max=cavity_dmin_max,
                cavity_dmax_min=cavity_dmax_min,
                cavity_dmax_max=cavity_dmax_max,
                num_cavities_min=num_cavities_min,
                num_cavities_max=num_cavities_max,
                chains_min=chains_min,
                chains_max=chains_max,
                residues_min=residues_min,
                residues_max=residues_max,
                seq_unique_chains_min=seq_unique_chains_min,
                seq_unique_chains_max=seq_unique_chains_max,
                frac_alpha_min=frac_alpha_min,
                frac_alpha_max=frac_alpha_max,
                frac_beta_min=frac_beta_min,
                frac_beta_max=frac_beta_max,
                frac_coil_min=frac_coil_min,
                frac_coil_max=frac_coil_max,
                pore_circularity_min=pore_circularity_min,
                pore_circularity_max=pore_circularity_max,
                pore_uniformity_min=pore_uniformity_min,
                pore_uniformity_max=pore_uniformity_max,
                pocket_circularity_min=pocket_circularity_min,
                pocket_circularity_max=pocket_circularity_max,
                pocket_uniformity_min=pocket_uniformity_min,
                pocket_uniformity_max=pocket_uniformity_max,
                cavity_circularity_min=cavity_circularity_min,
                cavity_circularity_max=cavity_circularity_max,
                cavity_uniformity_min=cavity_uniformity_min,
                cavity_uniformity_max=cavity_uniformity_max,
                limit=limit,
                offset=offset,
                sort_by=sort_by,
                sort_dir=sort_dir,
            )
        except FileNotFoundError as error:
            raise HTTPException(status_code=404, detail=str(error)) from error
        except ValueError as error:
            raise HTTPException(status_code=422, detail=str(error)) from error

        result["rows"] = [_decorate_record(row) for row in result["rows"]]
        return result

    @app.get("/api/hits/{structure_id}")
    def get_hit(structure_id: int) -> dict[str, Any]:
        return _decorate_record(_get_detail_or_404(app.state.db_path, structure_id))

    @app.get("/api/hits/{structure_id}/viewer-data")
    def get_viewer_data(structure_id: int) -> dict[str, Any]:
        detail = _get_detail_or_404(app.state.db_path, structure_id)
        structure_path = web_db.resolve_artifact_path(detail.get("annotated_cif_path"))
        voxel_resolution = detail.get("run_resolution")

        return {
            "structure_id": structure_id,
            "source_label": detail.get("source_label"),
            "voxel_resolution": voxel_resolution,
            "volume_surface": build_volume_surface_style(voxel_resolution),
            "structure_format": web_db.detect_structure_format(structure_path),
            "structure_url": (
                f"/files/structure/{structure_id}"
                if detail.get("annotated_cif_path")
                else None
            ),
            "annotation_url": (
                f"/files/annotation/{structure_id}"
                if detail.get("annotation_json_path")
                else None
            ),
            "thumbnails": {
                "x": (
                    f"/files/render/{structure_id}/x.png"
                    if detail.get("x_png_path")
                    else None
                ),
                "y": (
                    f"/files/render/{structure_id}/y.png"
                    if detail.get("y_png_path")
                    else None
                ),
                "z": (
                    f"/files/render/{structure_id}/z.png"
                    if detail.get("z_png_path")
                    else None
                ),
            },
        }

    @app.get("/files/structure/{structure_id}")
    def get_structure_file(structure_id: int) -> FileResponse:
        artifacts = _get_artifacts_or_404(app.state.db_path, structure_id)
        structure_path = web_db.resolve_artifact_path(artifacts.get("annotated_cif_path"))
        structure_format = web_db.detect_structure_format(structure_path)
        media_type = (
            "chemical/x-pdb" if structure_format == "pdb" else "chemical/x-cif"
        )
        return _file_response_or_404(
            raw_path=artifacts.get("annotated_cif_path"),
            description="Annotated structure",
            media_type=media_type,
        )

    @app.get("/files/annotation/{structure_id}")
    def get_annotation_file(structure_id: int) -> FileResponse:
        artifacts = _get_artifacts_or_404(app.state.db_path, structure_id)
        return _file_response_or_404(
            raw_path=artifacts.get("annotation_json_path"),
            description="Annotation JSON",
            media_type="application/json",
        )

    @app.get("/files/render/{structure_id}/{axis}.png")
    def get_render_file(structure_id: int, axis: str) -> FileResponse:
        if axis not in web_db.VALID_RENDER_AXES:
            raise HTTPException(status_code=404, detail=f"Unsupported render axis: {axis}")

        artifacts = _get_artifacts_or_404(app.state.db_path, structure_id)
        raw_path = artifacts.get(f"{axis}_png_path")
        return _file_response_or_404(
            raw_path=raw_path,
            description=f"{axis.upper()} thumbnail",
            media_type="image/png",
        )

    return app


def create_default_app() -> FastAPI:
    return create_app()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Serve a local web UI for browsing the volumizer gallery index."
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=DEFAULT_DB_PATH,
        help="SQLite gallery DB path (default: data/gallery.db).",
    )
    parser.add_argument(
        "--host",
        type=str,
        default="127.0.0.1",
        help="Host interface to bind (default: 127.0.0.1).",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Port to bind (default: 8000).",
    )
    parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable auto-reload during local development.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    os.environ[GALLERY_DB_ENV] = str(_resolve_db_path(args.db))

    import uvicorn

    if args.reload:
        uvicorn.run(
            "volumizer.web.app:create_default_app",
            factory=True,
            host=args.host,
            port=args.port,
            reload=True,
        )
    else:
        uvicorn.run(
            create_default_app(),
            host=args.host,
            port=args.port,
            reload=False,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
