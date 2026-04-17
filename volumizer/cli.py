"""
Command-line interface for volumizer.
"""

from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, as_completed, wait
from datetime import datetime, timezone
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version as package_version
import inspect
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
import traceback
from typing import Sequence

from volumizer import rcsb
from volumizer.constants import (
    DIRECT_SURFACE_MOUTH_MERGE_GAP_VOXELS,
    VALID_DIRECT_SURFACE_COMPONENT_CONNECTIVITY_MODES,
)


DEFAULT_METADATA_STORE_DIRNAME = "entry_metadata"
METADATA_RECORD_FORMAT_VERSION = 1
STATUS_RECORD_FORMAT_VERSION = 1
DEFAULT_CHECKPOINT_FILENAME = "run.checkpoint.json"
MANIFEST_FORMAT_VERSION = 1
BACKEND_ENV = "VOLUMIZER_BACKEND"
BACKEND_AUTO = "auto"
BACKEND_NATIVE = "native"
BACKEND_PYTHON = "python"
VALID_BACKENDS = (BACKEND_AUTO, BACKEND_NATIVE, BACKEND_PYTHON)
DEFAULT_ASSEMBLY_POLICY = "biological"
DEFAULT_CLUSTER_MAX_RESIDUES = 20000
VALID_ASSEMBLY_POLICIES = (
    "biological",
    "asymmetric",
    "auto",
)
_UNSET = object()
_ANALYSIS_WORKER_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "RAYON_NUM_THREADS",
)
DEFAULT_ANALYSIS_WORKER_TIMEOUT_SECONDS = 1800.0


class PostAssemblyResidueLimitExceeded(RuntimeError):
    """
    Runtime size guard for assemblies that expand well beyond deposited-entry metadata.
    """

    def __init__(
        self,
        *,
        actual_residues: int,
        max_residues: int,
        assembly_policy: str,
    ):
        self.actual_residues = int(actual_residues)
        self.max_residues = int(max_residues)
        self.assembly_policy = str(assembly_policy)
        super().__init__(
            "Structure exceeds post-assembly residue limit before annotation: "
            f"{self.actual_residues} > {self.max_residues} "
            f"(assembly_policy={self.assembly_policy})"
        )


@lru_cache(maxsize=1)
def _native_backend_module():
    from volumizer import native_backend

    return native_backend


@lru_cache(maxsize=1)
def _pdb_module():
    from volumizer import pdb

    return pdb


@lru_cache(maxsize=1)
def _utils_module():
    from volumizer import utils

    return utils


@lru_cache(maxsize=1)
def _volumizer_module():
    from volumizer import volumizer

    return volumizer


def _requested_backend_label() -> str:
    backend = os.getenv(BACKEND_ENV, BACKEND_AUTO).strip().lower()
    if backend not in VALID_BACKENDS:
        return BACKEND_AUTO
    return backend


def _get_post_assembly_max_residues(args: argparse.Namespace) -> int | None:
    if getattr(args, "command", None) != "cluster":
        return None

    max_residues = getattr(args, "cluster_max_residues", None)
    if max_residues is None:
        return None
    return int(max_residues)


def _get_analysis_worker_timeout_seconds(
    args: argparse.Namespace,
    *,
    use_isolated_workers: bool,
) -> float | None:
    if not use_isolated_workers:
        return None

    configured_timeout = getattr(args, "worker_timeout_seconds", None)
    if configured_timeout is None:
        return DEFAULT_ANALYSIS_WORKER_TIMEOUT_SECONDS
    if configured_timeout <= 0:
        return None
    return float(configured_timeout)


def _resolve_cli_version() -> str:
    try:
        return package_version("volumizer")
    except PackageNotFoundError:
        return "0+unknown"


def _utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _format_eta_seconds(eta_seconds: float | None) -> str:
    if eta_seconds is None or eta_seconds < 0:
        return "unknown"

    total_seconds = int(round(eta_seconds))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours > 0:
        return f"{hours}:{minutes:02d}:{seconds:02d}"
    return f"{minutes:02d}:{seconds:02d}"


def _emit_stage_progress(
    *,
    stage: str,
    completed: int,
    total: int,
    failed: int,
    started_at: float,
    last_emitted_at: float,
    interval: float,
    active: int | None = None,
    queued: int | None = None,
    active_sources: Sequence[str] | None = None,
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
    eta_seconds = (remaining / rate) if rate > 0 else None
    percent = (100.0 * completed / total) if total > 0 else 100.0
    suffix = ""
    if active is not None:
        suffix += f", active={active}"
    if queued is not None:
        suffix += f", queued={queued}"
    if active_sources:
        suffix += f", active_sources={','.join(active_sources)}"

    print(
        (
            f"{stage} progress: {completed}/{total} ({percent:.1f}%), "
            f"failed={failed}, rate={rate:.2f}/s, "
            f"eta={_format_eta_seconds(eta_seconds)}{suffix}"
        ),
        file=sys.stderr,
    )
    return now


def _add_common_analysis_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where annotated CIF/JSON outputs and run summary are written.",
    )
    parser.add_argument(
        "--download-dir",
        type=Path,
        default=None,
        help="Directory for downloaded RCSB structures (defaults to <output-dir>/downloads).",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=3.0,
        help="Voxel resolution in Angstroms (default: 3.0).",
    )
    parser.add_argument(
        "--min-voxels",
        type=int,
        default=4,
        help="Minimum voxels cutoff for reporting volumes (default: 4).",
    )
    parser.add_argument(
        "--min-volume",
        type=float,
        default=None,
        help="Optional minimum volume cutoff for reporting.",
    )
    parser.add_argument(
        "--include-hubs",
        action="store_true",
        help=(
            "Include hub volumes in written annotation JSON/CIF outputs. "
            "Default omits hubs."
        ),
    )
    parser.add_argument(
        "--backend",
        choices=sorted(VALID_BACKENDS),
        default=None,
        help="Override backend mode for this run (python|auto|native).",
    )
    parser.add_argument(
        "--surface-connectivity",
        choices=VALID_DIRECT_SURFACE_COMPONENT_CONNECTIVITY_MODES,
        default="custom18",
        help=(
            "Direct-surface smoothing connectivity for classifying "
            "cavity/pocket/pore/hub mouths (default: custom18)."
        ),
    )
    parser.add_argument(
        "--merge-mouth-gap-voxels",
        type=int,
        default=DIRECT_SURFACE_MOUTH_MERGE_GAP_VOXELS,
        help=(
            "Merge very close raw pore mouths into a single human-facing mouth. "
            "Negative values disable the override (default: 1)."
        ),
    )
    parser.add_argument(
        "--assembly-policy",
        choices=VALID_ASSEMBLY_POLICIES,
        default=DEFAULT_ASSEMBLY_POLICY,
        help=(
            "Structure assembly policy for load stage: "
            "biological (default), asymmetric, or auto."
        ),
    )
    parser.add_argument(
        "--keep-non-protein",
        action="store_true",
        help="Keep non-protein residues during structure cleaning.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Parallel worker count for analysis and cluster network steps (default: 1).",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=60.0,
        help="Network timeout in seconds for RCSB requests (default: 60).",
    )
    parser.add_argument(
        "--worker-timeout-seconds",
        type=float,
        default=DEFAULT_ANALYSIS_WORKER_TIMEOUT_SECONDS,
        help=(
            "Max seconds per isolated analysis worker attempt before it is "
            "terminated and treated as failed. Applies to cluster/native "
            f"isolated workers. (default: {int(DEFAULT_ANALYSIS_WORKER_TIMEOUT_SECONDS)}, "
            "set <= 0 to disable)"
        ),
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=2,
        help="Retry count for transient RCSB network errors (default: 2).",
    )
    parser.add_argument(
        "--retry-delay",
        type=float,
        default=1.0,
        help="Base delay in seconds between retries (default: 1.0).",
    )

    write_group = parser.add_mutually_exclusive_group()
    write_group.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing downloaded/output files.",
    )
    write_group.add_argument(
        "--reannotate",
        action="store_true",
        help=(
            "Re-run analysis from existing local input files, overwriting "
            "annotated outputs without refreshing existing downloads."
        ),
    )
    write_group.add_argument(
        "--resume",
        action="store_true",
        help="Skip structures that already have annotated CIF + JSON outputs.",
    )

    checkpoint_group = parser.add_mutually_exclusive_group()
    checkpoint_group.add_argument(
        "--checkpoint",
        type=Path,
        default=None,
        help=(
            "Checkpoint JSON path "
            f"(default: <output-dir>/{DEFAULT_CHECKPOINT_FILENAME})."
        ),
    )
    checkpoint_group.add_argument(
        "--no-checkpoint",
        action="store_true",
        help="Disable checkpoint persistence.",
    )

    parser.add_argument(
        "--progress-jsonl",
        type=Path,
        default=None,
        help="Optional JSONL path for structured per-event progress output.",
    )
    parser.add_argument(
        "--failures-manifest",
        type=Path,
        default=None,
        help=(
            "Optional manifest path for failed structures; output is "
            "replayable with `analyze --manifest`."
        ),
    )
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=30.0,
        help=(
            "Seconds between human-readable progress/ETA updates "
            "(default: 30, set <= 0 to disable)."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Resolve and filter structure inputs, then write run summary without "
            "downloading structure files or running volume analysis."
        ),
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop on first structure-level processing failure.",
    )


def _add_cluster_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--cluster-identity",
        type=int,
        required=True,
        choices=sorted(rcsb.VALID_CLUSTER_IDENTITIES),
        help=(
            "RCSB sequence identity threshold for representative cluster download "
            "(e.g. 30, 50, 90)."
        ),
    )
    parser.add_argument(
        "--max-structures",
        type=int,
        default=None,
        help="Optional cap for selected cluster representatives.",
    )
    parser.add_argument(
        "--num-shards",
        type=int,
        default=None,
        help="Optional number of deterministic cluster shards.",
    )
    parser.add_argument(
        "--shard-index",
        type=int,
        default=None,
        help="0-based shard index to process when --num-shards is set.",
    )
    parser.add_argument(
        "--write-manifest",
        type=Path,
        default=None,
        help=(
            "Optional JSON manifest output path for selected IDs and "
            "filter rejections."
        ),
    )

    cluster_method_group = parser.add_mutually_exclusive_group()
    cluster_method_group.add_argument(
        "--cluster-method",
        action="append",
        default=None,
        metavar="METHOD",
        help=(
            "Allowed method filter(s). Repeatable. "
            "Supported aliases: xray|x-ray, em|cryo-em, nmr, neutron. "
            "Default: xray + em."
        ),
    )
    cluster_method_group.add_argument(
        "--cluster-allow-all-methods",
        action="store_true",
        help="Disable method filtering for cluster selection.",
    )
    parser.add_argument(
        "--cluster-max-resolution",
        type=float,
        default=None,
        help=(
            "Optional max best-resolution (Angstrom) filter. "
            "Entries missing resolution are excluded when set."
        ),
    )
    parser.add_argument(
        "--cluster-max-residues",
        type=int,
        default=DEFAULT_CLUSTER_MAX_RESIDUES,
        help=(
            "Max total deposited polymer residues per entry. "
            f"Entries exceeding this are skipped. "
            f"(default: {DEFAULT_CLUSTER_MAX_RESIDUES})"
        ),
    )
    parser.add_argument(
        "--metadata-cache",
        type=Path,
        default=None,
        help=(
            "Legacy option name for the per-entry metadata-store directory "
            f"(default: <output-dir>/{DEFAULT_METADATA_STORE_DIRNAME})."
        ),
    )
    parser.add_argument(
        "--no-metadata-cache",
        action="store_true",
        help="Disable metadata-store read/write for cluster runs.",
    )


def build_parser() -> argparse.ArgumentParser:
    """
    Build CLI argument parser.
    """
    parser = argparse.ArgumentParser(
        prog="volumizer",
        description=(
            "Analyze one or more structures and write annotated CIF + JSON outputs."
        ),
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {_resolve_cli_version()}",
        help="Show installed volumizer version and exit.",
    )

    subparsers = parser.add_subparsers(dest="command")

    analyze_parser = subparsers.add_parser(
        "analyze",
        help="Analyze a local structure file or one PDB ID.",
    )
    analyze_source_group = analyze_parser.add_mutually_exclusive_group(required=True)
    analyze_source_group.add_argument(
        "--input",
        type=Path,
        help="Path to a local input structure file (.pdb, .cif, .mmtf, etc.).",
    )
    analyze_source_group.add_argument(
        "--pdb-id",
        type=str,
        help="Single PDB ID to download from RCSB and analyze.",
    )
    analyze_source_group.add_argument(
        "--manifest",
        type=Path,
        help=(
            "Path to a manifest JSON containing a `structures` list "
            "(from `cluster --write-manifest` or custom)."
        ),
    )
    analyze_source_group.add_argument(
        "--from-summary",
        type=Path,
        help=(
            "Path to an existing run.summary.json to replay selected "
            "structures."
        ),
    )
    analyze_parser.add_argument(
        "--only",
        choices=("failed", "skipped", "planned", "all"),
        default="failed",
        help=(
            "Used with --from-summary to select which entries to replay "
            "(default: failed)."
        ),
    )
    _add_common_analysis_args(analyze_parser)
    analyze_parser.set_defaults(
        command="analyze",
        cluster_identity=None,
        max_structures=None,
        num_shards=None,
        shard_index=None,
        cluster_method=None,
        cluster_allow_all_methods=False,
        cluster_max_resolution=None,
        cluster_max_residues=None,
        metadata_cache=None,
        no_metadata_cache=False,
        write_manifest=None,
    )

    cluster_parser = subparsers.add_parser(
        "cluster",
        help="Select and analyze cluster representatives from RCSB.",
    )
    _add_cluster_args(cluster_parser)
    _add_common_analysis_args(cluster_parser)
    cluster_parser.set_defaults(
        command="cluster",
        input=None,
        pdb_id=None,
        manifest=None,
        from_summary=None,
        only="failed",
    )

    cache_parser = subparsers.add_parser(
        "cache",
        help="Inspect or mutate metadata cache files.",
    )
    cache_subparsers = cache_parser.add_subparsers(dest="cache_command")
    cache_subparsers.required = True

    cache_inspect_parser = cache_subparsers.add_parser(
        "inspect",
        help="Print metadata-store summary.",
    )
    cache_inspect_parser.add_argument(
        "--metadata-cache",
        type=Path,
        required=True,
        help="Path to metadata-store directory.",
    )

    cache_clear_negative_parser = cache_subparsers.add_parser(
        "clear-negative",
        help="Remove negative metadata-store records in place.",
    )
    cache_clear_negative_parser.add_argument(
        "--metadata-cache",
        type=Path,
        required=True,
        help="Path to metadata-store directory.",
    )

    return parser


def _normalize_argv_for_subcommands(argv_list: list[str]) -> list[str]:
    if len(argv_list) == 0:
        return argv_list

    known_subcommands = {"analyze", "cluster", "cache", "-h", "--help", "-V", "--version"}
    first = argv_list[0]
    if first in known_subcommands:
        return argv_list

    if first.startswith("-"):
        cluster_markers = {
            "--cluster-identity",
            "--max-structures",
            "--num-shards",
            "--shard-index",
            "--cluster-method",
            "--cluster-allow-all-methods",
            "--cluster-max-resolution",
            "--cluster-max-residues",
            "--metadata-cache",
            "--no-metadata-cache",
            "--write-manifest",
        }
        inferred_command = (
            "cluster"
            if any(marker in argv_list for marker in cluster_markers)
            else "analyze"
        )
        return [inferred_command, *argv_list]

    return argv_list


def _sanitize_label(label: str) -> str:
    sanitized = "".join(
        char if (char.isalnum() or char in {"-", "_"}) else "_" for char in label
    )
    sanitized = sanitized.strip("_").lower()
    return sanitized if len(sanitized) > 0 else "structure"


def _output_paths_for_label(output_dir: Path, source_label: str) -> tuple[Path, Path]:
    structure_output_path = output_dir / f"{source_label}.annotated.cif"
    annotation_output_path = output_dir / f"{source_label}.annotation.json"
    return structure_output_path, annotation_output_path


def _status_record_path_for_label(output_dir: Path, source_label: str) -> Path:
    return output_dir / f"{source_label}.status.json"


def _write_status_record(
    output_dir: Path,
    source_label: str,
    payload: dict,
) -> Path:
    status_record_path = _status_record_path_for_label(output_dir, source_label)
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_path = status_record_path.with_name(f".{status_record_path.name}.tmp")
    try:
        with temp_path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
        os.replace(temp_path, status_record_path)
    finally:
        if temp_path.exists():
            temp_path.unlink(missing_ok=True)
    return status_record_path


def _remove_status_record(output_dir: Path, source_label: str) -> None:
    status_record_path = _status_record_path_for_label(output_dir, source_label)
    if not status_record_path.exists():
        return
    if not status_record_path.is_file():
        raise FileExistsError(
            f"Status record path exists but is not a file: {status_record_path}"
        )
    status_record_path.unlink()


def _load_status_record(
    output_dir: Path,
    source_label: str,
) -> tuple[Path, dict | None]:
    status_record_path = _status_record_path_for_label(output_dir, source_label)
    if not status_record_path.is_file():
        return status_record_path, None
    if status_record_path.stat().st_size <= 0:
        return status_record_path, None
    try:
        payload = json.loads(status_record_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return status_record_path, None
    if not isinstance(payload, dict):
        return status_record_path, None
    return status_record_path, payload


def _resolved_path_matches(raw_path: str, expected_path: Path) -> bool:
    try:
        return Path(raw_path).expanduser().resolve() == expected_path.expanduser().resolve()
    except OSError:
        return False


def _get_output_state(output_dir: Path, source_label: str) -> tuple[Path, Path, bool, bool]:
    structure_output_path, annotation_output_path = _output_paths_for_label(
        output_dir,
        source_label,
    )
    structure_exists = structure_output_path.exists()
    annotation_exists = annotation_output_path.exists()
    return (
        structure_output_path,
        annotation_output_path,
        structure_exists,
        annotation_exists,
    )


def _validate_resume_outputs(
    output_dir: Path,
    source_label: str,
) -> dict[str, object]:
    structure_output_path, annotation_output_path, structure_exists, annotation_exists = (
        _get_output_state(
            output_dir,
            source_label,
        )
    )

    if not structure_exists and not annotation_exists:
        return {
            "status": "missing",
            "structure_output_path": structure_output_path,
            "annotation_output_path": annotation_output_path,
            "payload": None,
            "reason": None,
        }

    if structure_exists != annotation_exists:
        present_outputs: list[str] = []
        if structure_exists:
            present_outputs.append(structure_output_path.name)
        if annotation_exists:
            present_outputs.append(annotation_output_path.name)
        return {
            "status": "partial",
            "structure_output_path": structure_output_path,
            "annotation_output_path": annotation_output_path,
            "payload": None,
            "reason": ", ".join(present_outputs),
        }

    for label, path in (
        ("structure output", structure_output_path),
        ("annotation output", annotation_output_path),
    ):
        if not path.is_file():
            return {
                "status": "invalid",
                "structure_output_path": structure_output_path,
                "annotation_output_path": annotation_output_path,
                "payload": None,
                "reason": f"{label} is not a regular file: {path}",
            }

        if path.stat().st_size <= 0:
            return {
                "status": "invalid",
                "structure_output_path": structure_output_path,
                "annotation_output_path": annotation_output_path,
                "payload": None,
                "reason": f"{label} is empty: {path}",
            }

    try:
        payload = json.loads(annotation_output_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        return {
            "status": "invalid",
            "structure_output_path": structure_output_path,
            "annotation_output_path": annotation_output_path,
            "payload": None,
            "reason": f"annotation JSON is invalid: {error}",
        }

    if not isinstance(payload, dict):
        return {
            "status": "invalid",
            "structure_output_path": structure_output_path,
            "annotation_output_path": annotation_output_path,
            "payload": None,
            "reason": "annotation JSON root must be an object",
        }

    required_keys = (
        "source",
        "input_path",
        "output_structure_cif",
        "num_volumes",
        "volumes",
    )
    missing_keys = [key for key in required_keys if key not in payload]
    if missing_keys:
        return {
            "status": "invalid",
            "structure_output_path": structure_output_path,
            "annotation_output_path": annotation_output_path,
            "payload": payload,
            "reason": "annotation JSON missing required keys: "
            + ", ".join(missing_keys),
        }

    if payload["source"] != source_label:
        return {
            "status": "invalid",
            "structure_output_path": structure_output_path,
            "annotation_output_path": annotation_output_path,
            "payload": payload,
            "reason": (
                f"annotation source mismatch: expected {source_label!r}, "
                f"found {payload['source']!r}"
            ),
        }

    raw_structure_output = payload["output_structure_cif"]
    if not isinstance(raw_structure_output, str) or len(raw_structure_output.strip()) == 0:
        return {
            "status": "invalid",
            "structure_output_path": structure_output_path,
            "annotation_output_path": annotation_output_path,
            "payload": payload,
            "reason": "annotation output_structure_cif must be a non-empty string",
        }

    if Path(raw_structure_output).expanduser().resolve() != structure_output_path.resolve():
        return {
            "status": "invalid",
            "structure_output_path": structure_output_path,
            "annotation_output_path": annotation_output_path,
            "payload": payload,
            "reason": (
                "annotation output_structure_cif does not match expected path: "
                f"{raw_structure_output!r}"
            ),
        }

    volumes = payload["volumes"]
    if not isinstance(volumes, list):
        return {
            "status": "invalid",
            "structure_output_path": structure_output_path,
            "annotation_output_path": annotation_output_path,
            "payload": payload,
            "reason": "annotation volumes must be a list",
        }

    num_volumes = payload["num_volumes"]
    if isinstance(num_volumes, bool) or not isinstance(num_volumes, int):
        return {
            "status": "invalid",
            "structure_output_path": structure_output_path,
            "annotation_output_path": annotation_output_path,
            "payload": payload,
            "reason": "annotation num_volumes must be an integer",
        }

    if num_volumes != len(volumes):
        return {
            "status": "invalid",
            "structure_output_path": structure_output_path,
            "annotation_output_path": annotation_output_path,
            "payload": payload,
            "reason": (
                "annotation num_volumes does not match volumes length: "
                f"{num_volumes} != {len(volumes)}"
            ),
        }

    return {
        "status": "valid",
        "structure_output_path": structure_output_path,
        "annotation_output_path": annotation_output_path,
        "payload": payload,
        "reason": None,
    }


def _get_resume_action(
    *,
    source_label: str,
    input_path: Path,
    output_dir: Path,
    assembly_policy: str,
    max_residues: int | None,
) -> dict[str, object]:
    validation = _validate_resume_outputs(
        output_dir,
        source_label,
    )
    status = validation["status"]
    if status == "valid":
        return {
            "action": "skip",
            "reason_kind": "valid",
            "skip_entry": _build_resume_skip_entry(
                source_label=source_label,
                input_path=input_path,
                output_dir=output_dir,
                annotation_payload=validation["payload"],
            ),
            "overwrite_existing_outputs": False,
            "log_message": None,
        }

    if status == "partial":
        return {
            "action": "analyze",
            "reason_kind": "partial",
            "skip_entry": None,
            "overwrite_existing_outputs": False,
            "log_message": (
                f"re-analyzing {source_label}: partial outputs found (--resume): "
                f"{validation['reason']}"
            ),
        }

    if status == "invalid":
        return {
            "action": "analyze",
            "reason_kind": "invalid",
            "skip_entry": None,
            "overwrite_existing_outputs": True,
            "log_message": (
                f"re-analyzing {source_label}: invalid existing outputs (--resume): "
                f"{validation['reason']}"
            ),
        }

    cached_terminal_skip_entry = _get_cached_terminal_skip_entry(
        source_label=source_label,
        input_path=input_path,
        output_dir=output_dir,
        assembly_policy=assembly_policy,
        max_residues=max_residues,
    )
    if cached_terminal_skip_entry is not None:
        return {
            "action": "skip",
            "reason_kind": "terminal_skip",
            "skip_entry": cached_terminal_skip_entry,
            "overwrite_existing_outputs": False,
            "log_message": (
                "cached terminal skip (--resume): "
                f"post_assembly_residue_limit "
                f"{cached_terminal_skip_entry['actual_residues']} > "
                f"{cached_terminal_skip_entry['max_residues']}"
            ),
        }

    return {
        "action": "analyze",
        "reason_kind": "missing",
        "skip_entry": None,
        "overwrite_existing_outputs": False,
        "log_message": None,
    }


def _format_resume_reason_samples(samples: list[str], total: int) -> str:
    if total == 0 or len(samples) == 0:
        return ""

    suffix = ", ..." if total > len(samples) else ""
    return ", ".join(samples) + suffix


def _emit_resume_scan_summary(
    *,
    skipped_valid: int,
    skipped_terminal: int,
    reanalyze_partial: int,
    reanalyze_invalid: int,
    terminal_samples: list[str],
    partial_samples: list[str],
    invalid_samples: list[str],
) -> None:
    if (
        skipped_valid == 0
        and skipped_terminal == 0
        and reanalyze_partial == 0
        and reanalyze_invalid == 0
    ):
        return

    print(
        (
            "resume scan: "
            f"skipped_valid={skipped_valid}, "
            f"skipped_terminal={skipped_terminal}, "
            f"reanalyze_partial={reanalyze_partial}, "
            f"reanalyze_invalid={reanalyze_invalid}"
        ),
        file=sys.stderr,
    )

    terminal_summary = _format_resume_reason_samples(terminal_samples, skipped_terminal)
    if terminal_summary:
        print(
            f"resume terminal examples: {terminal_summary}",
            file=sys.stderr,
        )

    partial_summary = _format_resume_reason_samples(partial_samples, reanalyze_partial)
    if partial_summary:
        print(
            f"resume partial examples: {partial_summary}",
            file=sys.stderr,
        )

    invalid_summary = _format_resume_reason_samples(invalid_samples, reanalyze_invalid)
    if invalid_summary:
        print(
            f"resume invalid examples: {invalid_summary}",
            file=sys.stderr,
        )


def _emit_resume_scan_progress(
    *,
    processed: int,
    total: int,
    skipped_valid: int,
    skipped_terminal: int,
    reanalyze_partial: int,
    reanalyze_invalid: int,
) -> None:
    percent = (100.0 * processed / total) if total > 0 else 100.0
    print(
        (
            f"resume scan progress: {processed}/{total} ({percent:.1f}%), "
            f"skipped_valid={skipped_valid}, "
            f"skipped_terminal={skipped_terminal}, "
            f"reanalyze_partial={reanalyze_partial}, "
            f"reanalyze_invalid={reanalyze_invalid}"
        ),
        file=sys.stderr,
    )


def _prepare_pending_structures(
    *,
    structures: list[tuple[str, Path]],
    args: argparse.Namespace,
    output_dir: Path,
    tracker: _RunTracker,
) -> list[tuple[str, Path, bool]]:
    pending_structures: list[tuple[str, Path, bool]] = []
    skipped_valid = 0
    skipped_terminal = 0
    reanalyze_partial = 0
    reanalyze_invalid = 0
    terminal_samples: list[str] = []
    partial_samples: list[str] = []
    invalid_samples: list[str] = []
    pending_skip_persists = 0
    progress_interval = float(args.progress_interval)
    progress_started_at = time.monotonic()
    last_progress_emit = progress_started_at
    total_structures = len(structures)
    runtime_max_residues = _get_post_assembly_max_residues(args)

    if args.resume and total_structures > 0:
        print(
            f"resume scan: validating existing outputs for {total_structures} structures...",
            file=sys.stderr,
        )

    for index, (source_label, input_path) in enumerate(structures, start=1):
        overwrite_existing_outputs = False

        if args.resume:
            resume_action = _get_resume_action(
                source_label=source_label,
                input_path=input_path,
                output_dir=output_dir,
                assembly_policy=str(args.assembly_policy),
                max_residues=runtime_max_residues,
            )
            reason_kind = str(resume_action["reason_kind"])

            if resume_action["action"] == "skip":
                tracker.mark_skipped(
                    _enrich_structure_entry(
                        dict(resume_action["skip_entry"]),
                        args=args,
                        source_label=source_label,
                    ),
                    persist=False,
                )
                pending_skip_persists += 1
                if reason_kind == "terminal_skip":
                    skipped_terminal += 1
                    log_message = resume_action["log_message"]
                    if log_message is not None and len(terminal_samples) < 3:
                        terminal_samples.append(
                            f"{source_label} ({log_message.removeprefix('cached terminal skip (--resume): ')})"
                        )
                else:
                    skipped_valid += 1
            else:
                overwrite_existing_outputs = bool(
                    resume_action["overwrite_existing_outputs"]
                )
                log_message = resume_action["log_message"]
                if reason_kind == "partial":
                    reanalyze_partial += 1
                    if log_message is not None and len(partial_samples) < 3:
                        partial_samples.append(
                            log_message.removeprefix(
                                f"re-analyzing {source_label}: partial outputs found (--resume): "
                            )
                        )
                        partial_samples[-1] = f"{source_label} ({partial_samples[-1]})"
                elif reason_kind == "invalid":
                    reanalyze_invalid += 1
                    if log_message is not None and len(invalid_samples) < 3:
                        invalid_samples.append(
                            log_message.removeprefix(
                                f"re-analyzing {source_label}: invalid existing outputs (--resume): "
                            )
                        )
                        invalid_samples[-1] = f"{source_label} ({invalid_samples[-1]})"

                pending_structures.append(
                    (source_label, input_path, overwrite_existing_outputs)
                )
        else:
            pending_structures.append(
                (source_label, input_path, overwrite_existing_outputs)
            )

        if args.resume and progress_interval > 0:
            now = time.monotonic()
            if (now - last_progress_emit) >= progress_interval:
                _emit_resume_scan_progress(
                    processed=index,
                    total=total_structures,
                    skipped_valid=skipped_valid,
                    skipped_terminal=skipped_terminal,
                    reanalyze_partial=reanalyze_partial,
                    reanalyze_invalid=reanalyze_invalid,
                )
                if pending_skip_persists > 0:
                    tracker.persist_checkpoint()
                    pending_skip_persists = 0
                last_progress_emit = now

    if pending_skip_persists > 0:
        tracker.persist_checkpoint()

    if args.resume:
        _emit_resume_scan_summary(
            skipped_valid=skipped_valid,
            skipped_terminal=skipped_terminal,
            reanalyze_partial=reanalyze_partial,
            reanalyze_invalid=reanalyze_invalid,
            terminal_samples=terminal_samples,
            partial_samples=partial_samples,
            invalid_samples=invalid_samples,
        )

    return pending_structures


def _cleanup_partial_outputs(
    structure_output_path: Path,
    annotation_output_path: Path,
) -> list[Path]:
    structure_exists = structure_output_path.exists()
    annotation_exists = annotation_output_path.exists()
    if structure_exists == annotation_exists:
        return []

    removed: list[Path] = []
    for path in (structure_output_path, annotation_output_path):
        if not path.exists():
            continue
        if not path.is_file():
            raise FileExistsError(
                f"Partial output path exists but is not a file: {path}"
            )
        path.unlink()
        removed.append(path)
    return removed


def _callable_accepts_keyword(fn, keyword: str) -> bool:
    try:
        signature = inspect.signature(fn)
    except (TypeError, ValueError):
        return True

    if keyword in signature.parameters:
        return True

    return any(
        parameter.kind == inspect.Parameter.VAR_KEYWORD
        for parameter in signature.parameters.values()
    )


def _format_exception_message(error: BaseException) -> str:
    message = str(error).strip()
    error_type = type(error).__name__
    if message:
        return f"{error_type}: {message}"

    repr_message = repr(error).strip()
    if repr_message:
        return repr_message
    return error_type


def _build_error_entry(
    *,
    source_label: str,
    input_path: Path,
    error: BaseException,
) -> dict[str, str]:
    return {
        "source": source_label,
        "input_path": str(input_path),
        "error": _format_exception_message(error),
        "error_type": type(error).__name__,
        "error_repr": repr(error),
        "traceback": "".join(
            traceback.format_exception(type(error), error, error.__traceback__)
        ),
    }


def _build_resume_skip_entry(
    source_label: str,
    input_path: Path,
    output_dir: Path,
    annotation_payload: dict | None = None,
) -> dict:
    structure_output_path, annotation_output_path = _output_paths_for_label(
        output_dir,
        source_label,
    )

    payload = annotation_payload if annotation_payload is not None else {}
    if annotation_payload is None and annotation_output_path.is_file():
        try:
            payload = json.loads(annotation_output_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            payload = {}

    return {
        "source": source_label,
        "pdb_id": _infer_pdb_id_for_result(source_label, input_path),
        "input_path": str(input_path),
        "structure_output": str(structure_output_path),
        "annotation_output": str(annotation_output_path),
        "reason": "resume_existing_outputs",
        "num_volumes": payload.get("num_volumes"),
        "largest_type": payload.get("largest_type"),
        "largest_volume": payload.get("largest_volume"),
    }


def _build_post_assembly_residue_limit_skip_entry(
    *,
    source_label: str,
    input_path: Path,
    output_dir: Path,
    error: PostAssemblyResidueLimitExceeded,
) -> dict:
    structure_output_path, annotation_output_path = _output_paths_for_label(
        output_dir,
        source_label,
    )
    return {
        "source": source_label,
        "pdb_id": _infer_pdb_id_for_result(source_label, input_path),
        "input_path": str(input_path),
        "structure_output": str(structure_output_path),
        "annotation_output": str(annotation_output_path),
        "reason": "post_assembly_residue_limit",
        "actual_residues": error.actual_residues,
        "max_residues": error.max_residues,
        "assembly_policy": error.assembly_policy,
    }


def _build_post_assembly_residue_limit_status_payload(
    *,
    source_label: str,
    input_path: Path,
    output_dir: Path,
    error: PostAssemblyResidueLimitExceeded,
) -> dict:
    skip_entry = _build_post_assembly_residue_limit_skip_entry(
        source_label=source_label,
        input_path=input_path,
        output_dir=output_dir,
        error=error,
    )
    return {
        "status_record_format": STATUS_RECORD_FORMAT_VERSION,
        "kind": "terminal_skip",
        **skip_entry,
    }


def _write_post_assembly_residue_limit_status_record(
    *,
    source_label: str,
    input_path: Path,
    output_dir: Path,
    error: PostAssemblyResidueLimitExceeded,
) -> Path:
    return _write_status_record(
        output_dir,
        source_label,
        _build_post_assembly_residue_limit_status_payload(
            source_label=source_label,
            input_path=input_path,
            output_dir=output_dir,
            error=error,
        ),
    )


def _record_post_assembly_residue_limit_skip(
    *,
    tracker: _RunTracker,
    source_label: str,
    input_path: Path,
    output_dir: Path,
    error: PostAssemblyResidueLimitExceeded,
    extra_entry_fields: dict[str, object] | None = None,
) -> dict:
    _write_post_assembly_residue_limit_status_record(
        source_label=source_label,
        input_path=input_path,
        output_dir=output_dir,
        error=error,
    )
    skip_entry = _build_post_assembly_residue_limit_skip_entry(
        source_label=source_label,
        input_path=input_path,
        output_dir=output_dir,
        error=error,
    )
    if extra_entry_fields:
        for key, value in extra_entry_fields.items():
            if value is None:
                continue
            skip_entry[key] = list(value) if isinstance(value, list) else value
    tracker.mark_skipped(skip_entry)
    return skip_entry


def _get_cached_terminal_skip_entry(
    *,
    source_label: str,
    input_path: Path,
    output_dir: Path,
    assembly_policy: str,
    max_residues: int | None,
) -> dict | None:
    if max_residues is None:
        return None

    structure_output_path, annotation_output_path = _output_paths_for_label(
        output_dir,
        source_label,
    )
    status_record_path, payload = _load_status_record(output_dir, source_label)
    if payload is None:
        return None
    if payload.get("kind") != "terminal_skip":
        return None
    if payload.get("reason") != "post_assembly_residue_limit":
        return None
    if payload.get("source") != source_label:
        return None
    if payload.get("assembly_policy") != assembly_policy:
        return None

    raw_input_path = payload.get("input_path")
    if not isinstance(raw_input_path, str) or not _resolved_path_matches(
        raw_input_path,
        input_path,
    ):
        return None

    raw_structure_output = payload.get("structure_output")
    if not isinstance(raw_structure_output, str) or not _resolved_path_matches(
        raw_structure_output,
        structure_output_path,
    ):
        return None

    raw_annotation_output = payload.get("annotation_output")
    if not isinstance(raw_annotation_output, str) or not _resolved_path_matches(
        raw_annotation_output,
        annotation_output_path,
    ):
        return None

    actual_residues = payload.get("actual_residues")
    if isinstance(actual_residues, bool) or not isinstance(actual_residues, int):
        return None
    if actual_residues <= int(max_residues):
        return None

    return {
        "source": source_label,
        "pdb_id": payload.get("pdb_id") or _infer_pdb_id_for_result(source_label, input_path),
        "input_path": str(input_path),
        "structure_output": str(structure_output_path),
        "annotation_output": str(annotation_output_path),
        "reason": "post_assembly_residue_limit",
        "actual_residues": actual_residues,
        "max_residues": int(max_residues),
        "assembly_policy": assembly_policy,
        "status_record": str(status_record_path),
        "status_kind": "terminal_skip",
    }


def _build_dry_run_plan_entry(
    source_label: str,
    input_path: Path,
    output_dir: Path,
) -> dict:
    structure_output_path, annotation_output_path = _output_paths_for_label(
        output_dir,
        source_label,
    )
    return {
        "source": source_label,
        "pdb_id": _infer_pdb_id_for_result(source_label, input_path),
        "input_path": str(input_path),
        "structure_output": str(structure_output_path),
        "annotation_output": str(annotation_output_path),
    }


def _resolve_cluster_method_filters(args: argparse.Namespace) -> list[str] | None:
    if args.cluster_allow_all_methods:
        return None

    if args.cluster_method is None or len(args.cluster_method) == 0:
        return list(rcsb.DEFAULT_CLUSTER_METHOD_FILTERS)

    normalized_methods = []
    for raw_method in args.cluster_method:
        normalized_methods.append(rcsb.normalize_method_filter_name(raw_method))

    return list(dict.fromkeys(normalized_methods))


def _normalize_cluster_member_pdb_ids(
    raw_value: object,
    *,
    strict: bool = False,
    context: str = "cluster_member_pdb_ids",
) -> list[str] | None:
    if raw_value is None:
        return None

    if not isinstance(raw_value, list):
        if strict:
            raise ValueError(f"{context} must be a list of PDB IDs.")
        return None

    normalized_ids: list[str] = []
    seen: set[str] = set()
    for index, raw_item in enumerate(raw_value, start=1):
        if not isinstance(raw_item, str):
            if strict:
                raise ValueError(f"{context}[{index}] must be a string PDB ID.")
            continue
        try:
            normalized_id = rcsb.normalize_pdb_id(raw_item)
        except ValueError as error:
            if strict:
                raise ValueError(f"{context}[{index}] is invalid: {raw_item!r}") from error
            continue
        if normalized_id in seen:
            continue
        normalized_ids.append(normalized_id)
        seen.add(normalized_id)

    if len(normalized_ids) == 0:
        return None
    return normalized_ids


def _set_structure_metadata(
    args: argparse.Namespace,
    metadata_by_source: dict[str, dict[str, object]],
) -> None:
    setattr(args, "_structure_metadata_by_source", metadata_by_source)


def _get_structure_metadata(
    args: argparse.Namespace,
    source_label: str,
) -> dict[str, object]:
    metadata_by_source = getattr(args, "_structure_metadata_by_source", None)
    if not isinstance(metadata_by_source, dict):
        return {}
    metadata = metadata_by_source.get(source_label)
    if not isinstance(metadata, dict):
        return {}
    return metadata


def _build_structure_metadata_by_source(
    entries: list[dict[str, object]],
) -> dict[str, dict[str, object]]:
    metadata_by_source: dict[str, dict[str, object]] = {}
    for entry in entries:
        source_label = entry.get("source")
        if not isinstance(source_label, str):
            continue

        cluster_member_pdb_ids = _normalize_cluster_member_pdb_ids(
            entry.get("cluster_member_pdb_ids"),
            strict=False,
        )
        if cluster_member_pdb_ids is None:
            continue

        metadata_by_source[source_label] = {
            "cluster_member_pdb_ids": cluster_member_pdb_ids,
        }

    return metadata_by_source


def _enrich_structure_entry(
    entry: dict[str, object],
    *,
    args: argparse.Namespace,
    source_label: str,
) -> dict[str, object]:
    metadata = _get_structure_metadata(args, source_label)
    if len(metadata) == 0:
        return entry

    enriched = dict(entry)
    for key, value in metadata.items():
        if value is None:
            continue
        if isinstance(value, list):
            enriched[key] = list(value)
        else:
            enriched[key] = value
    return enriched


def _save_manifest(manifest_path: Path, payload: dict) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _load_manifest_structures(manifest_path: Path) -> list[dict]:
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Manifest does not exist: {manifest_path}")

    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid manifest JSON: {manifest_path}: {error}") from error

    if not isinstance(payload, dict):
        raise ValueError(f"Manifest root must be a JSON object: {manifest_path}")

    raw_structures = payload.get("structures")
    if not isinstance(raw_structures, list) or len(raw_structures) == 0:
        raise ValueError(
            "Manifest must include a non-empty `structures` list: "
            f"{manifest_path}"
        )

    manifest_structures: list[dict] = []
    for index, raw_entry in enumerate(raw_structures, start=1):
        if not isinstance(raw_entry, dict):
            raise ValueError(
                f"Manifest structure entry #{index} must be an object: {manifest_path}"
            )

        raw_source = raw_entry.get("source")
        raw_pdb_id = raw_entry.get("pdb_id")
        raw_input_path = raw_entry.get("input_path")

        if raw_pdb_id is None and raw_input_path is None:
            raise ValueError(
                "Manifest structure entries require `pdb_id` or `input_path`: "
                f"{manifest_path} entry #{index}"
            )

        source_label: str | None = None
        if raw_source is not None:
            if not isinstance(raw_source, str):
                raise ValueError(
                    f"Manifest `source` must be a string: {manifest_path} entry #{index}"
                )
            source_label = _sanitize_label(raw_source)

        entry: dict[str, str | Path | list[str]] = {}
        if raw_pdb_id is not None:
            if not isinstance(raw_pdb_id, str):
                raise ValueError(
                    f"Manifest `pdb_id` must be a string: {manifest_path} entry #{index}"
                )
            normalized_id = rcsb.normalize_pdb_id(raw_pdb_id)
            entry["pdb_id"] = normalized_id
            if source_label is None:
                source_label = _sanitize_label(normalized_id)

        if raw_input_path is not None:
            if not isinstance(raw_input_path, str):
                raise ValueError(
                    f"Manifest `input_path` must be a string: {manifest_path} entry #{index}"
                )
            input_path = Path(raw_input_path)
            if not input_path.is_absolute():
                input_path = manifest_path.parent / input_path
            entry["input_path"] = input_path
            if source_label is None:
                source_label = _sanitize_label(input_path.stem)

        cluster_member_pdb_ids = _normalize_cluster_member_pdb_ids(
            raw_entry.get("cluster_member_pdb_ids"),
            strict=True,
            context=(
                "Manifest `cluster_member_pdb_ids` "
                f"{manifest_path} entry #{index}"
            ),
        )
        if cluster_member_pdb_ids is not None:
            entry["cluster_member_pdb_ids"] = cluster_member_pdb_ids

        entry["source"] = source_label if source_label is not None else "structure"
        manifest_structures.append(entry)

    return manifest_structures


def _resolve_manifest_structures(
    manifest_structures: list[dict],
    args: argparse.Namespace,
    download_dir: Path,
    output_dir: Path,
) -> list[tuple[str, Path]]:
    structures: list[tuple[str, Path]] = []
    for entry in manifest_structures:
        source_label = str(entry["source"])
        input_path = entry.get("input_path")
        if isinstance(input_path, Path):
            if not input_path.is_file():
                raise FileNotFoundError(
                    f"Manifest input structure does not exist: {input_path}"
                )
            structures.append((source_label, input_path))
            continue

        normalized_id = str(entry["pdb_id"])
        if args.dry_run:
            structures.append((source_label, download_dir / f"{normalized_id}.cif"))
            continue

        structure_path = rcsb.download_structure_cif(
            normalized_id,
            output_dir=download_dir,
            overwrite=args.overwrite,
            timeout=args.timeout,
            retries=args.retries,
            retry_delay=args.retry_delay,
        )
        structures.append((source_label, structure_path))

    return structures


def _resolve_summary_input_path(raw_input_path: str, summary_path: Path) -> Path:
    input_path = Path(raw_input_path)
    if input_path.is_absolute():
        return input_path

    if input_path.is_file():
        return input_path

    return summary_path.parent / input_path


def _load_structure_entries_from_summary(
    summary_path: Path,
    only: str,
) -> list[dict[str, object]]:
    if not summary_path.is_file():
        raise FileNotFoundError(f"Summary file does not exist: {summary_path}")

    try:
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid summary JSON: {summary_path}: {error}") from error

    if not isinstance(payload, dict):
        raise ValueError(f"Summary root must be a JSON object: {summary_path}")

    section_by_only = {
        "failed": "errors",
        "skipped": "skipped",
        "planned": "planned",
    }
    if only == "all":
        selected_sections = ("errors", "skipped", "planned", "results")
    else:
        selected_sections = (section_by_only[only],)

    structures: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()

    for section in selected_sections:
        section_entries = payload.get(section, [])
        if not isinstance(section_entries, list):
            raise ValueError(
                f"Summary section `{section}` must be a list: {summary_path}"
            )

        for entry in section_entries:
            if not isinstance(entry, dict):
                continue

            raw_source = entry.get("source")
            raw_input_path = entry.get("input_path")
            if not isinstance(raw_source, str) or not isinstance(raw_input_path, str):
                continue

            source_label = _sanitize_label(raw_source)
            input_path = _resolve_summary_input_path(raw_input_path, summary_path)
            key = (source_label, str(input_path))
            if key in seen:
                continue
            seen.add(key)
            structure_entry: dict[str, object] = {
                "source": source_label,
                "input_path": input_path,
            }
            cluster_member_pdb_ids = _normalize_cluster_member_pdb_ids(
                entry.get("cluster_member_pdb_ids"),
                strict=False,
            )
            if cluster_member_pdb_ids is not None:
                structure_entry["cluster_member_pdb_ids"] = cluster_member_pdb_ids
            structures.append(structure_entry)

    if len(structures) == 0:
        raise RuntimeError(
            f"No structures selected from summary for --only {only}: {summary_path}"
        )

    return structures


def _load_structures_from_summary(
    summary_path: Path,
    only: str,
) -> list[tuple[str, Path]]:
    structure_entries = _load_structure_entries_from_summary(summary_path, only)
    return [
        (str(entry["source"]), Path(entry["input_path"]))
        for entry in structure_entries
    ]


def _write_cluster_manifest(
    manifest_path: Path,
    args: argparse.Namespace,
    selected_ids: list[str],
    representative_to_member_ids: dict[str, list[str]],
    cluster_method_filters: list[str] | None,
    rejection_entries: list[dict],
    examined_count: int,
    filtered_by_method: int,
    filtered_by_resolution: int,
    filtered_missing_resolution: int,
    filtered_by_residue_count: int,
    metadata_errors: int,
    metadata_hits: int,
    negative_metadata_hits: int,
    metadata_misses: int,
    metadata_updates: int,
    negative_metadata_updates: int,
) -> None:
    payload = {
        "manifest_format": MANIFEST_FORMAT_VERSION,
        "created_at": _utc_timestamp(),
        "source_command": "cluster",
        "selection": {
            "cluster_identity": args.cluster_identity,
            "max_structures": args.max_structures,
            "num_shards": args.num_shards,
            "shard_index": args.shard_index,
            "cluster_method_filters": cluster_method_filters,
            "cluster_allow_all_methods": args.cluster_allow_all_methods,
            "cluster_max_resolution": args.cluster_max_resolution,
            "cluster_max_residues": args.cluster_max_residues,
        },
        "summary": {
            "selected": len(selected_ids),
            "examined": examined_count,
            "filtered_by_method": filtered_by_method,
            "filtered_by_resolution": filtered_by_resolution,
            "missing_resolution": filtered_missing_resolution,
            "filtered_by_residue_count": filtered_by_residue_count,
            "metadata_errors": metadata_errors,
            "metadata_hits": metadata_hits,
            "negative_metadata_hits": negative_metadata_hits,
            "metadata_misses": metadata_misses,
            "metadata_updates": metadata_updates,
            "negative_metadata_updates": negative_metadata_updates,
        },
        "structures": [],
        "rejections": rejection_entries,
    }
    for representative_id in selected_ids:
        structure_entry: dict[str, object] = {
            "source": _sanitize_label(representative_id),
            "pdb_id": representative_id,
        }
        cluster_member_pdb_ids = representative_to_member_ids.get(representative_id)
        if cluster_member_pdb_ids is not None:
            structure_entry["cluster_member_pdb_ids"] = list(cluster_member_pdb_ids)
        payload["structures"].append(structure_entry)
    _save_manifest(manifest_path, payload)
    print(f"wrote manifest: {manifest_path}", file=sys.stderr)


def _write_failures_manifest(
    manifest_path: Path,
    command: str,
    errors: list[dict],
    summary_path: Path | None = None,
) -> None:
    structures: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()

    for entry in errors:
        raw_source = entry.get("source")
        raw_input_path = entry.get("input_path")
        if not isinstance(raw_source, str) or not isinstance(raw_input_path, str):
            continue

        source_label = _sanitize_label(raw_source)
        normalized_input_path = str(Path(raw_input_path).resolve(strict=False))
        dedupe_key = (source_label, normalized_input_path)
        if dedupe_key in seen:
            continue
        seen.add(dedupe_key)

        structures.append(
            {
                "source": source_label,
                "input_path": normalized_input_path,
            }
        )

    payload = {
        "manifest_format": MANIFEST_FORMAT_VERSION,
        "created_at": _utc_timestamp(),
        "source_command": command,
        "source_summary": str(summary_path) if summary_path is not None else None,
        "num_failed": len(structures),
        "structures": structures,
    }
    _save_manifest(manifest_path, payload)
    print(
        f"wrote failures manifest: {manifest_path} (failures={len(structures)})",
        file=sys.stderr,
    )


def _resolve_metadata_cache_path(
    args: argparse.Namespace,
    output_dir: Path,
) -> Path | None:
    if args.no_metadata_cache:
        return None

    if args.metadata_cache is not None:
        return Path(args.metadata_cache)

    return output_dir / DEFAULT_METADATA_STORE_DIRNAME


def _resolve_checkpoint_path(
    args: argparse.Namespace,
    output_dir: Path,
) -> Path | None:
    if args.no_checkpoint:
        return None

    if args.checkpoint is not None:
        return Path(args.checkpoint)

    return output_dir / DEFAULT_CHECKPOINT_FILENAME


def _resolve_progress_jsonl_path(args: argparse.Namespace) -> Path | None:
    if args.progress_jsonl is None:
        return None
    return Path(args.progress_jsonl)


def _apply_cluster_shard(
    representative_ids: list[str],
    num_shards: int | None,
    shard_index: int | None,
) -> list[str]:
    if num_shards is None or shard_index is None:
        return representative_ids

    return [
        representative_id
        for index, representative_id in enumerate(representative_ids)
        if (index % num_shards) == shard_index
    ]


def _make_checkpoint_signature(
    args: argparse.Namespace,
    output_dir: Path,
    download_dir: Path,
    cluster_method_filters: list[str] | None,
    metadata_cache_path: Path | None,
) -> dict:
    return {
        "command": args.command,
        "input": str(args.input) if args.input is not None else None,
        "pdb_id": args.pdb_id,
        "manifest": str(args.manifest) if args.manifest is not None else None,
        "from_summary": (
            str(args.from_summary) if args.from_summary is not None else None
        ),
        "from_summary_only": args.only if args.from_summary is not None else None,
        "cluster_identity": args.cluster_identity,
        "max_structures": args.max_structures,
        "num_shards": args.num_shards,
        "shard_index": args.shard_index,
        "cluster_method_filters": cluster_method_filters,
        "cluster_allow_all_methods": args.cluster_allow_all_methods,
        "cluster_max_resolution": args.cluster_max_resolution,
        "cluster_max_residues": args.cluster_max_residues,
        "metadata_cache": str(metadata_cache_path) if metadata_cache_path else None,
        "write_manifest": str(args.write_manifest) if args.write_manifest is not None else None,
        "failures_manifest": (
            str(args.failures_manifest) if args.failures_manifest is not None else None
        ),
        "output_dir": str(output_dir),
        "download_dir": str(download_dir),
        "resolution": args.resolution,
        "min_voxels": args.min_voxels,
        "min_volume": args.min_volume,
        "surface_connectivity": args.surface_connectivity,
        "merge_mouth_gap_voxels": args.merge_mouth_gap_voxels,
        "backend": args.backend,
        "assembly_policy": args.assembly_policy,
        "keep_non_protein": args.keep_non_protein,
        "worker_timeout_seconds": args.worker_timeout_seconds,
        "progress_interval": args.progress_interval,
        "dry_run": args.dry_run,
    }


def _metadata_record_path(metadata_store_dir: Path | None, pdb_id: str) -> Path | None:
    if metadata_store_dir is None:
        return None
    normalized_id = rcsb.normalize_pdb_id(pdb_id)
    return metadata_store_dir / f"{normalized_id}.json"


def _parse_metadata_record_payload(payload: object) -> tuple[dict | None, dict | None]:
    if not isinstance(payload, dict):
        return None, None

    record_kind = payload.get("kind")
    if record_kind == "entry":
        entry_metadata = payload.get("entry_metadata")
        if isinstance(entry_metadata, dict):
            return entry_metadata, None
        return None, None

    if record_kind == "negative":
        negative_entry = payload.get("negative_entry")
        if isinstance(negative_entry, dict):
            return None, negative_entry
        return None, None

    if payload.get("reason") == "permanent_metadata_error":
        return None, payload

    return payload, None


def _load_metadata_record(
    metadata_store_dir: Path | None,
    pdb_id: str,
) -> tuple[dict | None, dict | None]:
    record_path = _metadata_record_path(metadata_store_dir, pdb_id)
    if record_path is None or not record_path.is_file():
        return None, None

    try:
        payload = json.loads(record_path.read_text(encoding="utf-8"))
    except UnicodeDecodeError:
        try:
            payload = json.loads(record_path.read_text(encoding="utf-8", errors="replace"))
        except (json.JSONDecodeError, OSError):
            return None, None
    except (json.JSONDecodeError, OSError):
        return None, None

    return _parse_metadata_record_payload(payload)


def _write_metadata_record(
    metadata_store_dir: Path | None,
    pdb_id: str,
    *,
    entry_metadata: dict | None = None,
    negative_entry: dict | None = None,
) -> None:
    if metadata_store_dir is None:
        return
    if (entry_metadata is None) == (negative_entry is None):
        raise ValueError("Provide exactly one of entry_metadata or negative_entry.")

    record_path = _metadata_record_path(metadata_store_dir, pdb_id)
    if record_path is None:
        return

    metadata_store_dir.mkdir(parents=True, exist_ok=True)
    normalized_id = rcsb.normalize_pdb_id(pdb_id)
    payload = {
        "metadata_record_format": METADATA_RECORD_FORMAT_VERSION,
        "pdb_id": normalized_id,
    }
    if entry_metadata is not None:
        payload["kind"] = "entry"
        payload["entry_metadata"] = entry_metadata
    else:
        payload["kind"] = "negative"
        payload["negative_entry"] = negative_entry

    temp_path = record_path.with_name(f".{record_path.name}.tmp")
    try:
        with temp_path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, separators=(",", ":"))
            handle.write("\n")
        os.replace(temp_path, record_path)
    finally:
        if temp_path.exists():
            temp_path.unlink(missing_ok=True)


def _load_metadata_cache(cache_path: Path | None) -> tuple[dict[str, dict], dict[str, dict]]:
    entries: dict[str, dict] = {}
    negative_entries: dict[str, dict] = {}
    if cache_path is None or not cache_path.is_dir():
        return entries, negative_entries

    for record_path in sorted(cache_path.glob("*.json")):
        try:
            normalized_id = rcsb.normalize_pdb_id(record_path.stem)
        except ValueError:
            continue
        entry_metadata, negative_entry = _load_metadata_record(cache_path, normalized_id)
        if entry_metadata is not None:
            entries[normalized_id] = entry_metadata
        elif negative_entry is not None:
            negative_entries[normalized_id] = negative_entry

    return entries, negative_entries


def _clear_negative_metadata_records(cache_path: Path | None) -> int:
    if cache_path is None or not cache_path.is_dir():
        return 0

    removed = 0
    for record_path in sorted(cache_path.glob("*.json")):
        try:
            payload = json.loads(record_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError, UnicodeDecodeError):
            continue
        _, negative_entry = _parse_metadata_record_payload(payload)
        if negative_entry is None:
            continue
        record_path.unlink(missing_ok=True)
        removed += 1

    return removed


class _RunTracker:
    """
    Tracks in-memory run state and persists checkpoint/progress events.
    """

    def __init__(
        self,
        checkpoint_path: Path | None,
        progress_jsonl_path: Path | None,
        signature: dict,
        resume: bool,
    ):
        self.checkpoint_path = checkpoint_path
        self.progress_jsonl_path = progress_jsonl_path
        self.signature = signature
        self.resume = bool(resume)

        self._entries: dict[str, dict[str, dict]] = {
            "results": {},
            "errors": {},
            "skipped": {},
            "planned": {},
        }
        self._order: dict[str, list[str]] = {
            "results": [],
            "errors": [],
            "skipped": [],
            "planned": [],
        }

        self._init_progress_stream()

    def _init_progress_stream(self) -> None:
        if self.progress_jsonl_path is None:
            return

        self.progress_jsonl_path.parent.mkdir(parents=True, exist_ok=True)
        self.progress_jsonl_path.write_text("", encoding="utf-8")

    def _remove_source_from_kind(self, kind: str, source: str) -> None:
        if source in self._entries[kind]:
            del self._entries[kind][source]
        if source in self._order[kind]:
            self._order[kind].remove(source)

    def _remove_source_everywhere(self, source: str) -> None:
        for kind in ("results", "errors", "skipped", "planned"):
            self._remove_source_from_kind(kind, source)

    def _set_entry(self, kind: str, entry: dict) -> None:
        source = entry["source"]
        self._remove_source_everywhere(source)
        self._entries[kind][source] = entry
        self._order[kind].append(source)

    def _set_non_terminal_entry(self, kind: str, entry: dict) -> None:
        source = entry["source"]
        for terminal_kind in ("results", "errors", "skipped"):
            self._remove_source_from_kind(terminal_kind, source)
        self._remove_source_from_kind(kind, source)
        self._entries[kind][source] = entry
        self._order[kind].append(source)

    @property
    def results(self) -> list[dict]:
        return [self._entries["results"][s] for s in self._order["results"]]

    @property
    def errors(self) -> list[dict]:
        return [self._entries["errors"][s] for s in self._order["errors"]]

    @property
    def skipped(self) -> list[dict]:
        return [self._entries["skipped"][s] for s in self._order["skipped"]]

    @property
    def planned(self) -> list[dict]:
        return [self._entries["planned"][s] for s in self._order["planned"]]

    def has_result(self, source: str) -> bool:
        return source in self._entries["results"]

    def emit_event(self, event: str, **payload: object) -> None:
        if self.progress_jsonl_path is None:
            return

        event_record = {
            "ts": _utc_timestamp(),
            "event": event,
        }
        event_record.update(payload)
        with self.progress_jsonl_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(event_record) + "\n")

    def persist_checkpoint(self) -> None:
        if self.checkpoint_path is None:
            return

        self.checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "checkpoint_format": 1,
            "updated_at": _utc_timestamp(),
            "signature": self.signature,
            "results": self.results,
            "errors": self.errors,
            "skipped": self.skipped,
            "planned": self.planned,
        }
        self.checkpoint_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    def mark_run_start(self, total_structures: int, dry_run: bool) -> None:
        self.emit_event(
            "run_started",
            total_structures=total_structures,
            dry_run=dry_run,
            resume=self.resume,
        )
        self.persist_checkpoint()

    def mark_run_complete(self, exit_code: int) -> None:
        self.emit_event(
            "run_completed",
            exit_code=exit_code,
            num_processed=len(self.results),
            num_failed=len(self.errors),
            num_skipped=len(self.skipped),
            num_planned=len(self.planned),
        )
        self.persist_checkpoint()

    def mark_planned(self, entry: dict) -> None:
        self._set_non_terminal_entry("planned", entry)
        self.emit_event(
            "structure_planned",
            source=entry["source"],
            input_path=entry["input_path"],
        )
        self.persist_checkpoint()

    def mark_skipped(self, entry: dict, persist: bool = True) -> None:
        self._set_entry("skipped", entry)
        self.emit_event(
            "structure_skipped",
            source=entry["source"],
            input_path=entry["input_path"],
            reason=entry.get("reason"),
        )
        if persist:
            self.persist_checkpoint()

    def mark_started(self, source: str, input_path: Path) -> None:
        self.emit_event(
            "structure_started",
            source=source,
            input_path=str(input_path),
        )

    def mark_result(self, entry: dict, persist: bool = True) -> None:
        self._set_entry("results", entry)
        self.emit_event(
            "structure_succeeded",
            source=entry["source"],
            input_path=entry["input_path"],
            num_volumes=entry.get("num_volumes"),
        )
        if persist:
            self.persist_checkpoint()

    def mark_error(self, entry: dict, persist: bool = True) -> None:
        self._set_entry("errors", entry)
        self.emit_event(
            "structure_failed",
            source=entry["source"],
            input_path=entry["input_path"],
            error=entry["error"],
        )
        if persist:
            self.persist_checkpoint()


def _extract_status_code_from_error_message(error: Exception) -> int | None:
    match = re.search(r"\b([1-5][0-9]{2})\b", str(error))
    if match is None:
        return None
    try:
        return int(match.group(1))
    except ValueError:
        return None


def _is_permanent_metadata_error(error: Exception) -> tuple[bool, int | None]:
    if isinstance(error, rcsb.RCSBFetchError):
        status_code = error.status_code
        is_permanent = bool(error.permanent)
        if status_code in rcsb.TERMINAL_HTTP_STATUS_CODES:
            is_permanent = True
        return is_permanent, status_code

    status_code = _extract_status_code_from_error_message(error)
    if status_code in rcsb.TERMINAL_HTTP_STATUS_CODES:
        return True, status_code

    return False, status_code


def _build_negative_cache_entry(
    error: Exception,
    status_code: int | None,
) -> dict:
    return {
        "reason": "permanent_metadata_error",
        "status_code": status_code,
        "error": _format_exception_message(error),
    }


def _fetch_metadata_for_ids(
    representative_ids: list[str],
    args: argparse.Namespace,
    suppress_progress: bool = False,
) -> tuple[dict[str, dict], int, dict[str, dict], list[dict]]:
    fetched: dict[str, dict] = {}
    errors = 0
    negative_updates: dict[str, dict] = {}
    error_entries: list[dict] = []

    total_ids = len(representative_ids)
    if total_ids == 0:
        return fetched, errors, negative_updates, error_entries

    progress_interval = float(args.progress_interval) if not suppress_progress else 0
    processed = 0
    progress_started_at = time.monotonic()
    last_progress_emit = progress_started_at

    if args.jobs <= 1:
        for representative_id in representative_ids:
            try:
                fetched[representative_id] = rcsb.fetch_entry_metadata(
                    representative_id,
                    timeout=args.timeout,
                    retries=args.retries,
                    retry_delay=args.retry_delay,
                )
            except Exception as error:
                errors += 1
                is_permanent, status_code = _is_permanent_metadata_error(error)
                if is_permanent:
                    negative_updates[representative_id] = _build_negative_cache_entry(
                        error,
                        status_code,
                    )
                error_entries.append(
                    {
                        "pdb_id": representative_id,
                        "reason": "metadata_error",
                        "status_code": status_code,
                        "permanent": is_permanent,
                        "error": _format_exception_message(error),
                        "error_type": type(error).__name__,
                        "error_repr": repr(error),
                        "traceback": "".join(
                            traceback.format_exception(
                                type(error),
                                error,
                                error.__traceback__,
                            )
                        ),
                    }
                )
                print(
                    f"metadata error for {representative_id}: "
                    f"{_format_exception_message(error)}",
                    file=sys.stderr,
                )
                processed += 1
                last_progress_emit = _emit_stage_progress(
                    stage="metadata",
                    completed=processed,
                    total=total_ids,
                    failed=errors,
                    started_at=progress_started_at,
                    last_emitted_at=last_progress_emit,
                    interval=progress_interval,
                )
                if args.fail_fast:
                    raise
                continue

            processed += 1
            last_progress_emit = _emit_stage_progress(
                stage="metadata",
                completed=processed,
                total=total_ids,
                failed=errors,
                started_at=progress_started_at,
                last_emitted_at=last_progress_emit,
                interval=progress_interval,
            )

        _emit_stage_progress(
            stage="metadata",
            completed=processed,
            total=total_ids,
            failed=errors,
            started_at=progress_started_at,
            last_emitted_at=last_progress_emit,
            interval=progress_interval,
            force=not suppress_progress,
        )
        return fetched, errors, negative_updates, error_entries

    max_workers = min(args.jobs, total_ids)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_id = {
            executor.submit(
                rcsb.fetch_entry_metadata,
                representative_id,
                timeout=args.timeout,
                retries=args.retries,
                retry_delay=args.retry_delay,
            ): representative_id
            for representative_id in representative_ids
        }

        for future in as_completed(future_to_id):
            representative_id = future_to_id[future]
            try:
                fetched[representative_id] = future.result()
            except Exception as error:
                errors += 1
                is_permanent, status_code = _is_permanent_metadata_error(error)
                if is_permanent:
                    negative_updates[representative_id] = _build_negative_cache_entry(
                        error,
                        status_code,
                    )
                error_entries.append(
                    {
                        "pdb_id": representative_id,
                        "reason": "metadata_error",
                        "status_code": status_code,
                        "permanent": is_permanent,
                        "error": _format_exception_message(error),
                        "error_type": type(error).__name__,
                        "error_repr": repr(error),
                        "traceback": "".join(
                            traceback.format_exception(
                                type(error),
                                error,
                                error.__traceback__,
                            )
                        ),
                    }
                )
                print(
                    f"metadata error for {representative_id}: "
                    f"{_format_exception_message(error)}",
                    file=sys.stderr,
                )
                processed += 1
                last_progress_emit = _emit_stage_progress(
                    stage="metadata",
                    completed=processed,
                    total=total_ids,
                    failed=errors,
                    started_at=progress_started_at,
                    last_emitted_at=last_progress_emit,
                    interval=progress_interval,
                )
                if args.fail_fast:
                    raise
                continue

            processed += 1
            last_progress_emit = _emit_stage_progress(
                stage="metadata",
                completed=processed,
                total=total_ids,
                failed=errors,
                started_at=progress_started_at,
                last_emitted_at=last_progress_emit,
                interval=progress_interval,
            )

    _emit_stage_progress(
        stage="metadata",
        completed=processed,
        total=total_ids,
        failed=errors,
        started_at=progress_started_at,
        last_emitted_at=last_progress_emit,
        interval=progress_interval,
        force=not suppress_progress,
    )
    return fetched, errors, negative_updates, error_entries


def _download_cluster_structures(
    selected_ids: list[str],
    args: argparse.Namespace,
    download_dir: Path,
    output_dir: Path,
) -> dict[str, Path]:
    del output_dir
    ids_to_download = []
    downloaded_paths: dict[str, Path] = {}
    for representative_id in selected_ids:
        existing_path = download_dir / f"{representative_id}.cif"
        if existing_path.is_file() and not args.overwrite:
            downloaded_paths[representative_id] = existing_path
        else:
            ids_to_download.append(representative_id)

    if len(ids_to_download) == 0:
        return downloaded_paths

    progress_interval = float(args.progress_interval)
    total_downloads = len(ids_to_download)
    completed = 0
    progress_started_at = time.monotonic()
    last_progress_emit = progress_started_at

    if args.jobs <= 1:
        for index, representative_id in enumerate(ids_to_download, start=1):
            print(
                f"[{index}/{total_downloads}] downloading {representative_id}...",
                file=sys.stderr,
            )
            downloaded_paths[representative_id] = rcsb.download_structure_cif(
                representative_id,
                output_dir=download_dir,
                overwrite=args.overwrite,
                timeout=args.timeout,
                retries=args.retries,
                retry_delay=args.retry_delay,
            )
            completed += 1
            last_progress_emit = _emit_stage_progress(
                stage="download",
                completed=completed,
                total=total_downloads,
                failed=0,
                started_at=progress_started_at,
                last_emitted_at=last_progress_emit,
                interval=progress_interval,
            )

        _emit_stage_progress(
            stage="download",
            completed=completed,
            total=total_downloads,
            failed=0,
            started_at=progress_started_at,
            last_emitted_at=last_progress_emit,
            interval=progress_interval,
            force=True,
        )
        return downloaded_paths

    print(
        f"downloading {total_downloads} structures with {args.jobs} workers...",
        file=sys.stderr,
    )
    with ThreadPoolExecutor(max_workers=min(args.jobs, total_downloads)) as executor:
        future_to_id = {
            executor.submit(
                rcsb.download_structure_cif,
                representative_id,
                output_dir=download_dir,
                overwrite=args.overwrite,
                timeout=args.timeout,
                retries=args.retries,
                retry_delay=args.retry_delay,
            ): representative_id
            for representative_id in ids_to_download
        }

        for future in as_completed(future_to_id):
            representative_id = future_to_id[future]
            downloaded_paths[representative_id] = future.result()
            completed += 1
            last_progress_emit = _emit_stage_progress(
                stage="download",
                completed=completed,
                total=total_downloads,
                failed=0,
                started_at=progress_started_at,
                last_emitted_at=last_progress_emit,
                interval=progress_interval,
            )

    _emit_stage_progress(
        stage="download",
        completed=completed,
        total=total_downloads,
        failed=0,
        started_at=progress_started_at,
        last_emitted_at=last_progress_emit,
        interval=progress_interval,
        force=True,
    )
    return downloaded_paths


def resolve_input_structures(
    args: argparse.Namespace,
    download_dir: Path,
    output_dir: Path,
    metadata_cache_path: Path | None | object = _UNSET,
) -> list[tuple[str, Path]]:
    """
    Resolve CLI source arguments into labeled local structure paths.
    """
    _set_structure_metadata(args, {})

    if args.input is not None:
        input_path = Path(args.input)
        if not input_path.is_file():
            raise FileNotFoundError(f"Input structure does not exist: {input_path}")
        return [(_sanitize_label(input_path.stem), input_path)]

    if args.from_summary is not None:
        summary_path = Path(args.from_summary)
        summary_entries = _load_structure_entries_from_summary(summary_path, args.only)
        _set_structure_metadata(args, _build_structure_metadata_by_source(summary_entries))
        summary_structures = [
            (str(entry["source"]), Path(entry["input_path"]))
            for entry in summary_entries
        ]
        for _, input_path in summary_structures:
            if not input_path.is_file():
                raise FileNotFoundError(
                    "Summary replay input structure does not exist: "
                    f"{input_path}"
                )
        return summary_structures

    if args.manifest is not None:
        manifest_path = Path(args.manifest)
        manifest_structures = _load_manifest_structures(manifest_path)
        _set_structure_metadata(args, _build_structure_metadata_by_source(manifest_structures))
        return _resolve_manifest_structures(
            manifest_structures,
            args,
            download_dir,
            output_dir,
        )

    if args.pdb_id is not None:
        normalized_id = rcsb.normalize_pdb_id(args.pdb_id)
        source_label = _sanitize_label(normalized_id)

        if args.dry_run:
            return [(source_label, download_dir / f"{normalized_id}.cif")]

        structure_path = rcsb.download_structure_cif(
            normalized_id,
            output_dir=download_dir,
            overwrite=args.overwrite,
            timeout=args.timeout,
            retries=args.retries,
            retry_delay=args.retry_delay,
        )
        return [(source_label, structure_path)]

    if args.cluster_identity is not None:
        cluster_method_filters = _resolve_cluster_method_filters(args)
        if metadata_cache_path is _UNSET:
            resolved_metadata_cache_path = _resolve_metadata_cache_path(args, output_dir)
        else:
            resolved_metadata_cache_path = metadata_cache_path
        if resolved_metadata_cache_path is None:
            print("metadata store: disabled", file=sys.stderr)
        else:
            print(f"metadata store: {resolved_metadata_cache_path}", file=sys.stderr)

        print(
            f"fetching cluster representatives for identity={args.cluster_identity}...",
            file=sys.stderr,
        )
        representative_to_member_ids = rcsb.fetch_cluster_representative_member_entry_ids(
            identity=args.cluster_identity,
            max_structures=None,
            timeout=args.timeout,
            retries=args.retries,
            retry_delay=args.retry_delay,
        )
        representative_ids = list(representative_to_member_ids.keys())
        total_representative_ids = len(representative_ids)
        if total_representative_ids == 0:
            raise RuntimeError(
                f"No representative PDB IDs found for cluster identity {args.cluster_identity}."
            )

        representative_ids = _apply_cluster_shard(
            representative_ids,
            args.num_shards,
            args.shard_index,
        )
        if args.num_shards is not None:
            print(
                (
                    "cluster shard selection: "
                    f"shard_index={args.shard_index}, "
                    f"num_shards={args.num_shards}, "
                    f"shard_representatives={len(representative_ids)}, "
                    f"total_representatives={total_representative_ids}"
                ),
                file=sys.stderr,
            )
        if len(representative_ids) == 0:
            raise RuntimeError(
                "No representative PDB IDs in selected shard. "
                f"identity={args.cluster_identity}, "
                f"shard_index={args.shard_index}, "
                f"num_shards={args.num_shards}, "
                f"total_representatives={total_representative_ids}."
            )

        target_selection = len(representative_ids)
        if args.max_structures is not None and args.max_structures > 0:
            target_selection = args.max_structures

        selected_ids: list[str] = []
        examined_count = 0
        filtered_by_method = 0
        filtered_by_resolution = 0
        filtered_missing_resolution = 0
        filtered_by_residue_count = 0
        metadata_errors = 0
        metadata_hits = 0
        negative_metadata_hits = 0
        metadata_misses = 0
        metadata_updates = 0
        negative_metadata_updates = 0
        manifest_rejections: list[dict] = []

        batch_size = max(1, int(args.jobs))
        total_representatives = len(representative_ids)
        progress_interval = float(args.progress_interval)
        meta_progress_started_at = time.monotonic()
        meta_last_progress_emit = meta_progress_started_at
        meta_processed = 0

        for batch_start in range(0, total_representatives, batch_size):
            if len(selected_ids) >= target_selection:
                break

            batch_ids = representative_ids[batch_start : batch_start + batch_size]
            batch_metadata: dict[str, dict] = {}
            missing_ids: list[str] = []

            for representative_id in batch_ids:
                cached_metadata, cached_negative = _load_metadata_record(
                    resolved_metadata_cache_path,
                    representative_id,
                )
                if isinstance(cached_metadata, dict):
                    batch_metadata[representative_id] = cached_metadata
                    metadata_hits += 1
                    continue

                if isinstance(cached_negative, dict):
                    negative_metadata_hits += 1
                    manifest_rejections.append(
                        {
                            "pdb_id": representative_id,
                            "reason": "negative_metadata_store",
                            "status_code": cached_negative.get("status_code"),
                            "detail": cached_negative.get("reason"),
                        }
                    )
                    continue

                metadata_misses += 1
                missing_ids.append(representative_id)

            (
                fetched_metadata,
                batch_errors,
                batch_negative_entries,
                batch_error_entries,
            ) = _fetch_metadata_for_ids(
                missing_ids,
                args,
                suppress_progress=True,
            )
            metadata_errors += batch_errors
            negative_metadata_updates += len(batch_negative_entries)
            manifest_rejections.extend(batch_error_entries)

            for representative_id, metadata in fetched_metadata.items():
                batch_metadata[representative_id] = metadata
                metadata_updates += 1
                _write_metadata_record(
                    resolved_metadata_cache_path,
                    representative_id,
                    entry_metadata=metadata,
                )

            for representative_id, negative_entry in batch_negative_entries.items():
                _write_metadata_record(
                    resolved_metadata_cache_path,
                    representative_id,
                    negative_entry=negative_entry,
                )

            for representative_id in batch_ids:
                if len(selected_ids) >= target_selection:
                    break

                examined_count += 1
                entry_metadata = batch_metadata.get(representative_id)
                if entry_metadata is None:
                    continue

                passes_filters, rejection_reason = rcsb.entry_passes_filters(
                    entry_metadata,
                    allowed_method_filters=cluster_method_filters,
                    max_resolution=args.cluster_max_resolution,
                    max_residues=args.cluster_max_residues,
                )
                if not passes_filters:
                    if rejection_reason == "experimental_method":
                        filtered_by_method += 1
                    elif rejection_reason == "resolution":
                        filtered_by_resolution += 1
                    elif rejection_reason == "missing_resolution":
                        filtered_missing_resolution += 1
                    elif rejection_reason == "residue_count":
                        filtered_by_residue_count += 1
                    manifest_rejections.append(
                        {
                            "pdb_id": representative_id,
                            "reason": rejection_reason,
                        }
                    )
                    continue

                selected_ids.append(representative_id)

            meta_processed = batch_start + len(batch_ids)
            meta_last_progress_emit = _emit_stage_progress(
                stage="metadata",
                completed=meta_processed,
                total=total_representatives,
                failed=metadata_errors,
                started_at=meta_progress_started_at,
                last_emitted_at=meta_last_progress_emit,
                interval=progress_interval,
            )

        _emit_stage_progress(
            stage="metadata",
            completed=meta_processed,
            total=total_representatives,
            failed=metadata_errors,
            started_at=meta_progress_started_at,
            last_emitted_at=meta_last_progress_emit,
            interval=progress_interval,
            force=True,
        )

        if args.write_manifest is not None:
            _write_cluster_manifest(
                manifest_path=Path(args.write_manifest),
                args=args,
                selected_ids=selected_ids,
                representative_to_member_ids=representative_to_member_ids,
                cluster_method_filters=cluster_method_filters,
                rejection_entries=manifest_rejections,
                examined_count=examined_count,
                filtered_by_method=filtered_by_method,
                filtered_by_resolution=filtered_by_resolution,
                filtered_missing_resolution=filtered_missing_resolution,
                filtered_by_residue_count=filtered_by_residue_count,
                metadata_errors=metadata_errors,
                metadata_hits=metadata_hits,
                negative_metadata_hits=negative_metadata_hits,
                metadata_misses=metadata_misses,
                metadata_updates=metadata_updates,
                negative_metadata_updates=negative_metadata_updates,
            )

        if len(selected_ids) == 0:
            raise RuntimeError(
                "No structures matched cluster filters. "
                f"Examined={examined_count}, "
                f"filtered_by_method={filtered_by_method}, "
                f"filtered_by_resolution={filtered_by_resolution}, "
                f"missing_resolution={filtered_missing_resolution}, "
                f"filtered_by_residue_count={filtered_by_residue_count}, "
                f"metadata_errors={metadata_errors}."
            )

        method_label = (
            "all"
            if cluster_method_filters is None
            else ",".join(cluster_method_filters)
        )
        print(
            "cluster selection summary: "
            f"selected={len(selected_ids)}, "
            f"examined={examined_count}, "
            f"shard_index={args.shard_index}, "
            f"num_shards={args.num_shards}, "
            f"methods={method_label}, "
            f"max_resolution={args.cluster_max_resolution}, "
            f"max_residues={args.cluster_max_residues}, "
            f"metadata_hits={metadata_hits}, "
            f"negative_metadata_hits={negative_metadata_hits}, "
            f"metadata_misses={metadata_misses}, "
            f"metadata_updates={metadata_updates}, "
            f"negative_metadata_updates={negative_metadata_updates}, "
            f"filtered_by_method={filtered_by_method}, "
            f"filtered_by_resolution={filtered_by_resolution}, "
            f"missing_resolution={filtered_missing_resolution}, "
            f"filtered_by_residue_count={filtered_by_residue_count}, "
            f"metadata_errors={metadata_errors}",
            file=sys.stderr,
        )

        selected_structure_entries: list[dict[str, object]] = []
        for representative_id in selected_ids:
            structure_entry: dict[str, object] = {
                "source": _sanitize_label(representative_id),
            }
            cluster_member_pdb_ids = representative_to_member_ids.get(representative_id)
            if cluster_member_pdb_ids is not None:
                structure_entry["cluster_member_pdb_ids"] = list(cluster_member_pdb_ids)
            selected_structure_entries.append(structure_entry)
        _set_structure_metadata(
            args,
            _build_structure_metadata_by_source(selected_structure_entries),
        )

        if args.dry_run:
            return [
                (_sanitize_label(representative_id), download_dir / f"{representative_id}.cif")
                for representative_id in selected_ids
            ]

        downloaded_paths = _download_cluster_structures(
            selected_ids,
            args,
            download_dir,
            output_dir,
        )

        return [
            (_sanitize_label(representative_id), downloaded_paths[representative_id])
            for representative_id in selected_ids
        ]

    raise RuntimeError("No input source selected.")


def _build_annotation_payload(
    source_label: str,
    input_path: Path,
    structure_output_path: Path,
    annotation_df,
    prepared_structure=None,
) -> dict:
    utils = _utils_module()
    volumes = json.loads(annotation_df.to_json(orient="records"))
    payload = {
        "source": source_label,
        "input_path": str(input_path),
        "output_structure_cif": str(structure_output_path),
        "backend": utils.get_active_backend(),
        "resolution": utils.VOXEL_SIZE,
        "num_volumes": len(volumes),
        "largest_type": volumes[0]["type"] if volumes else None,
        "largest_volume": volumes[0]["volume"] if volumes else None,
        "volumes": volumes,
    }
    if prepared_structure is not None:
        pdb = _pdb_module()
        payload.update(pdb.compute_sse_fractions(prepared_structure))
    return payload


def _filter_annotation_dataframe_for_output(annotation_df, include_hubs: bool):
    if include_hubs or "type" not in annotation_df.columns:
        return annotation_df

    keep_mask = annotation_df["type"].astype(str).str.lower() != "hub"
    if bool(keep_mask.all()):
        return annotation_df
    return annotation_df.loc[keep_mask].reset_index(drop=True)


def _filter_annotation_structure_for_output(annotation_structure, include_hubs: bool):
    if include_hubs or not hasattr(annotation_structure, "res_name"):
        return annotation_structure

    keep_indices = [
        index
        for index, res_name in enumerate(annotation_structure.res_name)
        if str(res_name) != "HUB"
    ]
    if len(keep_indices) == len(annotation_structure):
        return annotation_structure
    return annotation_structure[keep_indices]


def _infer_pdb_id_for_result(source_label: str, input_path: Path) -> str | None:
    for candidate in (source_label, input_path.stem):
        try:
            return rcsb.normalize_pdb_id(candidate)
        except ValueError:
            continue
    return None


def analyze_structure_file(
    source_label: str,
    input_path: Path,
    output_dir: Path,
    min_voxels: int,
    min_volume: float | None,
    overwrite: bool,
    assembly_policy: str = DEFAULT_ASSEMBLY_POLICY,
    max_residues: int | None = None,
    include_hubs: bool = False,
) -> dict:
    """
    Run volumizer on one structure file and write outputs.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    structure_output_path, annotation_output_path = _output_paths_for_label(
        output_dir,
        source_label,
    )

    _cleanup_partial_outputs(
        structure_output_path=structure_output_path,
        annotation_output_path=annotation_output_path,
    )
    _remove_status_record(output_dir, source_label)

    if not overwrite and structure_output_path.exists():
        raise FileExistsError(
            f"Output structure already exists: {structure_output_path}. Use --overwrite."
        )
    if not overwrite and annotation_output_path.exists():
        raise FileExistsError(
            f"Output annotation already exists: {annotation_output_path}. Use --overwrite."
        )

    pdb = _pdb_module()
    volumizer = _volumizer_module()

    input_structure = pdb.load_structure(
        input_path,
        assembly_policy=assembly_policy,
    )
    prepared_structure = volumizer.prepare_pdb_structure(input_structure)
    if max_residues is not None:
        prepared_residue_count = pdb.get_structure_residue_count(prepared_structure)
        if prepared_residue_count > int(max_residues):
            raise PostAssemblyResidueLimitExceeded(
                actual_residues=prepared_residue_count,
                max_residues=int(max_residues),
                assembly_policy=assembly_policy,
            )
    annotation_df, annotation_structure = volumizer.annotate_structure_volumes(
        prepared_structure,
        min_voxels=min_voxels,
        min_volume=min_volume,
    )
    output_annotation_df = _filter_annotation_dataframe_for_output(
        annotation_df,
        include_hubs=include_hubs,
    )
    output_annotation_structure = _filter_annotation_structure_for_output(
        annotation_structure,
        include_hubs=include_hubs,
    )

    pdb.ensure_b_factor_annotation(prepared_structure)
    pdb.ensure_b_factor_annotation(output_annotation_structure)
    combined_structure = prepared_structure + output_annotation_structure
    try:
        pdb.save_structure(combined_structure, structure_output_path)

        payload = _build_annotation_payload(
            source_label,
            input_path,
            structure_output_path,
            output_annotation_df,
            prepared_structure=prepared_structure,
        )
        annotation_output_path.write_text(
            json.dumps(payload, indent=2),
            encoding="utf-8",
        )
    except Exception:
        _cleanup_partial_outputs(
            structure_output_path=structure_output_path,
            annotation_output_path=annotation_output_path,
        )
        raise

    return {
        "source": source_label,
        "pdb_id": _infer_pdb_id_for_result(source_label, input_path),
        "input_path": str(input_path),
        "structure_output": str(structure_output_path),
        "annotation_output": str(annotation_output_path),
        "num_volumes": payload["num_volumes"],
        "largest_type": payload["largest_type"],
        "largest_volume": payload["largest_volume"],
    }


def _should_use_isolated_analysis_workers(
    *,
    args: argparse.Namespace,
    active_backend: str | None,
) -> bool:
    return args.command == "cluster" and active_backend == BACKEND_NATIVE


def _invoke_analysis_callable(
    fn,
    *,
    source_label: str,
    input_path: Path,
    output_dir: Path,
    min_voxels: int,
    min_volume: float | None,
    overwrite: bool,
    assembly_policy: str,
    max_residues: int | None = None,
    include_hubs: bool = False,
) -> dict:
    kwargs = {
        "source_label": source_label,
        "input_path": input_path,
        "output_dir": output_dir,
        "min_voxels": min_voxels,
        "min_volume": min_volume,
        "overwrite": overwrite,
        "assembly_policy": assembly_policy,
    }
    if max_residues is not None and _callable_accepts_keyword(fn, "max_residues"):
        kwargs["max_residues"] = int(max_residues)
    if _callable_accepts_keyword(fn, "include_hubs"):
        kwargs["include_hubs"] = bool(include_hubs)
    return fn(**kwargs)


def _build_analysis_worker_env() -> dict[str, str]:
    env = os.environ.copy()
    env.setdefault("PYTHONFAULTHANDLER", "1")
    for variable in _ANALYSIS_WORKER_THREAD_ENV_VARS:
        env.setdefault(variable, "1")
    return env


def _build_analysis_worker_command(
    *,
    source_label: str,
    input_path: Path,
    output_dir: Path,
    min_voxels: int,
    min_volume: float | None,
    overwrite: bool,
    assembly_policy: str,
    resolution: float,
    keep_non_protein: bool,
    backend: str | None,
    surface_connectivity: str,
    merge_mouth_gap_voxels: int,
    max_residues: int | None = None,
    include_hubs: bool = False,
) -> list[str]:
    command = [
        sys.executable,
        "-m",
        "volumizer._analysis_worker",
        "--source-label",
        source_label,
        "--input-path",
        str(input_path),
        "--output-dir",
        str(output_dir),
        "--min-voxels",
        str(int(min_voxels)),
        "--assembly-policy",
        assembly_policy,
        "--resolution",
        str(float(resolution)),
        "--surface-connectivity",
        surface_connectivity,
        "--merge-mouth-gap-voxels",
        str(int(merge_mouth_gap_voxels)),
    ]
    if min_volume is not None:
        command.extend(["--min-volume", str(float(min_volume))])
    if overwrite:
        command.append("--overwrite")
    if keep_non_protein:
        command.append("--keep-non-protein")
    if backend:
        command.extend(["--backend", backend])
    if max_residues is not None:
        command.extend(["--max-residues", str(int(max_residues))])
    if include_hubs:
        command.append("--include-hubs")
    return command


def _parse_analysis_worker_success(
    *,
    completed: subprocess.CompletedProcess[str],
    source_label: str,
    assembly_policy: str,
) -> dict:
    stdout = completed.stdout.strip()
    try:
        payload = json.loads(stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(
            f"analysis worker returned invalid JSON output for {source_label}"
        ) from error
    if not isinstance(payload, dict):
        raise RuntimeError("analysis worker returned non-object JSON output")
    if payload.get("status") == "skipped_post_assembly_residue_limit":
        raise PostAssemblyResidueLimitExceeded(
            actual_residues=int(payload["actual_residues"]),
            max_residues=int(payload["max_residues"]),
            assembly_policy=str(payload.get("assembly_policy") or assembly_policy),
        )
    return payload


def _coerce_subprocess_output_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def _format_worker_failure_message(
    *,
    source_label: str,
    completed: subprocess.CompletedProcess[str],
) -> str:
    stdout = completed.stdout.strip()
    stderr = completed.stderr.strip()
    if completed.returncode < 0:
        signal_number = -completed.returncode
        try:
            signal_name = signal.Signals(signal_number).name
            signal_label = f"{signal_number} ({signal_name})"
        except ValueError:
            signal_label = str(signal_number)
        detail = stderr or stdout or "worker terminated without stderr/stdout output"
        return (
            f"analysis worker terminated by signal {signal_label} for "
            f"{source_label}: {detail}"
        )

    detail = stderr or stdout or "worker exited without stderr/stdout output"
    return (
        f"analysis worker failed with exit code {completed.returncode} "
        f"for {source_label}: {detail}"
    )


def _format_worker_timeout_message(
    *,
    source_label: str,
    error: subprocess.TimeoutExpired,
) -> str:
    stdout = _coerce_subprocess_output_text(error.stdout).strip()
    stderr = _coerce_subprocess_output_text(error.stderr).strip()
    detail = stderr or stdout or "worker exceeded timeout without stderr/stdout output"
    return (
        f"analysis worker timed out after {float(error.timeout):.1f}s "
        f"for {source_label}: {detail}"
    )


def _run_isolated_analysis_worker(
    *,
    source_label: str,
    input_path: Path,
    output_dir: Path,
    min_voxels: int,
    min_volume: float | None,
    overwrite: bool,
    assembly_policy: str,
    resolution: float,
    keep_non_protein: bool,
    backend: str | None,
    surface_connectivity: str,
    merge_mouth_gap_voxels: int,
    max_residues: int | None = None,
    worker_timeout_seconds: float | None = None,
    include_hubs: bool = False,
) -> dict:
    env = _build_analysis_worker_env()
    backend_attempts: list[str | None] = [backend]
    if backend == BACKEND_NATIVE:
        backend_attempts.append(BACKEND_NATIVE)
        backend_attempts.append(BACKEND_PYTHON)

    failure_messages: list[str] = []

    for attempt_index, attempt_backend in enumerate(backend_attempts, start=1):
        try:
            completed = subprocess.run(
                _build_analysis_worker_command(
                    source_label=source_label,
                    input_path=input_path,
                    output_dir=output_dir,
                    min_voxels=min_voxels,
                    min_volume=min_volume,
                    overwrite=overwrite or attempt_index > 1,
                    assembly_policy=assembly_policy,
                    resolution=resolution,
                    keep_non_protein=keep_non_protein,
                    backend=attempt_backend,
                    surface_connectivity=surface_connectivity,
                    merge_mouth_gap_voxels=merge_mouth_gap_voxels,
                    max_residues=max_residues,
                    include_hubs=include_hubs,
                ),
                capture_output=True,
                text=True,
                check=False,
                env=env,
                timeout=worker_timeout_seconds,
            )
        except subprocess.TimeoutExpired as error:
            failure_messages.append(
                f"attempt {attempt_index}"
                + (
                    f" backend={attempt_backend}"
                    if attempt_backend is not None
                    else ""
                )
                + ": "
                + _format_worker_timeout_message(
                    source_label=source_label,
                    error=error,
                )
            )

            should_retry_timeout = (
                attempt_backend == BACKEND_NATIVE
                and attempt_index < len(backend_attempts)
            )
            if should_retry_timeout:
                continue
            break

        if completed.returncode == 0:
            return _parse_analysis_worker_success(
                completed=completed,
                source_label=source_label,
                assembly_policy=assembly_policy,
            )

        failure_messages.append(
            f"attempt {attempt_index}"
            + (
                f" backend={attempt_backend}"
                if attempt_backend is not None
                else ""
            )
            + ": "
            + _format_worker_failure_message(
                source_label=source_label,
                completed=completed,
            )
        )

        should_retry_signal = (
            completed.returncode < 0
            and attempt_backend == BACKEND_NATIVE
            and attempt_index < len(backend_attempts)
        )
        if should_retry_signal:
            continue

        break

    raise RuntimeError("; ".join(failure_messages))


def _plan_dry_run(
    structures: list[tuple[str, Path]],
    args: argparse.Namespace,
    output_dir: Path,
    tracker: _RunTracker,
) -> None:
    pending_structures = _prepare_pending_structures(
        structures=structures,
        args=args,
        output_dir=output_dir,
        tracker=tracker,
    )
    planned_now = 0

    for source_label, input_path, _ in pending_structures:
        tracker.mark_planned(
            _enrich_structure_entry(
                _build_dry_run_plan_entry(
                    source_label=source_label,
                    input_path=input_path,
                    output_dir=output_dir,
                ),
                args=args,
                source_label=source_label,
            )
        )
        planned_now += 1

    if planned_now > 0:
        print(
            f"dry-run selected {planned_now} structure(s); no downloads/analysis executed.",
            file=sys.stderr,
        )


def _analyze_structures(
    structures: list[tuple[str, Path]],
    args: argparse.Namespace,
    output_dir: Path,
    tracker: _RunTracker,
    active_backend: str | None = None,
) -> int:
    pending_structures = _prepare_pending_structures(
        structures=structures,
        args=args,
        output_dir=output_dir,
        tracker=tracker,
    )

    analysis_workers = int(args.jobs)
    if args.fail_fast and analysis_workers > 1:
        print(
            "--fail-fast enabled; forcing analysis worker count to 1.",
            file=sys.stderr,
        )
        analysis_workers = 1

    if len(pending_structures) == 0:
        return analysis_workers

    total_pending = len(pending_structures)
    use_isolated_workers = _should_use_isolated_analysis_workers(
        args=args,
        active_backend=active_backend,
    )
    runtime_max_residues = _get_post_assembly_max_residues(args)
    worker_timeout_seconds = _get_analysis_worker_timeout_seconds(
        args,
        use_isolated_workers=use_isolated_workers,
    )
    progress_interval = float(args.progress_interval)
    progress_started_at = time.monotonic()
    last_progress_emit = progress_started_at

    def _active_source_labels(
        labels: Sequence[str],
    ) -> list[str] | None:
        if len(labels) == 0 or len(labels) > 3:
            return None
        return list(labels)

    if use_isolated_workers:
        worker_mode_message = "analysis worker mode: isolated subprocesses"
        if worker_timeout_seconds is not None:
            worker_mode_message += f" (timeout={worker_timeout_seconds:.0f}s)"
        print(worker_mode_message, file=sys.stderr)

    if analysis_workers <= 1:
        completed = 0
        failed = 0
        for source_label, input_path, overwrite_existing_outputs in pending_structures:
            tracker.mark_started(source_label, input_path)
            print(f"analyzing {source_label}: {input_path}", file=sys.stderr)
            try:
                if use_isolated_workers:
                    result = _run_isolated_analysis_worker(
                        source_label=source_label,
                        input_path=input_path,
                        output_dir=output_dir,
                        min_voxels=int(args.min_voxels),
                        min_volume=args.min_volume,
                        overwrite=bool(
                            args.overwrite
                            or args.reannotate
                            or overwrite_existing_outputs
                        ),
                        assembly_policy=str(args.assembly_policy),
                        resolution=float(args.resolution),
                        keep_non_protein=bool(args.keep_non_protein),
                        backend=active_backend,
                        surface_connectivity=str(args.surface_connectivity),
                        merge_mouth_gap_voxels=int(args.merge_mouth_gap_voxels),
                        max_residues=runtime_max_residues,
                        worker_timeout_seconds=worker_timeout_seconds,
                        include_hubs=bool(args.include_hubs),
                    )
                else:
                    result = _invoke_analysis_callable(
                        analyze_structure_file,
                        source_label=source_label,
                        input_path=input_path,
                        output_dir=output_dir,
                        min_voxels=int(args.min_voxels),
                        min_volume=args.min_volume,
                        overwrite=bool(
                            args.overwrite
                            or args.reannotate
                            or overwrite_existing_outputs
                        ),
                        assembly_policy=str(args.assembly_policy),
                        max_residues=runtime_max_residues,
                        include_hubs=bool(args.include_hubs),
                    )
                tracker.mark_result(
                    _enrich_structure_entry(
                        result,
                        args=args,
                        source_label=source_label,
                    )
                )
            except PostAssemblyResidueLimitExceeded as error:
                _record_post_assembly_residue_limit_skip(
                    tracker=tracker,
                    source_label=source_label,
                    input_path=input_path,
                    output_dir=output_dir,
                    error=error,
                    extra_entry_fields=_get_structure_metadata(args, source_label),
                )
                print(f"skipping {source_label}: {error}", file=sys.stderr)
                completed += 1
                last_progress_emit = _emit_stage_progress(
                    stage="analysis",
                    completed=completed,
                    total=total_pending,
                    failed=failed,
                    started_at=progress_started_at,
                    last_emitted_at=last_progress_emit,
                    interval=progress_interval,
                    active=1 if completed < total_pending else 0,
                    queued=max(total_pending - completed - 1, 0),
                )
                continue
            except Exception as error:  # pragma: no cover - exercised via CLI integration tests
                failed += 1
                error_entry = _enrich_structure_entry(
                    _build_error_entry(
                        source_label=source_label,
                        input_path=input_path,
                        error=error,
                    ),
                    args=args,
                    source_label=source_label,
                )
                tracker.mark_error(error_entry)
                print(
                    f"error for {source_label}: {error_entry['error']}",
                    file=sys.stderr,
                )
                completed += 1
                last_progress_emit = _emit_stage_progress(
                    stage="analysis",
                    completed=completed,
                    total=total_pending,
                    failed=failed,
                    started_at=progress_started_at,
                    last_emitted_at=last_progress_emit,
                    interval=progress_interval,
                    active=1 if completed < total_pending else 0,
                    queued=max(total_pending - completed - 1, 0),
                )
                if args.fail_fast:
                    break
                continue

            completed += 1
            last_progress_emit = _emit_stage_progress(
                stage="analysis",
                completed=completed,
                total=total_pending,
                failed=failed,
                started_at=progress_started_at,
                last_emitted_at=last_progress_emit,
                interval=progress_interval,
                active=1 if completed < total_pending else 0,
                queued=max(total_pending - completed - 1, 0),
            )

        _emit_stage_progress(
            stage="analysis",
            completed=completed,
            total=total_pending,
            failed=failed,
            started_at=progress_started_at,
            last_emitted_at=last_progress_emit,
            interval=progress_interval,
            active=0,
            queued=0,
            force=True,
        )
        return analysis_workers

    print(
        f"analyzing {total_pending} structures with {analysis_workers} workers...",
        file=sys.stderr,
    )
    with ThreadPoolExecutor(max_workers=analysis_workers) as executor:
        future_to_context = {
            executor.submit(
                (
                    _run_isolated_analysis_worker
                    if use_isolated_workers
                    else _invoke_analysis_callable
                ),
                **(
                    {
                        "source_label": source_label,
                        "input_path": input_path,
                        "output_dir": output_dir,
                        "min_voxels": int(args.min_voxels),
                        "min_volume": args.min_volume,
                        "overwrite": bool(
                            args.overwrite
                            or args.reannotate
                            or overwrite_existing_outputs
                        ),
                        "assembly_policy": str(args.assembly_policy),
                        "resolution": float(args.resolution),
                        "keep_non_protein": bool(args.keep_non_protein),
                        "backend": active_backend,
                        "surface_connectivity": str(args.surface_connectivity),
                        "merge_mouth_gap_voxels": int(args.merge_mouth_gap_voxels),
                        "max_residues": runtime_max_residues,
                        "worker_timeout_seconds": worker_timeout_seconds,
                        "include_hubs": bool(args.include_hubs),
                    }
                    if use_isolated_workers
                    else {
                        "fn": analyze_structure_file,
                        "source_label": source_label,
                        "input_path": input_path,
                        "output_dir": output_dir,
                        "min_voxels": int(args.min_voxels),
                        "min_volume": args.min_volume,
                        "overwrite": bool(
                            args.overwrite
                            or args.reannotate
                            or overwrite_existing_outputs
                        ),
                        "assembly_policy": str(args.assembly_policy),
                        "max_residues": runtime_max_residues,
                        "include_hubs": bool(args.include_hubs),
                    }
                ),
            ): (index, source_label, input_path)
            for index, (
                source_label,
                input_path,
                overwrite_existing_outputs,
            ) in enumerate(pending_structures)
        }

        for _, source_label, input_path in future_to_context.values():
            tracker.mark_started(source_label, input_path)

        ordered_results: dict[int, dict] = {}
        ordered_skips: dict[int, dict] = {}
        ordered_errors: dict[int, dict] = {}
        completed = 0
        failed = 0
        pending_futures = set(future_to_context)
        last_progress_emit = _emit_stage_progress(
            stage="analysis",
            completed=completed,
            total=total_pending,
            failed=failed,
            started_at=progress_started_at,
            last_emitted_at=last_progress_emit,
            interval=progress_interval,
            active=min(analysis_workers, total_pending),
            queued=max(total_pending - min(analysis_workers, total_pending), 0),
            active_sources=_active_source_labels(
                [source_label for _, source_label, _ in future_to_context.values()]
            ),
            force=True,
        )

        wait_timeout = max(0.1, progress_interval) if progress_interval > 0 else None

        while pending_futures:
            done, not_done = wait(
                pending_futures,
                timeout=wait_timeout,
                return_when=FIRST_COMPLETED,
            )
            if len(done) == 0:
                active = min(len(pending_futures), analysis_workers)
                queued = max(total_pending - completed - active, 0)
                active_sources = _active_source_labels(
                    [future_to_context[future][1] for future in pending_futures]
                )
                last_progress_emit = _emit_stage_progress(
                    stage="analysis",
                    completed=completed,
                    total=total_pending,
                    failed=failed,
                    started_at=progress_started_at,
                    last_emitted_at=last_progress_emit,
                    interval=progress_interval,
                    active=active,
                    queued=queued,
                    active_sources=active_sources,
                    force=True,
                )
                continue

            pending_futures = set(not_done)
            for future in done:
                index, source_label, input_path = future_to_context[future]
                try:
                    ordered_results[index] = _enrich_structure_entry(
                        future.result(),
                        args=args,
                        source_label=source_label,
                    )
                except PostAssemblyResidueLimitExceeded as error:
                    ordered_skips[index] = _enrich_structure_entry(
                        _build_post_assembly_residue_limit_skip_entry(
                            source_label=source_label,
                            input_path=input_path,
                            output_dir=output_dir,
                            error=error,
                        ),
                        args=args,
                        source_label=source_label,
                    )
                    _write_post_assembly_residue_limit_status_record(
                        source_label=source_label,
                        input_path=input_path,
                        output_dir=output_dir,
                        error=error,
                    )
                    print(f"skipping {source_label}: {error}", file=sys.stderr)
                except Exception as error:  # pragma: no cover - exercised via CLI integration tests
                    failed += 1
                    ordered_errors[index] = _enrich_structure_entry(
                        _build_error_entry(
                            source_label=source_label,
                            input_path=input_path,
                            error=error,
                        ),
                        args=args,
                        source_label=source_label,
                    )
                    print(
                        f"error for {source_label}: {ordered_errors[index]['error']}",
                        file=sys.stderr,
                    )

                completed += 1

            active = min(len(pending_futures), analysis_workers)
            queued = max(total_pending - completed - active, 0)
            active_sources = _active_source_labels(
                [future_to_context[future][1] for future in pending_futures]
            )
            last_progress_emit = _emit_stage_progress(
                stage="analysis",
                completed=completed,
                total=total_pending,
                failed=failed,
                started_at=progress_started_at,
                last_emitted_at=last_progress_emit,
                interval=progress_interval,
                active=active,
                queued=queued,
                active_sources=active_sources,
            )

    _emit_stage_progress(
        stage="analysis",
        completed=completed,
        total=total_pending,
        failed=failed,
        started_at=progress_started_at,
        last_emitted_at=last_progress_emit,
        interval=progress_interval,
        active=0,
        queued=0,
        force=True,
    )

    for index in range(total_pending):
        if index in ordered_results:
            tracker.mark_result(ordered_results[index], persist=False)
        elif index in ordered_skips:
            tracker.mark_skipped(ordered_skips[index], persist=False)
        elif index in ordered_errors:
            tracker.mark_error(ordered_errors[index], persist=False)

    tracker.persist_checkpoint()

    return analysis_workers


def _run_cache_command(args: argparse.Namespace) -> int:
    cache_path = Path(args.metadata_cache)

    if args.cache_command == "inspect":
        entries, negative_entries = _load_metadata_cache(cache_path)
        payload = {
            "metadata_store": str(cache_path),
            "exists": cache_path.is_dir(),
            "entries": len(entries),
            "negative_entries": len(negative_entries),
        }
        print(json.dumps(payload, indent=2))
        return 0

    if args.cache_command == "clear-negative":
        removed = _clear_negative_metadata_records(cache_path)
        print(
            f"cleared {removed} negative metadata record{'s' if removed != 1 else ''}: {cache_path}",
            file=sys.stderr,
        )
        return 0

    raise ValueError(f"Unknown cache command: {args.cache_command}")


def _run_analysis_command(args: argparse.Namespace) -> int:
    if args.jobs < 1:
        raise ValueError("--jobs must be >= 1.")
    if args.retries < 0:
        raise ValueError("--retries must be >= 0.")
    if args.retry_delay < 0:
        raise ValueError("--retry-delay must be >= 0.")
    if args.worker_timeout_seconds < 0:
        raise ValueError("--worker-timeout-seconds must be >= 0.")
    if args.progress_interval < 0:
        raise ValueError("--progress-interval must be >= 0.")
    if args.assembly_policy not in VALID_ASSEMBLY_POLICIES:
        raise ValueError(
            "--assembly-policy must be one of: "
            + ", ".join(VALID_ASSEMBLY_POLICIES)
        )

    if args.command == "cluster":
        has_num_shards = args.num_shards is not None
        has_shard_index = args.shard_index is not None
        if has_num_shards != has_shard_index:
            raise ValueError("Use --num-shards and --shard-index together.")
        if has_num_shards:
            if args.num_shards < 1:
                raise ValueError("--num-shards must be >= 1.")
            if args.shard_index < 0:
                raise ValueError("--shard-index must be >= 0.")
            if args.shard_index >= args.num_shards:
                raise ValueError("--shard-index must be < --num-shards.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    download_dir = (
        Path(args.download_dir)
        if args.download_dir is not None
        else output_dir / "downloads"
    )

    cluster_method_filters = None
    metadata_cache_path = None
    if args.command == "cluster":
        cluster_method_filters = _resolve_cluster_method_filters(args)
        metadata_cache_path = _resolve_metadata_cache_path(args, output_dir)

    checkpoint_path = _resolve_checkpoint_path(args, output_dir)
    progress_jsonl_path = _resolve_progress_jsonl_path(args)

    tracker = _RunTracker(
        checkpoint_path=checkpoint_path,
        progress_jsonl_path=progress_jsonl_path,
        signature=_make_checkpoint_signature(
            args,
            output_dir=output_dir,
            download_dir=download_dir,
            cluster_method_filters=cluster_method_filters,
            metadata_cache_path=metadata_cache_path,
        ),
        resume=args.resume,
    )

    if args.command == "cluster":
        structures = resolve_input_structures(
            args,
            download_dir,
            output_dir,
            metadata_cache_path=metadata_cache_path,
        )
    else:
        structures = resolve_input_structures(args, download_dir, output_dir)
    tracker.mark_run_start(total_structures=len(structures), dry_run=args.dry_run)

    active_backend = None
    if args.dry_run:
        _plan_dry_run(structures, args, output_dir, tracker)
        analysis_workers = 0
    else:
        native_backend = _native_backend_module()
        utils = _utils_module()

        if args.backend is not None:
            os.environ[native_backend.BACKEND_ENV] = args.backend
            native_backend.clear_backend_cache()

        utils.set_resolution(float(args.resolution))
        utils.set_non_protein(bool(args.keep_non_protein))
        utils.set_surface_component_connectivity_mode(
            str(args.surface_connectivity)
        )
        utils.set_surface_mouth_merge_gap_voxels(
            int(args.merge_mouth_gap_voxels)
        )

        active_backend = native_backend.active_backend()
        print(f"[volumizer] backend: {active_backend}", file=sys.stderr)
        analysis_workers = _analyze_structures(
            structures,
            args,
            output_dir,
            tracker,
            active_backend=active_backend,
        )

    summary = {
        "config": {
            "command": args.command,
            "input": str(args.input) if args.input is not None else None,
            "pdb_id": args.pdb_id,
            "manifest": str(args.manifest) if args.manifest is not None else None,
            "from_summary": (
                str(args.from_summary) if args.from_summary is not None else None
            ),
            "from_summary_only": args.only if args.from_summary is not None else None,
            "cluster_identity": args.cluster_identity,
            "max_structures": args.max_structures,
            "num_shards": args.num_shards,
            "shard_index": args.shard_index,
            "cluster_method_filters": cluster_method_filters,
            "cluster_allow_all_methods": args.cluster_allow_all_methods,
            "cluster_max_resolution": args.cluster_max_resolution,
            "cluster_max_residues": args.cluster_max_residues,
            "metadata_cache": str(metadata_cache_path) if metadata_cache_path else None,
            "no_metadata_cache": args.no_metadata_cache,
            "write_manifest": (
                str(args.write_manifest) if args.write_manifest is not None else None
            ),
            "failures_manifest": (
                str(args.failures_manifest) if args.failures_manifest is not None else None
            ),
            "checkpoint": str(checkpoint_path) if checkpoint_path else None,
            "no_checkpoint": args.no_checkpoint,
            "progress_jsonl": str(progress_jsonl_path) if progress_jsonl_path else None,
            "progress_interval": args.progress_interval,
            "output_dir": str(output_dir),
            "download_dir": str(download_dir),
            "resolution": args.resolution,
            "assembly_policy": args.assembly_policy,
            "min_voxels": args.min_voxels,
            "min_volume": args.min_volume,
            "surface_connectivity": args.surface_connectivity,
            "merge_mouth_gap_voxels": args.merge_mouth_gap_voxels,
            "backend": (
                active_backend
                if active_backend is not None
                else (args.backend if args.backend is not None else _requested_backend_label())
            ),
            "keep_non_protein": args.keep_non_protein,
            "include_hubs": args.include_hubs,
            "jobs": args.jobs,
            "analysis_workers": analysis_workers,
            "timeout": args.timeout,
            "worker_timeout_seconds": args.worker_timeout_seconds,
            "retries": args.retries,
            "retry_delay": args.retry_delay,
            "overwrite": args.overwrite,
            "reannotate": args.reannotate,
            "resume": args.resume,
            "dry_run": args.dry_run,
        },
        "num_processed": len(tracker.results),
        "num_failed": len(tracker.errors),
        "num_skipped": len(tracker.skipped),
        "num_planned": len(tracker.planned),
        "results": tracker.results,
        "errors": tracker.errors,
        "skipped": tracker.skipped,
        "planned": tracker.planned,
    }

    summary_path = output_dir / "run.summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote summary: {summary_path}", file=sys.stderr)

    if args.failures_manifest is not None:
        _write_failures_manifest(
            manifest_path=Path(args.failures_manifest),
            command=args.command,
            errors=tracker.errors,
            summary_path=summary_path,
        )

    exit_code = 1 if len(tracker.errors) > 0 else 0
    tracker.mark_run_complete(exit_code)
    return exit_code


def run_cli(args: argparse.Namespace) -> int:
    """
    Execute CLI command and return process status code.
    """
    if args.command == "cache":
        return _run_cache_command(args)
    if args.command in {"analyze", "cluster"}:
        return _run_analysis_command(args)
    raise ValueError("No subcommand selected. Use: analyze, cluster, or cache.")


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()

    raw_argv = list(sys.argv[1:] if argv is None else argv)
    normalized_argv = _normalize_argv_for_subcommands(raw_argv)

    args = parser.parse_args(normalized_argv)
    try:
        return run_cli(args)
    except ValueError as error:
        parser.error(str(error))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
