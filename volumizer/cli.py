"""
Command-line interface for volumizer.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import sys
from typing import Sequence

from volumizer import native_backend, pdb, rcsb, utils, volumizer


DEFAULT_METADATA_CACHE_FILENAME = "entry_metadata_cache.json"
METADATA_CACHE_FORMAT_VERSION = 2
DEFAULT_CHECKPOINT_FILENAME = "run.checkpoint.json"


def _utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


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

    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--input",
        type=Path,
        help="Path to a local input structure file (.pdb, .cif, .mmtf, etc.).",
    )
    source_group.add_argument(
        "--pdb-id",
        type=str,
        help="Single PDB ID to download from RCSB and analyze.",
    )
    source_group.add_argument(
        "--cluster-identity",
        type=int,
        choices=sorted(rcsb.VALID_CLUSTER_IDENTITIES),
        help=(
            "RCSB sequence identity threshold for representative cluster download "
            "(e.g. 30, 50, 90)."
        ),
    )

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
        "--max-structures",
        type=int,
        default=None,
        help="Optional cap when using --cluster-identity.",
    )

    cluster_method_group = parser.add_mutually_exclusive_group()
    cluster_method_group.add_argument(
        "--cluster-method",
        action="append",
        default=None,
        metavar="METHOD",
        help=(
            "Allowed method filter(s) for --cluster-identity. Repeatable. "
            "Supported aliases: xray|x-ray, em|cryo-em, nmr, neutron. "
            "Default: xray + em."
        ),
    )
    cluster_method_group.add_argument(
        "--cluster-allow-all-methods",
        action="store_true",
        help="Disable method filtering when using --cluster-identity.",
    )
    parser.add_argument(
        "--cluster-max-resolution",
        type=float,
        default=None,
        help=(
            "Optional max best-resolution (Angstrom) filter for --cluster-identity. "
            "Entries missing resolution are excluded when this is set."
        ),
    )
    parser.add_argument(
        "--metadata-cache",
        type=Path,
        default=None,
        help=(
            "Path for cluster entry-metadata cache JSON "
            f"(default: <output-dir>/{DEFAULT_METADATA_CACHE_FILENAME})."
        ),
    )
    parser.add_argument(
        "--no-metadata-cache",
        action="store_true",
        help="Disable metadata-cache read/write for --cluster-identity runs.",
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
        default=2,
        help="Minimum voxels cutoff for reporting volumes (default: 2).",
    )
    parser.add_argument(
        "--min-volume",
        type=float,
        default=None,
        help="Optional minimum volume cutoff for reporting.",
    )
    parser.add_argument(
        "--backend",
        choices=sorted(native_backend.VALID_BACKENDS),
        default=None,
        help="Override backend mode for this run (python|auto|native).",
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
        "--resume",
        action="store_true",
        help="Skip structures that already have both output files.",
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

    return parser


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


def _has_complete_outputs(output_dir: Path, source_label: str) -> bool:
    structure_output_path, annotation_output_path = _output_paths_for_label(
        output_dir,
        source_label,
    )
    return structure_output_path.exists() and annotation_output_path.exists()


def _build_resume_skip_entry(
    source_label: str,
    input_path: Path,
    output_dir: Path,
) -> dict:
    structure_output_path, annotation_output_path = _output_paths_for_label(
        output_dir,
        source_label,
    )

    payload = {}
    if annotation_output_path.is_file():
        try:
            payload = json.loads(annotation_output_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            payload = {}

    return {
        "source": source_label,
        "input_path": str(input_path),
        "structure_output": str(structure_output_path),
        "annotation_output": str(annotation_output_path),
        "reason": "resume_existing_outputs",
        "num_volumes": payload.get("num_volumes"),
        "largest_type": payload.get("largest_type"),
        "largest_volume": payload.get("largest_volume"),
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


def _resolve_metadata_cache_path(
    args: argparse.Namespace,
    output_dir: Path,
) -> Path | None:
    if args.no_metadata_cache:
        return None

    if args.metadata_cache is not None:
        return Path(args.metadata_cache)

    return output_dir / DEFAULT_METADATA_CACHE_FILENAME


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


def _make_checkpoint_signature(
    args: argparse.Namespace,
    output_dir: Path,
    download_dir: Path,
    cluster_method_filters: list[str] | None,
    metadata_cache_path: Path | None,
) -> dict:
    return {
        "input": str(args.input) if args.input is not None else None,
        "pdb_id": args.pdb_id,
        "cluster_identity": args.cluster_identity,
        "max_structures": args.max_structures,
        "cluster_method_filters": cluster_method_filters,
        "cluster_allow_all_methods": args.cluster_allow_all_methods,
        "cluster_max_resolution": args.cluster_max_resolution,
        "metadata_cache": str(metadata_cache_path) if metadata_cache_path else None,
        "output_dir": str(output_dir),
        "download_dir": str(download_dir),
        "resolution": args.resolution,
        "min_voxels": args.min_voxels,
        "min_volume": args.min_volume,
        "backend": args.backend,
        "keep_non_protein": args.keep_non_protein,
        "dry_run": args.dry_run,
    }


def _load_metadata_cache(cache_path: Path | None) -> tuple[dict[str, dict], dict[str, dict]]:
    if cache_path is None or not cache_path.is_file():
        return {}, {}

    try:
        raw_payload = json.loads(cache_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}, {}

    if not isinstance(raw_payload, dict):
        return {}, {}

    if "entries" in raw_payload or "negative_entries" in raw_payload:
        raw_entries = raw_payload.get("entries")
        raw_negative_entries = raw_payload.get("negative_entries")
    else:
        raw_entries = raw_payload
        raw_negative_entries = {}

    if not isinstance(raw_entries, dict):
        raw_entries = {}
    if not isinstance(raw_negative_entries, dict):
        raw_negative_entries = {}

    entries: dict[str, dict] = {}
    negative_entries: dict[str, dict] = {}

    for key, value in raw_entries.items():
        if not isinstance(key, str) or not isinstance(value, dict):
            continue
        try:
            normalized_key = rcsb.normalize_pdb_id(key)
        except ValueError:
            continue
        entries[normalized_key] = value

    for key, value in raw_negative_entries.items():
        if not isinstance(key, str) or not isinstance(value, dict):
            continue
        try:
            normalized_key = rcsb.normalize_pdb_id(key)
        except ValueError:
            continue
        negative_entries[normalized_key] = value

    return entries, negative_entries


def _save_metadata_cache(
    cache_path: Path,
    entries: dict[str, dict],
    negative_entries: dict[str, dict],
) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "cache_format": METADATA_CACHE_FORMAT_VERSION,
        "entries": dict(sorted(entries.items())),
        "negative_entries": dict(sorted(negative_entries.items())),
    }
    cache_path.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )


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
        self._load_checkpoint_if_resuming()

    def _init_progress_stream(self) -> None:
        if self.progress_jsonl_path is None:
            return

        self.progress_jsonl_path.parent.mkdir(parents=True, exist_ok=True)
        if not self.resume:
            self.progress_jsonl_path.write_text("", encoding="utf-8")

    def _load_checkpoint_if_resuming(self) -> None:
        if not self.resume:
            return
        if self.checkpoint_path is None or not self.checkpoint_path.is_file():
            return

        try:
            payload = json.loads(self.checkpoint_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            return

        if not isinstance(payload, dict):
            return

        checkpoint_signature = payload.get("signature")
        if checkpoint_signature != self.signature:
            print(
                "checkpoint signature mismatch; starting with empty run-state.",
                file=sys.stderr,
            )
            return

        for kind in ("results", "errors", "skipped", "planned"):
            raw_list = payload.get(kind, [])
            if not isinstance(raw_list, list):
                continue
            for entry in raw_list:
                if not isinstance(entry, dict):
                    continue
                source = entry.get("source")
                if not isinstance(source, str):
                    continue
                self._set_entry(kind, entry)

        loaded_count = (
            len(self.results)
            + len(self.errors)
            + len(self.skipped)
            + len(self.planned)
        )
        if loaded_count > 0:
            print(
                (
                    "loaded checkpoint state: "
                    f"results={len(self.results)}, "
                    f"errors={len(self.errors)}, "
                    f"skipped={len(self.skipped)}, "
                    f"planned={len(self.planned)}"
                ),
                file=sys.stderr,
            )

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

    def mark_skipped(self, entry: dict) -> None:
        self._set_entry("skipped", entry)
        self.emit_event(
            "structure_skipped",
            source=entry["source"],
            input_path=entry["input_path"],
            reason=entry.get("reason"),
        )
        self.persist_checkpoint()

    def mark_started(self, source: str, input_path: Path) -> None:
        self.emit_event(
            "structure_started",
            source=source,
            input_path=str(input_path),
        )

    def mark_result(self, entry: dict) -> None:
        self._set_entry("results", entry)
        self.emit_event(
            "structure_succeeded",
            source=entry["source"],
            input_path=entry["input_path"],
            num_volumes=entry.get("num_volumes"),
        )
        self.persist_checkpoint()

    def mark_error(self, entry: dict) -> None:
        self._set_entry("errors", entry)
        self.emit_event(
            "structure_failed",
            source=entry["source"],
            input_path=entry["input_path"],
            error=entry["error"],
        )
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
        "error": str(error),
    }


def _fetch_metadata_for_ids(
    representative_ids: list[str],
    args: argparse.Namespace,
    negative_entries: dict[str, dict],
) -> tuple[dict[str, dict], int, int]:
    fetched: dict[str, dict] = {}
    errors = 0
    negative_updates = 0

    if len(representative_ids) == 0:
        return fetched, errors, negative_updates

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
                    negative_entries[representative_id] = _build_negative_cache_entry(
                        error,
                        status_code,
                    )
                    negative_updates += 1
                print(
                    f"metadata error for {representative_id}: {error}",
                    file=sys.stderr,
                )
                if args.fail_fast:
                    raise
        return fetched, errors, negative_updates

    max_workers = min(args.jobs, len(representative_ids))
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
                    negative_entries[representative_id] = _build_negative_cache_entry(
                        error,
                        status_code,
                    )
                    negative_updates += 1
                print(
                    f"metadata error for {representative_id}: {error}",
                    file=sys.stderr,
                )
                if args.fail_fast:
                    raise

    return fetched, errors, negative_updates


def _download_cluster_structures(
    selected_ids: list[str],
    args: argparse.Namespace,
    download_dir: Path,
    output_dir: Path,
) -> dict[str, Path]:
    ids_to_download = []
    for representative_id in selected_ids:
        source_label = _sanitize_label(representative_id)
        if args.resume and _has_complete_outputs(output_dir, source_label):
            continue
        ids_to_download.append(representative_id)

    if len(ids_to_download) == 0:
        return {}

    downloaded_paths: dict[str, Path] = {}

    if args.jobs <= 1:
        for index, representative_id in enumerate(ids_to_download, start=1):
            print(
                f"[{index}/{len(ids_to_download)}] downloading {representative_id}...",
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
        return downloaded_paths

    print(
        f"downloading {len(ids_to_download)} structures with {args.jobs} workers...",
        file=sys.stderr,
    )
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
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

        completed = 0
        for future in as_completed(future_to_id):
            representative_id = future_to_id[future]
            downloaded_paths[representative_id] = future.result()
            completed += 1
            print(
                f"[{completed}/{len(ids_to_download)}] downloaded {representative_id}",
                file=sys.stderr,
            )

    return downloaded_paths


def resolve_input_structures(
    args: argparse.Namespace,
    download_dir: Path,
    output_dir: Path,
) -> list[tuple[str, Path]]:
    """
    Resolve CLI source arguments into labeled local structure paths.
    """
    if args.input is not None:
        input_path = Path(args.input)
        if not input_path.is_file():
            raise FileNotFoundError(f"Input structure does not exist: {input_path}")
        return [(_sanitize_label(input_path.stem), input_path)]

    if args.pdb_id is not None:
        normalized_id = rcsb.normalize_pdb_id(args.pdb_id)
        source_label = _sanitize_label(normalized_id)

        if args.resume and _has_complete_outputs(output_dir, source_label):
            return [(source_label, download_dir / f"{normalized_id}.cif")]

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
        metadata_cache_path = _resolve_metadata_cache_path(args, output_dir)
        metadata_entries, negative_entries = _load_metadata_cache(metadata_cache_path)

        representative_ids = rcsb.fetch_cluster_representative_entry_ids(
            identity=args.cluster_identity,
            max_structures=None,
            timeout=args.timeout,
            include_non_pdb=False,
            retries=args.retries,
            retry_delay=args.retry_delay,
        )
        if len(representative_ids) == 0:
            raise RuntimeError(
                f"No representative PDB IDs found for cluster identity {args.cluster_identity}."
            )

        target_selection = len(representative_ids)
        if args.max_structures is not None and args.max_structures > 0:
            target_selection = args.max_structures

        selected_ids: list[str] = []
        examined_count = 0
        filtered_by_method = 0
        filtered_by_resolution = 0
        filtered_missing_resolution = 0
        metadata_errors = 0
        cache_hits = 0
        negative_cache_hits = 0
        cache_misses = 0
        cache_updates = 0
        negative_cache_updates = 0

        batch_size = max(1, int(args.jobs))
        for batch_start in range(0, len(representative_ids), batch_size):
            if len(selected_ids) >= target_selection:
                break

            batch_ids = representative_ids[batch_start : batch_start + batch_size]
            batch_metadata: dict[str, dict] = {}
            missing_ids: list[str] = []

            for representative_id in batch_ids:
                cached_metadata = metadata_entries.get(representative_id)
                if isinstance(cached_metadata, dict):
                    batch_metadata[representative_id] = cached_metadata
                    cache_hits += 1
                    continue

                cached_negative = negative_entries.get(representative_id)
                if isinstance(cached_negative, dict):
                    negative_cache_hits += 1
                    continue

                cache_misses += 1
                missing_ids.append(representative_id)

            (
                fetched_metadata,
                batch_errors,
                batch_negative_updates,
            ) = _fetch_metadata_for_ids(
                missing_ids,
                args,
                negative_entries,
            )
            metadata_errors += batch_errors
            negative_cache_updates += batch_negative_updates

            for representative_id, metadata in fetched_metadata.items():
                batch_metadata[representative_id] = metadata
                metadata_entries[representative_id] = metadata
                cache_updates += 1

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
                )
                if not passes_filters:
                    if rejection_reason == "experimental_method":
                        filtered_by_method += 1
                    elif rejection_reason == "resolution":
                        filtered_by_resolution += 1
                    elif rejection_reason == "missing_resolution":
                        filtered_missing_resolution += 1
                    continue

                selected_ids.append(representative_id)

        if metadata_cache_path is not None and (
            cache_updates > 0 or negative_cache_updates > 0
        ):
            _save_metadata_cache(
                metadata_cache_path,
                entries=metadata_entries,
                negative_entries=negative_entries,
            )

        if len(selected_ids) == 0:
            raise RuntimeError(
                "No structures matched cluster filters. "
                f"Examined={examined_count}, "
                f"filtered_by_method={filtered_by_method}, "
                f"filtered_by_resolution={filtered_by_resolution}, "
                f"missing_resolution={filtered_missing_resolution}, "
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
            f"methods={method_label}, "
            f"max_resolution={args.cluster_max_resolution}, "
            f"cache_hits={cache_hits}, "
            f"negative_cache_hits={negative_cache_hits}, "
            f"cache_misses={cache_misses}, "
            f"cache_updates={cache_updates}, "
            f"negative_cache_updates={negative_cache_updates}, "
            f"filtered_by_method={filtered_by_method}, "
            f"filtered_by_resolution={filtered_by_resolution}, "
            f"missing_resolution={filtered_missing_resolution}, "
            f"metadata_errors={metadata_errors}",
            file=sys.stderr,
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

        structures = []
        for representative_id in selected_ids:
            source_label = _sanitize_label(representative_id)
            if args.resume and _has_complete_outputs(output_dir, source_label):
                structures.append((source_label, download_dir / f"{representative_id}.cif"))
            else:
                structures.append((source_label, downloaded_paths[representative_id]))

        return structures

    raise RuntimeError("No input source selected.")


def _build_annotation_payload(
    source_label: str,
    input_path: Path,
    structure_output_path: Path,
    annotation_df,
) -> dict:
    volumes = json.loads(annotation_df.to_json(orient="records"))
    return {
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


def analyze_structure_file(
    source_label: str,
    input_path: Path,
    output_dir: Path,
    min_voxels: int,
    min_volume: float | None,
    overwrite: bool,
) -> dict:
    """
    Run volumizer on one structure file and write outputs.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    structure_output_path, annotation_output_path = _output_paths_for_label(
        output_dir,
        source_label,
    )

    if not overwrite and structure_output_path.exists():
        raise FileExistsError(
            f"Output structure already exists: {structure_output_path}. Use --overwrite."
        )
    if not overwrite and annotation_output_path.exists():
        raise FileExistsError(
            f"Output annotation already exists: {annotation_output_path}. Use --overwrite."
        )

    input_structure = pdb.load_structure(input_path)
    prepared_structure = volumizer.prepare_pdb_structure(input_structure)
    annotation_df, annotation_structure = volumizer.annotate_structure_volumes(
        prepared_structure,
        min_voxels=min_voxels,
        min_volume=min_volume,
    )

    combined_structure = prepared_structure + annotation_structure
    pdb.save_structure(combined_structure, structure_output_path)

    payload = _build_annotation_payload(
        source_label,
        input_path,
        structure_output_path,
        annotation_df,
    )
    annotation_output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return {
        "source": source_label,
        "input_path": str(input_path),
        "structure_output": str(structure_output_path),
        "annotation_output": str(annotation_output_path),
        "num_volumes": payload["num_volumes"],
        "largest_type": payload["largest_type"],
        "largest_volume": payload["largest_volume"],
    }


def _plan_dry_run(
    structures: list[tuple[str, Path]],
    args: argparse.Namespace,
    output_dir: Path,
    tracker: _RunTracker,
) -> None:
    planned_now = 0

    for source_label, input_path in structures:
        if args.resume and _has_complete_outputs(output_dir, source_label):
            if tracker.has_result(source_label):
                continue
            print(f"skipping {source_label}: outputs already exist (--resume)", file=sys.stderr)
            tracker.mark_skipped(
                _build_resume_skip_entry(
                    source_label=source_label,
                    input_path=input_path,
                    output_dir=output_dir,
                )
            )
            continue

        tracker.mark_planned(
            _build_dry_run_plan_entry(
                source_label=source_label,
                input_path=input_path,
                output_dir=output_dir,
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
) -> int:
    pending_structures: list[tuple[str, Path]] = []
    for source_label, input_path in structures:
        if args.resume and _has_complete_outputs(output_dir, source_label):
            if tracker.has_result(source_label):
                continue
            print(f"skipping {source_label}: outputs already exist (--resume)", file=sys.stderr)
            tracker.mark_skipped(
                _build_resume_skip_entry(
                    source_label=source_label,
                    input_path=input_path,
                    output_dir=output_dir,
                )
            )
            continue
        pending_structures.append((source_label, input_path))

    analysis_workers = int(args.jobs)
    if args.fail_fast and analysis_workers > 1:
        print(
            "--fail-fast enabled; forcing analysis worker count to 1.",
            file=sys.stderr,
        )
        analysis_workers = 1

    if len(pending_structures) == 0:
        return analysis_workers

    if analysis_workers <= 1:
        for source_label, input_path in pending_structures:
            tracker.mark_started(source_label, input_path)
            print(f"analyzing {source_label}: {input_path}", file=sys.stderr)
            try:
                tracker.mark_result(
                    analyze_structure_file(
                        source_label=source_label,
                        input_path=input_path,
                        output_dir=output_dir,
                        min_voxels=int(args.min_voxels),
                        min_volume=args.min_volume,
                        overwrite=bool(args.overwrite),
                    )
                )
            except Exception as error:  # pragma: no cover - exercised via CLI integration tests
                error_entry = {
                    "source": source_label,
                    "input_path": str(input_path),
                    "error": str(error),
                }
                tracker.mark_error(error_entry)
                print(f"error for {source_label}: {error}", file=sys.stderr)
                if args.fail_fast:
                    break

        return analysis_workers

    print(
        f"analyzing {len(pending_structures)} structures with {analysis_workers} workers...",
        file=sys.stderr,
    )
    with ThreadPoolExecutor(max_workers=analysis_workers) as executor:
        future_to_context = {
            executor.submit(
                analyze_structure_file,
                source_label=source_label,
                input_path=input_path,
                output_dir=output_dir,
                min_voxels=int(args.min_voxels),
                min_volume=args.min_volume,
                overwrite=bool(args.overwrite),
            ): (index, source_label, input_path)
            for index, (source_label, input_path) in enumerate(pending_structures)
        }

        for _, source_label, input_path in future_to_context.values():
            tracker.mark_started(source_label, input_path)

        ordered_results: dict[int, dict] = {}
        ordered_errors: dict[int, dict] = {}
        completed = 0

        for future in as_completed(future_to_context):
            index, source_label, input_path = future_to_context[future]
            completed += 1
            try:
                ordered_results[index] = future.result()
                print(
                    f"[{completed}/{len(pending_structures)}] completed {source_label}",
                    file=sys.stderr,
                )
            except Exception as error:  # pragma: no cover - exercised via CLI integration tests
                ordered_errors[index] = {
                    "source": source_label,
                    "input_path": str(input_path),
                    "error": str(error),
                }
                print(
                    f"[{completed}/{len(pending_structures)}] error for {source_label}: {error}",
                    file=sys.stderr,
                )

    for index in range(len(pending_structures)):
        if index in ordered_results:
            tracker.mark_result(ordered_results[index])
        elif index in ordered_errors:
            tracker.mark_error(ordered_errors[index])

    return analysis_workers


def run_cli(args: argparse.Namespace) -> int:
    """
    Execute CLI command and return process status code.
    """
    if args.jobs < 1:
        raise ValueError("--jobs must be >= 1.")
    if args.retries < 0:
        raise ValueError("--retries must be >= 0.")
    if args.retry_delay < 0:
        raise ValueError("--retry-delay must be >= 0.")

    if args.cluster_identity is None:
        if args.cluster_method is not None:
            raise ValueError("--cluster-method requires --cluster-identity.")
        if args.cluster_allow_all_methods:
            raise ValueError("--cluster-allow-all-methods requires --cluster-identity.")
        if args.cluster_max_resolution is not None:
            raise ValueError("--cluster-max-resolution requires --cluster-identity.")
        if args.metadata_cache is not None:
            raise ValueError("--metadata-cache requires --cluster-identity.")
        if args.no_metadata_cache:
            raise ValueError("--no-metadata-cache requires --cluster-identity.")

    if args.backend is not None:
        os.environ[native_backend.BACKEND_ENV] = args.backend
        native_backend.clear_backend_cache()

    utils.set_resolution(float(args.resolution))
    utils.set_non_protein(bool(args.keep_non_protein))

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    download_dir = (
        Path(args.download_dir)
        if args.download_dir is not None
        else output_dir / "downloads"
    )

    cluster_method_filters = None
    metadata_cache_path = None
    if args.cluster_identity is not None:
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

    structures = resolve_input_structures(args, download_dir, output_dir)
    tracker.mark_run_start(total_structures=len(structures), dry_run=args.dry_run)

    if args.dry_run:
        _plan_dry_run(structures, args, output_dir, tracker)
        analysis_workers = 0
    else:
        analysis_workers = _analyze_structures(structures, args, output_dir, tracker)

    summary = {
        "config": {
            "input": str(args.input) if args.input is not None else None,
            "pdb_id": args.pdb_id,
            "cluster_identity": args.cluster_identity,
            "max_structures": args.max_structures,
            "cluster_method_filters": cluster_method_filters,
            "cluster_allow_all_methods": args.cluster_allow_all_methods,
            "cluster_max_resolution": args.cluster_max_resolution,
            "metadata_cache": str(metadata_cache_path) if metadata_cache_path else None,
            "no_metadata_cache": args.no_metadata_cache,
            "checkpoint": str(checkpoint_path) if checkpoint_path else None,
            "no_checkpoint": args.no_checkpoint,
            "progress_jsonl": str(progress_jsonl_path) if progress_jsonl_path else None,
            "output_dir": str(output_dir),
            "download_dir": str(download_dir),
            "resolution": args.resolution,
            "min_voxels": args.min_voxels,
            "min_volume": args.min_volume,
            "backend": args.backend if args.backend is not None else utils.get_active_backend(),
            "keep_non_protein": args.keep_non_protein,
            "jobs": args.jobs,
            "analysis_workers": analysis_workers,
            "timeout": args.timeout,
            "retries": args.retries,
            "retry_delay": args.retry_delay,
            "overwrite": args.overwrite,
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

    exit_code = 1 if len(tracker.errors) > 0 else 0
    tracker.mark_run_complete(exit_code)
    return exit_code


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run_cli(args)
    except ValueError as error:
        parser.error(str(error))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
