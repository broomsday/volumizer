#!/usr/bin/env python3
"""
Delete corrupted downloaded CIF files from a volumizer run output directory
and clean stale checkpoint entries so the affected structures are re-downloaded
and re-analyzed on the next resume.

Corruption is detected by scanning each `<download-dir>/*.cif` for stray control
bytes (any byte < 0x20 other than tab/LF/CR). This matches the symptom of errors
like `ValueError: could not convert string to float: '405.\\x1792'` where a
truncated/garbled download left a non-printable character inside a data value.

The manifest produced by `volumizer cluster --write-manifest` lists structures
that must be downloaded; it has no per-entry download-status field. Once a
corrupt CIF is removed, the next run naturally re-fetches it (download is
skipped only when the local file already exists). Clearing matching entries from
the checkpoint's `results`/`errors`/`skipped`/`planned` sections ensures a
`--resume` rerun treats them as fresh work.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


CHECKPOINT_ENTRY_KINDS = ("results", "errors", "skipped", "planned")
DEFAULT_PROGRESS_INTERVAL_SECONDS = 2.0

_BAD_BYTE_PATTERN = re.compile(rb"[\x00-\x08\x0b\x0c\x0e-\x1f]")


def _default_jobs() -> int:
    return min(8, os.cpu_count() or 1)


def inspect_cif(path: Path) -> tuple[bool, str | None]:
    try:
        data = path.read_bytes()
    except OSError as error:
        return True, f"unreadable: {error}"

    if len(data) == 0:
        return True, "empty file"

    match = _BAD_BYTE_PATTERN.search(data)
    if match is not None:
        return True, f"contains control byte 0x{match.group()[0]:02X}"

    return False, None


def discover_downloads(download_dir: Path) -> list[Path]:
    if not download_dir.is_dir():
        raise FileNotFoundError(f"Download directory does not exist: {download_dir}")
    return sorted(download_dir.glob("*.cif"))


def source_label_for_cif(cif_path: Path) -> str:
    return cif_path.stem


def remove_sources_from_checkpoint(
    checkpoint_path: Path,
    sources_to_drop: set[str],
    dry_run: bool,
) -> dict[str, int]:
    removed_counts = dict.fromkeys(CHECKPOINT_ENTRY_KINDS, 0)

    if not checkpoint_path.is_file() or len(sources_to_drop) == 0:
        return removed_counts

    try:
        payload = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as error:
        print(
            f"warning: could not read checkpoint {checkpoint_path}: {error}",
            file=sys.stderr,
        )
        return removed_counts

    if not isinstance(payload, dict):
        return removed_counts

    changed = False
    for kind in CHECKPOINT_ENTRY_KINDS:
        entries = payload.get(kind)
        if not isinstance(entries, list):
            continue
        kept: list[dict] = []
        for entry in entries:
            if isinstance(entry, dict) and entry.get("source") in sources_to_drop:
                removed_counts[kind] += 1
                changed = True
                continue
            kept.append(entry)
        payload[kind] = kept

    if changed and not dry_run:
        checkpoint_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return removed_counts


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Delete corrupted CIFs from a volumizer run's downloads directory "
            "and clear matching checkpoint entries so they are re-fetched on "
            "the next --resume."
        ),
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Volumizer run output directory.",
    )
    parser.add_argument(
        "--download-dir",
        type=Path,
        default=None,
        help="Override downloads directory (default: <output-dir>/downloads).",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=None,
        help="Override checkpoint path (default: <output-dir>/run.checkpoint.json).",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=_default_jobs(),
        help=(
            "Parallel worker threads for scanning CIF files "
            f"(default: min(8, cpu_count) = {_default_jobs()})."
        ),
    )
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=DEFAULT_PROGRESS_INTERVAL_SECONDS,
        help=(
            "Seconds between progress updates during scanning "
            f"(default: {DEFAULT_PROGRESS_INTERVAL_SECONDS}). Set <=0 to disable."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report what would be removed without modifying files.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print every scanned file, not just corrupted ones.",
    )
    return parser.parse_args(argv)


def _emit_progress(
    processed: int,
    total: int,
    corrupt_count: int,
    start_time: float,
) -> None:
    elapsed = max(1e-9, time.monotonic() - start_time)
    rate = processed / elapsed
    percent = 100.0 * processed / total if total > 0 else 100.0
    remaining = max(0, total - processed)
    eta_seconds = remaining / rate if rate > 0 else 0.0
    print(
        f"scan progress: {processed}/{total} ({percent:.1f}%), "
        f"corrupt={corrupt_count}, rate={rate:.1f}/s, eta={eta_seconds:.0f}s",
        file=sys.stderr,
    )


def scan_downloads(
    cif_paths: list[Path],
    jobs: int,
    verbose: bool,
    progress_interval_seconds: float,
) -> list[tuple[Path, str]]:
    total = len(cif_paths)
    corrupt: list[tuple[Path, str]] = []
    start_time = time.monotonic()
    last_emit = start_time
    processed = 0
    workers = max(1, int(jobs))

    def _record(path: Path, is_bad: bool, reason: str | None) -> None:
        nonlocal last_emit
        if is_bad:
            corrupt.append((path, reason or "unknown"))
            print(f"corrupt: {path} ({reason})")
        elif verbose:
            print(f"ok: {path}")
        if progress_interval_seconds > 0:
            now = time.monotonic()
            if now - last_emit >= progress_interval_seconds:
                _emit_progress(processed, total, len(corrupt), start_time)
                last_emit = now

    if workers <= 1:
        for path in cif_paths:
            is_bad, reason = inspect_cif(path)
            processed += 1
            _record(path, is_bad, reason)
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_path = {
                executor.submit(inspect_cif, path): path for path in cif_paths
            }
            for future in as_completed(future_to_path):
                path = future_to_path[future]
                is_bad, reason = future.result()
                processed += 1
                _record(path, is_bad, reason)

    _emit_progress(processed, total, len(corrupt), start_time)
    return corrupt


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    output_dir: Path = args.output_dir
    if not output_dir.is_dir():
        print(f"error: output directory does not exist: {output_dir}", file=sys.stderr)
        return 2

    download_dir = args.download_dir if args.download_dir is not None else output_dir / "downloads"
    checkpoint_path = (
        args.checkpoint if args.checkpoint is not None else output_dir / "run.checkpoint.json"
    )

    try:
        cif_paths = discover_downloads(download_dir)
    except FileNotFoundError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    total = len(cif_paths)
    print(
        f"scanning {total} cif file(s) with {args.jobs} worker thread(s)...",
        file=sys.stderr,
    )
    corrupt_paths = scan_downloads(
        cif_paths,
        jobs=args.jobs,
        verbose=args.verbose,
        progress_interval_seconds=args.progress_interval,
    )

    action = "would delete" if args.dry_run else "deleting"
    for cif_path, _reason in corrupt_paths:
        print(f"{action}: {cif_path}")
        if not args.dry_run:
            try:
                cif_path.unlink()
            except OSError as error:
                print(f"warning: failed to delete {cif_path}: {error}", file=sys.stderr)

    sources_to_drop = {source_label_for_cif(path) for path, _ in corrupt_paths}
    removed_counts = remove_sources_from_checkpoint(
        checkpoint_path,
        sources_to_drop,
        dry_run=args.dry_run,
    )

    print(
        f"scanned {total} cif file(s); "
        f"corrupt={len(corrupt_paths)}; "
        f"checkpoint_removed="
        + ",".join(f"{kind}={removed_counts[kind]}" for kind in CHECKPOINT_ENTRY_KINDS)
        + (" (dry-run)" if args.dry_run else "")
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
