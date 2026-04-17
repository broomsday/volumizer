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
import sys
from pathlib import Path


ALLOWED_CONTROL_BYTES = frozenset({0x09, 0x0A, 0x0D})
CHECKPOINT_ENTRY_KINDS = ("results", "errors", "skipped", "planned")


def find_first_bad_byte(data: bytes) -> int | None:
    for byte in data:
        if byte < 0x20 and byte not in ALLOWED_CONTROL_BYTES:
            return byte
    return None


def inspect_cif(path: Path) -> tuple[bool, str | None]:
    try:
        data = path.read_bytes()
    except OSError as error:
        return True, f"unreadable: {error}"

    if len(data) == 0:
        return True, "empty file"

    bad_byte = find_first_bad_byte(data)
    if bad_byte is not None:
        return True, f"contains control byte 0x{bad_byte:02X}"

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
    corrupt_paths: list[tuple[Path, str]] = []

    for cif_path in cif_paths:
        is_bad, reason = inspect_cif(cif_path)
        if is_bad:
            corrupt_paths.append((cif_path, reason or "unknown"))
            print(f"corrupt: {cif_path} ({reason})")
        elif args.verbose:
            print(f"ok: {cif_path}")

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
