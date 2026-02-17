"""
Command-line interface for volumizer.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import sys
from typing import Sequence

from volumizer import native_backend, pdb, rcsb, utils, volumizer


DEFAULT_METADATA_CACHE_FILENAME = "entry_metadata_cache.json"


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


def _resolve_cluster_method_filters(args: argparse.Namespace) -> list[str] | None:
    if args.cluster_allow_all_methods:
        return None

    if args.cluster_method is None or len(args.cluster_method) == 0:
        return list(rcsb.DEFAULT_CLUSTER_METHOD_FILTERS)

    normalized_methods = []
    for raw_method in args.cluster_method:
        normalized_methods.append(rcsb.normalize_method_filter_name(raw_method))

    # De-duplicate while preserving user order.
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


def _load_metadata_cache(cache_path: Path | None) -> dict[str, dict]:
    if cache_path is None or not cache_path.is_file():
        return {}

    try:
        raw_payload = json.loads(cache_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}

    if not isinstance(raw_payload, dict):
        return {}

    cache: dict[str, dict] = {}
    for key, value in raw_payload.items():
        if not isinstance(key, str) or not isinstance(value, dict):
            continue
        try:
            normalized_key = rcsb.normalize_pdb_id(key)
        except ValueError:
            continue
        cache[normalized_key] = value

    return cache


def _save_metadata_cache(cache_path: Path, cache_payload: dict[str, dict]) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(
        json.dumps(dict(sorted(cache_payload.items())), indent=2),
        encoding="utf-8",
    )


def _fetch_metadata_for_ids(
    representative_ids: list[str],
    args: argparse.Namespace,
) -> tuple[dict[str, dict], int]:
    fetched: dict[str, dict] = {}
    errors = 0

    if len(representative_ids) == 0:
        return fetched, errors

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
                print(
                    f"metadata error for {representative_id}: {error}",
                    file=sys.stderr,
                )
                if args.fail_fast:
                    raise
        return fetched, errors

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
                print(
                    f"metadata error for {representative_id}: {error}",
                    file=sys.stderr,
                )
                if args.fail_fast:
                    raise

    return fetched, errors


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
        metadata_cache = _load_metadata_cache(metadata_cache_path)

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
        cache_misses = 0
        cache_updates = 0

        batch_size = max(1, int(args.jobs))
        for batch_start in range(0, len(representative_ids), batch_size):
            if len(selected_ids) >= target_selection:
                break

            batch_ids = representative_ids[batch_start : batch_start + batch_size]
            batch_metadata: dict[str, dict] = {}
            missing_ids: list[str] = []

            for representative_id in batch_ids:
                cached_metadata = metadata_cache.get(representative_id)
                if isinstance(cached_metadata, dict):
                    batch_metadata[representative_id] = cached_metadata
                    cache_hits += 1
                else:
                    cache_misses += 1
                    missing_ids.append(representative_id)

            fetched_metadata, batch_errors = _fetch_metadata_for_ids(missing_ids, args)
            metadata_errors += batch_errors
            for representative_id, metadata in fetched_metadata.items():
                batch_metadata[representative_id] = metadata
                metadata_cache[representative_id] = metadata
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

        if metadata_cache_path is not None and cache_updates > 0:
            _save_metadata_cache(metadata_cache_path, metadata_cache)

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
            f"cache_misses={cache_misses}, "
            f"cache_updates={cache_updates}, "
            f"filtered_by_method={filtered_by_method}, "
            f"filtered_by_resolution={filtered_by_resolution}, "
            f"missing_resolution={filtered_missing_resolution}, "
            f"metadata_errors={metadata_errors}",
            file=sys.stderr,
        )

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


def _analyze_structures(
    structures: list[tuple[str, Path]],
    args: argparse.Namespace,
    output_dir: Path,
) -> tuple[list[dict], list[dict], list[dict], int]:
    results: list[dict] = []
    errors: list[dict] = []
    skipped: list[dict] = []

    pending_structures: list[tuple[str, Path]] = []
    for source_label, input_path in structures:
        if args.resume and _has_complete_outputs(output_dir, source_label):
            print(f"skipping {source_label}: outputs already exist (--resume)", file=sys.stderr)
            skipped.append(
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
        return results, errors, skipped, analysis_workers

    if analysis_workers <= 1:
        for source_label, input_path in pending_structures:
            print(f"analyzing {source_label}: {input_path}", file=sys.stderr)
            try:
                result = analyze_structure_file(
                    source_label=source_label,
                    input_path=input_path,
                    output_dir=output_dir,
                    min_voxels=int(args.min_voxels),
                    min_volume=args.min_volume,
                    overwrite=bool(args.overwrite),
                )
                results.append(result)
            except Exception as error:  # pragma: no cover - exercised via CLI integration tests
                error_entry = {
                    "source": source_label,
                    "input_path": str(input_path),
                    "error": str(error),
                }
                errors.append(error_entry)
                print(f"error for {source_label}: {error}", file=sys.stderr)
                if args.fail_fast:
                    break

        return results, errors, skipped, analysis_workers

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
            results.append(ordered_results[index])
        elif index in ordered_errors:
            errors.append(ordered_errors[index])

    return results, errors, skipped, analysis_workers


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

    structures = resolve_input_structures(args, download_dir, output_dir)

    results, errors, skipped, analysis_workers = _analyze_structures(
        structures,
        args,
        output_dir,
    )

    cluster_method_filters = None
    metadata_cache_path = None
    if args.cluster_identity is not None:
        cluster_method_filters = _resolve_cluster_method_filters(args)
        metadata_cache_path = _resolve_metadata_cache_path(args, output_dir)

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
        },
        "num_processed": len(results),
        "num_failed": len(errors),
        "num_skipped": len(skipped),
        "results": results,
        "errors": errors,
        "skipped": skipped,
    }

    summary_path = output_dir / "run.summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote summary: {summary_path}", file=sys.stderr)

    if len(errors) > 0:
        return 1
    return 0


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
