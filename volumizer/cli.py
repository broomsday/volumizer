"""
Command-line interface for volumizer.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
from typing import Sequence

from volumizer import native_backend, pdb, rcsb, utils, volumizer


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
        "--timeout",
        type=float,
        default=60.0,
        help="Network timeout in seconds for RCSB downloads (default: 60).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing downloaded/output files.",
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


def resolve_input_structures(
    args: argparse.Namespace,
    download_dir: Path,
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
        structure_path = rcsb.download_structure_cif(
            normalized_id,
            output_dir=download_dir,
            overwrite=args.overwrite,
            timeout=args.timeout,
        )
        return [(_sanitize_label(normalized_id), structure_path)]

    if args.cluster_identity is not None:
        representative_ids = rcsb.fetch_cluster_representative_entry_ids(
            identity=args.cluster_identity,
            max_structures=args.max_structures,
            timeout=args.timeout,
            include_non_pdb=False,
        )
        if len(representative_ids) == 0:
            raise RuntimeError(
                f"No representative PDB IDs found for cluster identity {args.cluster_identity}."
            )

        structures = []
        for index, representative_id in enumerate(representative_ids, start=1):
            print(
                f"[{index}/{len(representative_ids)}] downloading {representative_id}...",
                file=sys.stderr,
            )
            structure_path = rcsb.download_structure_cif(
                representative_id,
                output_dir=download_dir,
                overwrite=args.overwrite,
                timeout=args.timeout,
            )
            structures.append((_sanitize_label(representative_id), structure_path))

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
    structure_output_path = output_dir / f"{source_label}.annotated.cif"
    annotation_output_path = output_dir / f"{source_label}.annotation.json"

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


def run_cli(args: argparse.Namespace) -> int:
    """
    Execute CLI command and return process status code.
    """
    if args.backend is not None:
        os.environ[native_backend.BACKEND_ENV] = args.backend
        native_backend.clear_backend_cache()

    utils.set_resolution(float(args.resolution))
    utils.set_non_protein(bool(args.keep_non_protein))

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    download_dir = Path(args.download_dir) if args.download_dir is not None else output_dir / "downloads"

    structures = resolve_input_structures(args, download_dir)

    results = []
    errors = []

    for source_label, input_path in structures:
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

    summary = {
        "config": {
            "input": str(args.input) if args.input is not None else None,
            "pdb_id": args.pdb_id,
            "cluster_identity": args.cluster_identity,
            "max_structures": args.max_structures,
            "output_dir": str(output_dir),
            "download_dir": str(download_dir),
            "resolution": args.resolution,
            "min_voxels": args.min_voxels,
            "min_volume": args.min_volume,
            "backend": args.backend if args.backend is not None else utils.get_active_backend(),
            "keep_non_protein": args.keep_non_protein,
        },
        "num_processed": len(results),
        "num_failed": len(errors),
        "results": results,
        "errors": errors,
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
    return run_cli(args)


if __name__ == "__main__":
    raise SystemExit(main())
