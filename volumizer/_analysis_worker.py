"""
Internal worker entrypoint for isolated structure analysis.
"""


import argparse
import json
import os
import sys
import traceback
from pathlib import Path

from volumizer.cli import (
    PostAssemblyResidueLimitExceeded,
    analyze_structure_file,
)
from volumizer.pdb import VALID_ASSEMBLY_POLICIES
from volumizer import native_backend, utils


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m volumizer._analysis_worker",
        description="Internal isolated worker for volumizer CLI analysis.",
    )
    parser.add_argument("--source-label", required=True)
    parser.add_argument("--input-path", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--min-voxels", type=int, required=True)
    parser.add_argument("--min-volume", type=float, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--assembly-policy",
        choices=VALID_ASSEMBLY_POLICIES,
        default="biological",
    )
    parser.add_argument("--resolution", type=float, required=True)
    parser.add_argument("--keep-non-protein", action="store_true")
    parser.add_argument("--backend", default=None)
    parser.add_argument("--max-residues", type=int, default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        if args.backend:
            os.environ[native_backend.BACKEND_ENV] = args.backend
            native_backend.clear_backend_cache()

        utils.set_resolution(float(args.resolution))
        utils.set_non_protein(bool(args.keep_non_protein))

        result = analyze_structure_file(
            source_label=str(args.source_label),
            input_path=Path(args.input_path),
            output_dir=Path(args.output_dir),
            min_voxels=int(args.min_voxels),
            min_volume=args.min_volume,
            overwrite=bool(args.overwrite),
            assembly_policy=str(args.assembly_policy),
            max_residues=args.max_residues,
        )
    except PostAssemblyResidueLimitExceeded as error:
        json.dump(
            {
                "status": "skipped_post_assembly_residue_limit",
                "actual_residues": error.actual_residues,
                "max_residues": error.max_residues,
                "assembly_policy": error.assembly_policy,
            },
            sys.stdout,
        )
        sys.stdout.write("\n")
        return 0
    except Exception:  # pragma: no cover - exercised via CLI subprocess integration
        traceback.print_exc(file=sys.stderr)
        return 1

    json.dump(result, sys.stdout)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
