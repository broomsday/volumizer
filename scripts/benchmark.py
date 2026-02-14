#!/usr/bin/env python3
"""
Benchmark the volumizer pipeline on representative structures.

Usage examples:
  python scripts/benchmark.py
  python scripts/benchmark.py --group medium --repeats 3 --output-json AGENTS/benchmark.latest.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

@dataclass(frozen=True)
class BenchmarkCase:
    category: str
    name: str
    path: Path


DEFAULT_CASES = [
    BenchmarkCase("small", "cavity", ROOT_DIR / "tests" / "pdbs" / "cavity.pdb"),
    BenchmarkCase("small", "pocket", ROOT_DIR / "tests" / "pdbs" / "pocket.pdb"),
    BenchmarkCase("small", "pore", ROOT_DIR / "tests" / "pdbs" / "pore.pdb"),
    BenchmarkCase("small", "hub", ROOT_DIR / "tests" / "pdbs" / "hub.pdb"),
    BenchmarkCase("medium", "4jpn", ROOT_DIR / "tests" / "pdbs" / "4jpn.pdb"),
    BenchmarkCase("large", "4jpp_assembly", ROOT_DIR / "tests" / "pdbs" / "4jpp.cif"),
]


def max_rss_kb() -> float | None:
    """
    Return the process max RSS in kilobytes, if available.
    """
    try:
        import resource

        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform == "darwin":
            return float(rss) / 1024.0
        return float(rss)
    except Exception:
        return None


def import_runtime_modules():
    """
    Import volumizer runtime modules with a clear error on missing deps.
    """
    try:
        from volumizer import utils, volumizer

        return utils, volumizer
    except ModuleNotFoundError as error:
        missing = error.name if error.name is not None else str(error)
        raise RuntimeError(
            "Missing dependency while importing volumizer runtime: "
            f"{missing}. Install project dependencies before benchmarking."
        ) from error


def run_single_case(case: BenchmarkCase, resolution: float) -> dict[str, Any]:
    """
    Run one volumizer benchmark case and return result metadata.
    """
    utils, volumizer = import_runtime_modules()
    utils.set_resolution(resolution)

    start = time.perf_counter()
    annotation_df, _, _ = volumizer.volumize_pdb(case.path)
    elapsed_seconds = time.perf_counter() - start

    if utils.using_native():
        backend = "native"
    elif utils.using_performant():
        backend = "ctypes-c"
    else:
        backend = "python"

    largest_type = None
    largest_volume = None
    if not annotation_df.empty:
        largest_type = str(annotation_df.iloc[0]["type"])
        largest_volume = float(annotation_df.iloc[0]["volume"])

    return {
        "category": case.category,
        "case_name": case.name,
        "path": str(case.path.relative_to(ROOT_DIR)),
        "resolution": resolution,
        "backend": backend,
        "elapsed_seconds": round(elapsed_seconds, 6),
        "max_rss_kb": max_rss_kb(),
        "num_detected_volumes": int(len(annotation_df)),
        "largest_type": largest_type,
        "largest_volume": largest_volume,
        "timestamp_utc": datetime.now(tz=timezone.utc).isoformat(),
    }


def summarize(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """
    Aggregate per-run results by case.
    """
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for row in results:
        key = (row["category"], row["case_name"], row["path"])
        grouped.setdefault(key, []).append(row)

    summary_rows = []
    for key in sorted(grouped.keys()):
        category, case_name, path = key
        runs = grouped[key]
        elapsed_values = [float(r["elapsed_seconds"]) for r in runs]
        rss_values = [float(r["max_rss_kb"]) for r in runs if r["max_rss_kb"] is not None]
        last = runs[-1]
        summary_rows.append(
            {
                "category": category,
                "case_name": case_name,
                "path": path,
                "runs": len(runs),
                "mean_seconds": round(sum(elapsed_values) / len(elapsed_values), 6),
                "min_seconds": round(min(elapsed_values), 6),
                "max_seconds": round(max(elapsed_values), 6),
                "max_rss_mb": round((max(rss_values) / 1024.0), 3) if rss_values else None,
                "backend": last["backend"],
                "largest_type": last["largest_type"],
                "largest_volume": last["largest_volume"],
                "num_detected_volumes": last["num_detected_volumes"],
            }
        )

    return summary_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark volumizer runtime and memory.")
    parser.add_argument(
        "--group",
        choices=["all", "small", "medium", "large"],
        default="all",
        help="Which benchmark group to run.",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=1,
        help="Number of runs per case.",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=2.0,
        help="Voxel resolution in Angstroms.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Optional path to write raw and summarized benchmark results as JSON.",
    )

    # Internal mode for per-case subprocess execution.
    parser.add_argument("--single-case", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--category", type=str, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--case-name", type=str, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--case-path", type=Path, default=None, help=argparse.SUPPRESS)

    return parser.parse_args()


def print_summary(summary_rows: list[dict[str, Any]]) -> None:
    """
    Emit a compact text table to stdout.
    """
    print(
        "category case runs mean_s min_s max_s max_rss_mb backend largest_type largest_volume num_volumes"
    )
    for row in summary_rows:
        largest_volume = (
            "None" if row["largest_volume"] is None else f"{row['largest_volume']:.3f}"
        )
        max_rss_mb = "None" if row["max_rss_mb"] is None else f"{row['max_rss_mb']:.3f}"
        print(
            f"{row['category']} {row['case_name']} {row['runs']} "
            f"{row['mean_seconds']:.6f} {row['min_seconds']:.6f} {row['max_seconds']:.6f} "
            f"{max_rss_mb} {row['backend']} {row['largest_type']} {largest_volume} "
            f"{row['num_detected_volumes']}"
        )


def main() -> int:
    args = parse_args()

    if args.repeats < 1:
        print("--repeats must be >= 1", file=sys.stderr)
        return 2

    if args.single_case:
        if args.category is None or args.case_name is None or args.case_path is None:
            print(
                "--single-case requires --category, --case-name, and --case-path",
                file=sys.stderr,
            )
            return 2

        case = BenchmarkCase(
            category=args.category, name=args.case_name, path=args.case_path
        )
        try:
            result = run_single_case(case, args.resolution)
        except RuntimeError as error:
            print(str(error), file=sys.stderr)
            return 1
        print(json.dumps(result))
        return 0

    selected_groups = {"small", "medium", "large"}
    if args.group != "all":
        selected_groups = {args.group}

    cases = [case for case in DEFAULT_CASES if case.category in selected_groups]
    existing_cases = [case for case in cases if case.path.is_file()]
    skipped_cases = [
        {
            "category": case.category,
            "name": case.name,
            "path": str(case.path.relative_to(ROOT_DIR)),
        }
        for case in cases
        if not case.path.is_file()
    ]

    if not existing_cases:
        print("No benchmark cases were found on disk.", file=sys.stderr)
        return 1

    all_results: list[dict[str, Any]] = []
    script_path = Path(__file__).resolve()
    for case in existing_cases:
        for repeat in range(args.repeats):
            cmd = [
                sys.executable,
                str(script_path),
                "--single-case",
                "--category",
                case.category,
                "--case-name",
                case.name,
                "--case-path",
                str(case.path),
                "--resolution",
                str(args.resolution),
            ]
            run = subprocess.run(
                cmd,
                cwd=ROOT_DIR,
                text=True,
                capture_output=True,
                check=False,
            )
            if run.returncode != 0:
                print(
                    f"Benchmark case failed for {case.name} repeat {repeat + 1}:",
                    file=sys.stderr,
                )
                print(run.stderr.strip(), file=sys.stderr)
                return run.returncode
            all_results.append(json.loads(run.stdout.strip()))

    summary_rows = summarize(all_results)
    print_summary(summary_rows)

    if args.output_json is not None:
        output = {
            "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
            "resolution": args.resolution,
            "repeats": args.repeats,
            "group": args.group,
            "results": all_results,
            "summary": summary_rows,
            "skipped_cases": skipped_cases,
        }
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(output, indent=2), encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
