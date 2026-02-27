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

VALID_ASSEMBLY_POLICIES = ("biological", "asymmetric", "auto")


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


def _normalize_stage_timings(
    stage_timings: dict[str, float], elapsed_seconds: float
) -> tuple[dict[str, float], float]:
    """
    Round stage timings and include an explicit residual for uninstrumented time.
    """
    normalized = {
        stage: round(float(seconds), 6)
        for stage, seconds in sorted(stage_timings.items())
    }
    stage_sum = float(sum(normalized.values()))
    residual = max(0.0, float(elapsed_seconds) - stage_sum)
    normalized["untracked"] = round(residual, 6)
    return normalized, round(stage_sum, 6)


def run_single_case(
    case: BenchmarkCase,
    resolution: float,
    assembly_policy: str = "biological",
) -> dict[str, Any]:
    """
    Run one volumizer benchmark case and return result metadata.
    """
    utils, volumizer = import_runtime_modules()
    utils.set_resolution(resolution)

    stage_timings: dict[str, float] = {}
    start = time.perf_counter()
    annotation_df, _, _ = volumizer.volumize_pdb(
        case.path,
        stage_timings=stage_timings,
        assembly_policy=assembly_policy,
    )
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

    normalized_stage_timings, tracked_stage_seconds = _normalize_stage_timings(
        stage_timings,
        elapsed_seconds,
    )

    return {
        "category": case.category,
        "case_name": case.name,
        "path": str(case.path.relative_to(ROOT_DIR)),
        "resolution": resolution,
        "assembly_policy": assembly_policy,
        "backend": backend,
        "elapsed_seconds": round(elapsed_seconds, 6),
        "tracked_stage_seconds": tracked_stage_seconds,
        "stage_timings_seconds": normalized_stage_timings,
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
    grouped: dict[tuple[str, str, str, str], list[dict[str, Any]]] = {}
    for row in results:
        key = (
            row["category"],
            row["case_name"],
            row["path"],
            row.get("assembly_policy", "biological"),
        )
        grouped.setdefault(key, []).append(row)

    summary_rows = []
    for key in sorted(grouped.keys()):
        category, case_name, path, assembly_policy = key
        runs = grouped[key]
        elapsed_values = [float(r["elapsed_seconds"]) for r in runs]
        mean_seconds = sum(elapsed_values) / len(elapsed_values)
        rss_values = [float(r["max_rss_kb"]) for r in runs if r["max_rss_kb"] is not None]
        last = runs[-1]

        stage_names = sorted(
            {
                stage
                for run in runs
                for stage in run.get("stage_timings_seconds", {}).keys()
            }
        )
        stage_mean_seconds: dict[str, float] = {}
        for stage in stage_names:
            stage_mean_seconds[stage] = round(
                sum(float(run.get("stage_timings_seconds", {}).get(stage, 0.0)) for run in runs)
                / len(runs),
                6,
            )

        stage_mean_percent_of_total: dict[str, float] = {}
        if mean_seconds > 0:
            for stage, value in stage_mean_seconds.items():
                stage_mean_percent_of_total[stage] = round((value / mean_seconds) * 100.0, 3)

        summary_rows.append(
            {
                "category": category,
                "case_name": case_name,
                "path": path,
                "runs": len(runs),
                "mean_seconds": round(mean_seconds, 6),
                "min_seconds": round(min(elapsed_values), 6),
                "max_seconds": round(max(elapsed_values), 6),
                "max_rss_mb": round((max(rss_values) / 1024.0), 3) if rss_values else None,
                "backend": last["backend"],
                "largest_type": last["largest_type"],
                "largest_volume": last["largest_volume"],
                "num_detected_volumes": last["num_detected_volumes"],
                "assembly_policy": assembly_policy,
                "stage_mean_seconds": stage_mean_seconds,
                "stage_mean_percent_of_total": stage_mean_percent_of_total,
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
        "--assembly-policy",
        choices=VALID_ASSEMBLY_POLICIES,
        default="biological",
        help="Assembly policy passed to load_structure (default: biological).",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Optional path to write raw and summarized benchmark results as JSON.",
    )
    parser.add_argument(
        "--print-stage-breakdown",
        action="store_true",
        help="Print per-case top stage-time contributors.",
    )
    parser.add_argument(
        "--stage-top-n",
        type=int,
        default=5,
        help="How many stages to show per case when printing stage breakdown.",
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
        "category case assembly_policy runs mean_s min_s max_s max_rss_mb backend largest_type largest_volume num_volumes"
    )
    for row in summary_rows:
        largest_volume = (
            "None" if row["largest_volume"] is None else f"{row['largest_volume']:.3f}"
        )
        max_rss_mb = "None" if row["max_rss_mb"] is None else f"{row['max_rss_mb']:.3f}"
        print(
            f"{row['category']} {row['case_name']} {row['assembly_policy']} {row['runs']} "
            f"{row['mean_seconds']:.6f} {row['min_seconds']:.6f} {row['max_seconds']:.6f} "
            f"{max_rss_mb} {row['backend']} {row['largest_type']} {largest_volume} "
            f"{row['num_detected_volumes']}"
        )


def print_stage_breakdown(summary_rows: list[dict[str, Any]], stage_top_n: int) -> None:
    """
    Emit a per-case breakdown of top stage contributors.
    """
    safe_top_n = max(1, int(stage_top_n))
    print("\ncase stage mean_s pct_total")
    for row in summary_rows:
        stage_percent = row.get("stage_mean_percent_of_total", {})
        stage_seconds = row.get("stage_mean_seconds", {})
        ordered = sorted(
            stage_seconds.keys(),
            key=lambda stage_name: float(stage_seconds.get(stage_name, 0.0)),
            reverse=True,
        )[:safe_top_n]

        for stage_name in ordered:
            print(
                f"{row['case_name']} {stage_name} "
                f"{float(stage_seconds.get(stage_name, 0.0)):.6f} "
                f"{float(stage_percent.get(stage_name, 0.0)):.3f}"
            )


def main() -> int:
    args = parse_args()

    if args.repeats < 1:
        print("--repeats must be >= 1", file=sys.stderr)
        return 2

    if args.stage_top_n < 1:
        print("--stage-top-n must be >= 1", file=sys.stderr)
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
            result = run_single_case(
                case,
                args.resolution,
                assembly_policy=args.assembly_policy,
            )
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
                "--assembly-policy",
                args.assembly_policy,
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
    if args.print_stage_breakdown:
        print_stage_breakdown(summary_rows, args.stage_top_n)

    if args.output_json is not None:
        output = {
            "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
            "resolution": args.resolution,
            "assembly_policy": args.assembly_policy,
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
