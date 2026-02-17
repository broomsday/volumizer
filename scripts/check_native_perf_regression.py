#!/usr/bin/env python3
"""
Run a lightweight native performance regression check against a same-run reference backend.

By default this compares native against the current python/ctypes backend on the medium group.
It is designed to be robust across different machines by using relative runtime ratios.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
BENCHMARK_SCRIPT = ROOT_DIR / "scripts" / "benchmark.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check native runtime ratio vs reference backend."
    )
    parser.add_argument(
        "--group",
        choices=["small", "medium", "large"],
        default="medium",
        help="Benchmark group to run for the check.",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=1,
        help="Repeats passed to scripts/benchmark.py.",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=2.0,
        help="Voxel resolution in Angstroms.",
    )
    parser.add_argument(
        "--max-ratio",
        type=float,
        default=1.10,
        help=(
            "Maximum allowed native/reference mean runtime ratio. "
            "Values >1 mean native slower than reference."
        ),
    )
    parser.add_argument(
        "--reference-backend",
        choices=["python", "native"],
        default="python",
        help="Backend mode used as the reference run.",
    )
    return parser.parse_args()


def run_benchmark(backend: str, group: str, repeats: int, resolution: float) -> dict:
    with tempfile.NamedTemporaryFile(
        mode="w+", suffix=".json", prefix=f"benchmark_{backend}_", delete=False
    ) as tmp_file:
        output_path = Path(tmp_file.name)

    cmd = [
        sys.executable,
        str(BENCHMARK_SCRIPT),
        "--group",
        group,
        "--repeats",
        str(repeats),
        "--resolution",
        str(resolution),
        "--output-json",
        str(output_path),
    ]
    env = os.environ.copy()
    env["VOLUMIZER_BACKEND"] = backend

    result = subprocess.run(
        cmd,
        cwd=ROOT_DIR,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"benchmark run failed for backend={backend}")

    try:
        return json.loads(output_path.read_text(encoding="utf-8"))
    finally:
        output_path.unlink(missing_ok=True)


def summary_by_case(benchmark_json: dict) -> dict[str, dict]:
    return {row["case_name"]: row for row in benchmark_json["summary"]}


def main() -> int:
    args = parse_args()
    if args.max_ratio <= 0:
        print("--max-ratio must be > 0", file=sys.stderr)
        return 2

    print(
        f"Running regression check group={args.group} repeats={args.repeats} "
        f"resolution={args.resolution}..."
    )
    reference_json = run_benchmark(
        backend=args.reference_backend,
        group=args.group,
        repeats=args.repeats,
        resolution=args.resolution,
    )
    native_json = run_benchmark(
        backend="native",
        group=args.group,
        repeats=args.repeats,
        resolution=args.resolution,
    )

    reference_rows = summary_by_case(reference_json)
    native_rows = summary_by_case(native_json)
    case_names = sorted(set(reference_rows).intersection(native_rows))
    if len(case_names) == 0:
        print("No overlapping benchmark cases to compare.", file=sys.stderr)
        return 1

    failures = []
    print("case native_s reference_s ratio pass")
    for case_name in case_names:
        native_s = float(native_rows[case_name]["mean_seconds"])
        reference_s = float(reference_rows[case_name]["mean_seconds"])
        ratio = native_s / reference_s if reference_s > 0 else float("inf")
        is_pass = ratio <= args.max_ratio
        marker = "yes" if is_pass else "no"
        print(
            f"{case_name} {native_s:.6f} {reference_s:.6f} "
            f"{ratio:.3f} {marker}"
        )
        if not is_pass:
            failures.append(
                (
                    case_name,
                    native_s,
                    reference_s,
                    ratio,
                )
            )

    if failures:
        print(
            "\nNative performance regression detected "
            f"(max allowed ratio={args.max_ratio:.3f}).",
            file=sys.stderr,
        )
        for case_name, native_s, reference_s, ratio in failures:
            print(
                f"- {case_name}: native={native_s:.6f}s, "
                f"reference={reference_s:.6f}s, ratio={ratio:.3f}",
                file=sys.stderr,
            )
        return 1

    print(
        "\nNative performance check passed "
        f"(all ratios <= {args.max_ratio:.3f})."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
