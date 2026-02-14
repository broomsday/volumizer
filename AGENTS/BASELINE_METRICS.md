# Volumizer Baseline Metrics

This file tracks baseline performance and memory for the current Python-first implementation.
Use it as the parity/performance reference while migrating to a Rust/C core.

## 1. Benchmark Command

Run from repository root:

```bash
bash src/compile_c_libs.sh
uv run --python 3.11 python scripts/benchmark.py --group all --repeats 3 --resolution 2.0 --output-json AGENTS/benchmark.baseline.json
```

Optional focused run:

```bash
bash src/compile_c_libs.sh
uv run --python 3.11 python scripts/benchmark.py --group medium --repeats 5 --resolution 2.0 --output-json AGENTS/benchmark.medium.json
```

## 2. Environment Snapshot

- Date: `2026-02-13 17:05:19 UTC`
- Machine: local Linux workstation
- OS: `Fedora Linux 43 (Xfce)`, kernel `Linux 6.18.8-200.fc43.x86_64`
- CPU: `Intel(R) Core(TM) i9-10900X CPU @ 3.70GHz` (`20` logical CPUs)
- RAM: `131563892 kB` (~`125.47 GiB`)
- Python version: `3.11.14`
- `volumizer` backend detected (`python` or `ctypes-c`): `ctypes-c`

## 3. Baseline Summary

Copy from benchmark stdout or `summary` in `AGENTS/benchmark.baseline.json`.

| category | case | runs | mean_s | min_s | max_s | max_rss_mb | backend | largest_type | largest_volume | num_volumes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| small | cavity | 3 | 0.169252 | 0.166973 | 0.172421 | 177.223 | ctypes-c | cavity | 184.000 | 4 |
| small | pocket | 3 | 0.169022 | 0.166937 | 0.172782 | 175.332 | ctypes-c | pocket | 840.000 | 2 |
| small | pore | 3 | 0.171256 | 0.170047 | 0.173176 | 172.594 | ctypes-c | pore | 1904.000 | 2 |
| small | hub | 3 | 0.182690 | 0.179612 | 0.188292 | 171.918 | ctypes-c | hub | 2176.000 | 2 |
| medium | 4jpn | 3 | 8.463614 | 8.435856 | 8.513956 | 698.324 | ctypes-c | pore | 45048.000 | 40 |
| large | 4jpp_assembly | 3 | 23.775894 | 23.663393 | 23.832774 | 1714.801 | ctypes-c | pore | 71160.000 | 66 |

## 4. Notes

- Keep resolution constant (`2.0 A`) when comparing before/after migration.
- If backend changes, record separate baseline blocks by backend.
- If hardware changes, do not compare absolute timing directly; compare only within the same host profile.
