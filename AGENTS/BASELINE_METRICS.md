# Volumizer Baseline Metrics

Baseline performance and memory data for all backends (python, ctypes-c, native) is documented in `AGENTS/PERFORMANCE_REPORT.md`.

## Benchmark Command

Run from repository root:

```bash
bash src/compile_c_libs.sh
uv run --python 3.11 python scripts/benchmark.py --group all --repeats 3 --resolution 2.0 --output-json benchmark.json
```

Optional focused run:

```bash
uv run --python 3.11 python scripts/benchmark.py --group medium --repeats 5 --resolution 2.0 --output-json benchmark.medium.json
```

## Notes

- Keep resolution constant (`2.0 A`) when comparing before/after changes.
- If backend changes, record separate baseline blocks by backend.
- If hardware changes, do not compare absolute timing directly; compare only within the same host profile.
