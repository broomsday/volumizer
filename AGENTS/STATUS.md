# Migration Status Snapshot (2026-02-15)

## Completed

- `uv` migration is in place (`pyproject.toml`, README workflow).
- Python compatibility pinned (`>=3.10,<3.12`) with `numpy<2.0`.
- C helper build script works from repo root.
- Phase 0 baseline/guardrails artifacts are in place:
- benchmark harness
- golden parity scaffolding
- baseline metrics document
- Rust native extension scaffold exists and builds via `maturin`.
- Native backend selection and fallback wiring are implemented.
- Native kernels implemented and integrated:
- fib sphere points
- neighbor voxel detection
- BFS component expansion
- buried-component classification
- Native parity/integration test scaffolding exists and is passing.

## Remaining

- Port exposed/buried solvent split to native.
- Reduce Python loop/dataframe overhead around native kernels where possible.
- Publish backend-vs-backend benchmark report in `AGENTS/PERFORMANCE_REPORT.md`.
- Add CI native build + parity test job.
- Finalize native packaging/wheel strategy and default backend policy.
- Plan and implement CLI/UI improvements for broader usability.
