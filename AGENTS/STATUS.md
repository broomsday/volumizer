# Migration Status Snapshot (2026-02-20)

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
- fib sphere batch expansion
- neighbor voxel detection
- BFS component expansion
- exposed/buried solvent split
- buried-component classification
- Native parity/integration test scaffolding exists and is passing.
- Backend comparison benchmark report published in `AGENTS/PERFORMANCE_REPORT.md`.
- Native optimization pass completed (batch fib expansion + neighbor-kernel rewrite) with updated benchmarks published.
- Orchestration-overhead reduction pass completed in Python integration layers (`fib_sphere.py`, `voxel.py`) with refreshed backend benchmarks published in `AGENTS/PERFORMANCE_REPORT.md`.
- New orchestration benchmark artifacts published:
- `AGENTS/benchmark.python.orchestration.json`
- `AGENTS/benchmark.ctypes-c.orchestration.json`
- `AGENTS/benchmark.native.orchestration.json`
- Stage-level benchmark profiling added (`scripts/benchmark.py` stage timing capture + summary + CLI stage breakdown output).
- Native classification split profiling (kernel vs mapping) completed with artifacts and report sections.
- Native mapping overhead reduction pass completed (fast axial-length kernel + mapping-path reductions), with medium/large repeat benchmarks showing substantial runtime improvements and mapping near-zero.
- Native prepare-structure overhead reduction pass completed (equivalent thin-SVD principal-axis alignment in `volumizer/align.py` + alignment parity tests), with medium/large repeat benchmarks showing additional runtime improvements and prepare stage near-zero.
- Native CI workflow added for parity tests and lightweight performance regression checks.
- CLI subcommand interface added (`analyze`, `cluster`, `cache`) for local files, PDB-ID download, cluster-representative runs, and metadata-cache inspection/maintenance while keeping legacy flag-only invocation compatibility.
- CLI resume mode and initial cluster selection filters added (default X-ray/cryo-EM method filter with optional method and resolution filters).
- CLI parallel worker support added (`--jobs`) with network retry controls and entry-metadata caching for large cluster runs.
- CLI dry-run mode added (`--dry-run`) plus negative metadata caching for permanent metadata failures (e.g. 404) to avoid repeat lookups.
- CLI checkpoint persistence added (default `run.checkpoint.json`) with resumable run-state recovery and optional JSONL progress events (`--progress-jsonl`).
- CLI manifest workflow added: cluster runs can write selected/rejected structure manifests (`--write-manifest`) and analyze runs can replay exact structure sets (`--manifest`).
- CLI summary replay added: analyze runs can replay prior run subsets via `--from-summary <run.summary.json> --only failed|skipped|planned|all` (default `failed`).
- CLI periodic progress/ETA output added via `--progress-interval` (default 30 seconds, disable with `<= 0`).
- CLI cluster sharding added via `--num-shards` + `--shard-index` for deterministic distributed processing of large representative lists.
- CLI failure export added via `--failures-manifest`, producing replayable analyze-manifest files from structure-level errors.

## Remaining

- Continue optimization focus on native classifier kernel cost, with stage-profile benchmark deltas and follow-on optimization of newly visible secondary stages.
- Finalize native packaging/wheel strategy and default backend policy.
- Continue CLI UX hardening: add validation subcommands and stronger runtime safety controls for very large jobs.
