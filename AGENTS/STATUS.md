# Migration Status Snapshot (2026-02-18)

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

- Reduce Python loop/dataframe overhead around native kernels where possible.
- Finalize native packaging/wheel strategy and default backend policy.
- Continue CLI UX hardening: richer progress/ETA output semantics and additional selection controls beyond the current baseline.
