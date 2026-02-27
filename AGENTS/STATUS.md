# Migration Status Snapshot (2026-02-27)

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
- Native classifier kernel sub-stage profiling added (Rust-side timing + Python stage integration), with repeat benchmark artifacts identifying BFS and remaining-set rebuild as dominant kernel costs.
- Native classifier kernel optimization pass completed (lookup-based BFS + remaining-set elimination), with repeat benchmarks showing substantial medium/large runtime and kernel-time reductions.
- Native classify-type kernel optimization pass completed (hash-lookup adjacency + subset-BFS acceleration), with additional major medium/large runtime and kernel-time reductions.
- Native first-shell exposed optimization pass completed (specialized native first-shell kernel + Python integration path selection), with medium/large repeats showing ~25-27% first-shell stage reductions and ~2.8-3.0% end-to-end gains.
- Native volumes-to-structure optimization pass completed (vectorized atom-array assignment + single-allocation group assembly), with medium/large repeats showing ~90-93% `volumes_to_structure` stage reductions and ~11-17% end-to-end gains.
- Native add-extra-points optimization pass completed (factorized element grouping + preallocated output assembly in native integration path), with medium/large repeats showing ~16-17% `add_extra_points` stage reductions and up to ~5% end-to-end gains.
- Native exposed/buried split optimization pass completed (low-copy selection API + precomputed flat-index output), with medium/large repeats showing ~30-32% `get_exposed_and_buried_voxels` stage reductions and ~2.7-2.8% end-to-end gains.
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
- `load_structure` assembly-policy controls are implemented (`biological` default, plus `asymmetric` and `auto`) with strict identity-only shortcutting for CIF `auto` mode.
- `load_structure` sub-stage timing instrumentation is implemented (`load_structure_parse_decode`, `load_structure_assembly_expand`, `load_structure_fallback`) and wired into stage-level benchmark reporting.
- Assembly-policy controls are wired through library and tooling surfaces: `volumizer.volumize_pdb()`, `volumizer.volumize_pdb_and_save()`, CLI (`--assembly-policy`), checkpoint/signature metadata, and benchmark harness (`scripts/benchmark.py`).
- Regression coverage added for load-policy behavior and timing-path semantics in `tests/test_pdb.py` and CLI test scaffolding updated for policy propagation.
- RCSB format-path audit completed:
- PDB-ID and cluster workflows are confirmed to download `.cif` inputs (not `.pdb`).
- BCIF entry availability at `models.rcsb.org` is validated for sampled representative IDs, enabling a safe format-evaluation track for `load_structure`.

## Remaining

- Primary optimization track: `load_structure` end-to-end time reduction with assembly-safe behavior preserved.
- Evaluate BCIF input-path adoption where compatible, with mmCIF fallback retained for assembly-expansion paths.
- Continue native classify-path kernel tuning as a parallel track.
- Finalize native packaging/wheel strategy and default backend policy.
- Continue CLI UX hardening: add validation subcommands and stronger runtime safety controls for very large jobs.
