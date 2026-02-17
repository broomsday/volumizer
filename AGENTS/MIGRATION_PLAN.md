# Volumizer Migration Plan (Python -> Rust/C Core)

## 0. Current Status (2026-02-16)

Phase status:
- `Phase 0` (baseline + guardrails): completed.
- `Phase 1` (native scaffold + FFI contract): completed.
- `Phase 2` (fib sphere native path): completed for kernel + integration, including first batch-optimization pass.
- `Phase 3` (neighbor/BFS native kernels): completed for kernel + integration.
- `Phase 4` (native classification): partially completed.

Completed in Phase 4:
- Native `classify_buried_components` is implemented in Rust.
- Native exposed/buried solvent split kernel is implemented in Rust.
- Python dispatch can consume native classifier output.
- End-to-end parity scaffolding for native vs Python pipeline is in place.

Remaining in Phase 4:
- Remove remaining Python orchestration loops from classification path where feasible.
- Benchmark native backend thoroughly and tune hotspots.

## 1. Objective

Migrate the computationally expensive parts of `volumizer` from Python into a native core (Rust preferred, C acceptable) while preserving current scientific behavior and improving install/use ergonomics.

## 2. Success Criteria

Functional parity targets:
- Existing tests pass with the native path enabled.
- Volume class assignment parity on current fixtures (`cavity`, `pocket`, `pore`, `hub`) remains unchanged.
- Numerical outputs stay within tolerance:
- `volume`: exact match for current fixtures.
- `x/y/z` axial lengths: absolute delta <= `0.01`.

Performance targets:
- End-to-end runtime speedup >= `3x` on medium proteins at `2.0 A` voxel size.
- Runtime speedup >= `5x` for large structures where BFS/component labeling dominates.
- Peak memory increase <= `1.5x` relative to current Python path.

Usability targets:
- One-command install path with native acceleration available by default on common platforms.
- First-class CLI for non-programmatic usage.

## 3. Recommended Technical Direction

Primary recommendation:
- Use Rust for the new core with `pyo3` + `maturin` to expose a Python extension module.

Why this path:
- Better memory safety than manual C for graph/grid traversal code.
- Strong performance for dense loops and contiguous array processing.
- Better packaging story for prebuilt wheels than manual `.so` compilation scripts.

Compatibility strategy:
- Keep the current Python API in `volumizer/volumizer.py`.
- Route hot paths to native extension when available.
- Keep a pure-Python fallback for portability and debugging.

## 4. Scope by Layer

Remain Python initially:
- Structure loading/cleaning (`volumizer/pdb.py`).
- Output dataframe and PDB annotation formatting.
- Public API and orchestration layer.

Move to native first:
- Atom shell expansion (`fib_sphere` generation + batch expansion).
- Neighbor detection between voxel sets.
- Connected-component BFS over buried voxels.
- Surface agglomeration counting used for cavity/pocket/pore/hub classification.

Potential second wave native work:
- Replace/avoid PyntCloud dependency with native voxelization from float coordinates.

## 5. Phased Plan

## Phase 0: Baseline and Guardrails

Deliverables:
- Benchmark harness for representative inputs:
- `small`: fixture proteins (`tests/pdbs/cavity.pdb`, etc.).
- `medium`: `tests/pdbs/4jpn.pdb`.
- `large`: at least one larger assembly fixture added to repo.
- Golden outputs captured as JSON for regression checks.

Tasks:
- Add benchmark script (`scripts/benchmark.py`) with wall time and peak memory logging.
- Add parity test helper that compares old/new outputs with explicit tolerances.
- Record baseline numbers in `AGENTS/BASELINE_METRICS.md`.

Exit criteria:
- Baseline metrics and goldens committed.
- CI job runs parity checks in addition to existing tests.

## Phase 1: Native Scaffolding and FFI Contract

Deliverables:
- New native package directory (example: `native/`).
- Extension module importable from Python (example: `volumizer_native`).
- Stable array-based function signatures documented in `AGENTS/FFI_CONTRACT.md`.

Tasks:
- Define contiguous-array interfaces:
- `float32[N,3]` for coordinates.
- `int32[M,3]` for voxel indices.
- Return buffers plus counts instead of per-voxel callbacks.
- Add runtime feature flag:
- `utils.NATIVE_BACKEND = "python" | "rust"` with auto-detect default.

Exit criteria:
- Extension builds locally and in CI.
- Python can call a trivial native function with numpy arrays.

## Phase 2: Port Fib Sphere Expansion

Deliverables:
- Native implementation of:
- sample estimation.
- Fibonacci sphere coordinate generation.
- bulk atom expansion by element/radius.

Tasks:
- Replace row-wise pandas iteration with batch numpy -> native call -> numpy return.
- Keep existing `add_extra_points()` public behavior and return schema unchanged.
- Add direct parity tests between Python and native outputs.

Exit criteria:
- `tests/test_fib_sphere.py` parity path passes.
- Benchmark shows clear speedup for expansion stage.

## Phase 3: Port Neighbor and BFS Kernels

Deliverables:
- Native implementations of:
- `get_neighbor_voxels`.
- `breadth_first_search`.

Tasks:
- Keep index semantics identical to current `voxel.py`.
- Add randomized property tests for graph connectivity equivalence.
- Integrate into `voxel.get_first_shell_exposed_voxels()` and component discovery path.

Exit criteria:
- `tests/test_voxel.py` passes under native path.
- BFS-heavy benchmarks show >= `3x` speedup for those stages.

## Phase 4: Port Full Volume Classification Kernel

Deliverables:
- Native function that accepts voxelized protein/solvent data and returns classified components:
- hub, pore, pocket, cavity, occluded.
- component voxel indices.
- surface-contact indices.
- centers and axial lengths (or raw indices if Python computes metrics initially).

Tasks:
- Move `get_exposed_and_buried_voxels()`, agglomeration, and type assignment into native code.
- Reduce Python looping to orchestration and final formatting only.
- Add deterministic ordering rules for component IDs to preserve stable outputs.

Exit criteria:
- End-to-end `volumize_pdb` parity on all fixtures.
- Overall medium-protein runtime >= `3x` improvement.

## Phase 5: Packaging, Install, and UX

Deliverables:
- Cross-platform wheel builds (Linux/macOS; Windows if feasible).
- Removal of manual `src/compile_c_libs.sh` requirement for standard installs.
- CLI command (example: `volumizer run ...`) with JSON and annotated PDB outputs.

Tasks:
- Add CLI entrypoint with:
- input path, output paths, resolution, min volume/voxels, keep-non-protein toggle.
- progress and timing output.
- Add install docs for:
- pip install with native wheels.
- source build fallback.

Exit criteria:
- New user can run a full analysis via CLI after a clean install without manual compilation.

## Phase 6: Stabilization and Release

Deliverables:
- Release candidate with migration notes and performance report.
- Deprecated-path policy for old ctypes backend.

Tasks:
- Run expanded regression suite.
- Add failure telemetry logs around backend selection and fallback reasons.
- Publish `AGENTS/PERFORMANCE_REPORT.md` with before/after numbers.

Exit criteria:
- Tagged release with native backend as default.
- Documented fallback path for unsupported environments.

## 6. Parallel Workstreams

Workstream A (Core algorithms):
- Owns native kernels and parity with classification semantics.

Workstream B (Python integration):
- Owns wrapper interfaces, fallback behavior, dataframe/PDB outputs.

Workstream C (Install + UX):
- Owns wheels, CLI, docs, release automation.

Recommended sequencing:
- Start A + B immediately.
- Start C at end of Phase 2, complete in Phases 5-6.

## 7. Risk Register and Mitigations

Risk:
- Silent scientific drift from algorithm refactors.
Mitigation:
- Golden fixtures, strict parity tests, randomized property checks.

Risk:
- FFI overhead removes performance gains.
Mitigation:
- Batch interfaces only, zero-copy numpy views where possible, avoid per-item calls.

Risk:
- Packaging complexity across platforms.
Mitigation:
- `maturin` wheels in CI and explicit fallback to Python backend.

Risk:
- ID/order instability breaks downstream tooling.
Mitigation:
- Deterministic sort keys for components and explicit test coverage for ordering.

## 8. Immediate Next Actions (Current)

1. Validate and tune the new native CI workflow (`.github/workflows/native-ci.yml`) across push/PR runs.
2. Define packaging plan for native wheels and default-backend policy (`python` vs `auto`).
3. Add release-oriented docs for native install/fallback behavior.
