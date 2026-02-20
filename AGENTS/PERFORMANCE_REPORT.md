# Volumizer Backend Performance Report (2026-02-16)

## 1. Scope

- Benchmark harness: `scripts/benchmark.py`
- Resolution: `2.0 A`
- Repeats: `3`
- Group: `all` (`small`, `medium`, `large`)
- Backends compared:
- `python`: `VOLUMIZER_BACKEND=python` with `src/voxel.so` and `src/fib_sphere.so` temporarily moved out
- `ctypes-c`: `VOLUMIZER_BACKEND=python` with compiled `.so` accelerators present
- `native` (pre-optimization): `VOLUMIZER_BACKEND=native`
- `native` (optimized): `VOLUMIZER_BACKEND=native` after fib-sphere batching/caching and neighbor-kernel rewrite

Artifacts:
- `AGENTS/benchmark.python.json`
- `AGENTS/benchmark.ctypes-c.json`
- `AGENTS/benchmark.native.json`
- `AGENTS/benchmark.native.optimized.json`

## 2. Environment Snapshot

- Date generated (UTC): `2026-02-16 17:39:59 UTC`
- OS: `Fedora Linux 43 (Xfce)`
- Kernel: `6.18.8-200.fc43.x86_64`
- CPU: `Intel(R) Core(TM) i9-10900X CPU @ 3.70GHz` (`20` logical CPUs)
- RAM: `131563892 kB` (~`125.47 GiB`)
- Python: `3.11.14`

## 3. Current Runtime Comparison (With Optimized Native)

Mean runtime in seconds:

| case | python_s | ctypes_s | native_optimized_s | native_vs_python | ctypes_vs_python | native_vs_ctypes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cavity | 0.394818 | 0.169685 | 0.104879 | 3.765x faster | 2.327x faster | 1.618x |
| pocket | 0.539714 | 0.205055 | 0.099812 | 5.407x faster | 2.632x faster | 2.054x |
| pore | 0.760463 | 0.189863 | 0.097166 | 7.826x faster | 4.005x faster | 1.954x |
| hub | 0.717067 | 0.189071 | 0.104784 | 6.843x faster | 3.793x faster | 1.804x |
| 4jpn | 83.107042 | 8.548185 | 2.888686 | 28.770x faster | 9.722x faster | 2.959x |
| 4jpp_assembly | 207.212250 | 24.393836 | 9.264783 | 22.366x faster | 8.494x faster | 2.633x |

Interpretation of `native_vs_ctypes`:
- values `< 1.0` mean native is slower than ctypes-c.
- values `> 1.0` mean native is faster than ctypes-c.

## 4. Optimization Delta (Native Before vs After)

| case | native_before_s | native_after_s | improvement |
| --- | ---: | ---: | ---: |
| cavity | 0.196238 | 0.104879 | 46.56% |
| pocket | 0.217411 | 0.099812 | 54.09% |
| pore | 0.259198 | 0.097166 | 62.51% |
| hub | 0.260075 | 0.104784 | 59.71% |
| 4jpn | 20.007940 | 2.888686 | 85.56% |
| 4jpp_assembly | 57.477134 | 9.264783 | 83.88% |

Optimization implemented:
- Added native batch kernel for Fibonacci sphere generation (`fibonacci_sphere_points_batch`).
- Switched native `add_extra_points` path to batched center expansion when available, with fallback for older native builds.
- Added caching for repeated fib sample/radius estimation.
- Replaced native neighbor detection O(`N*M`) scan with hash-set based 6-neighbor lookup.

## 5. Memory Snapshot

Peak RSS remains similar before/after optimization:
- medium (`4jpn`): ~`698-699 MB`
- large (`4jpp_assembly`): ~`1713-1715 MB`

No meaningful memory regression observed from this optimization pass.

## 6. Output Consistency Check

Across compared backend runs, each case matched on:
- `num_detected_volumes`
- `largest_type`
- `largest_volume`

This indicates no obvious high-level classification drift from the optimization.

## 7. Success Criteria Check (Migration Plan)

Against defined native-vs-python targets:
- Medium (`4jpn`) target `>= 3x`: **met** (`28.770x`)
- Large (`4jpp_assembly`) target `>= 5x`: **met** (`22.366x`)

## 8. Current Conclusion

- Native backend is substantially faster than pure Python and now significantly faster than ctypes-c across benchmark cases.
- Largest gains came from optimizing native neighbor detection, which had dominated runtime on medium/large structures.
- Next optimization work can shift from raw speed to stabilization tasks (CI, packaging, and preventing regressions with stage-level perf checks).

## 9. Orchestration Overhead Reduction Delta (2026-02-19)

Scope for this pass:
- Focus: Python orchestration overhead reductions around native kernels (`fib_sphere.py`, `voxel.py`).
- Benchmark harness: `scripts/benchmark.py`
- Resolution: `2.0 A`
- Repeats: `3`
- Group: `all` (`small`, `medium`, `large`)

Artifacts generated for this pass:
- `AGENTS/benchmark.python.orchestration.json`
- `AGENTS/benchmark.ctypes-c.orchestration.json`
- `AGENTS/benchmark.native.orchestration.json`

Environment snapshot:
- Date generated (UTC): `2026-02-19 20:42:35 UTC`
- OS/Kernel: `Linux 6.18.8-200.fc43.x86_64 x86_64 GNU/Linux`
- CPU: `Intel(R) Core(TM) i9-10900X CPU @ 3.70GHz` (`20` logical CPUs)
- RAM: `125.47 GiB`
- Python: `3.11.14`

### 9.1 Runtime Comparison (Current)

| case | python_s | ctypes_s | native_s | native_vs_python | ctypes_vs_python | native_vs_ctypes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cavity | 0.388639 | 0.166743 | 0.099065 | 3.923x faster | 2.331x faster | 1.683x |
| pocket | 0.524975 | 0.173382 | 0.089402 | 5.872x faster | 3.028x faster | 1.939x |
| pore | 0.713359 | 0.174720 | 0.085016 | 8.391x faster | 4.083x faster | 2.055x |
| hub | 0.716801 | 0.196823 | 0.095459 | 7.509x faster | 3.642x faster | 2.062x |
| 4jpn | 82.248926 | 8.613003 | 2.902769 | 28.335x faster | 9.549x faster | 2.967x |
| 4jpp_assembly | 206.498628 | 24.069907 | 8.778116 | 23.524x faster | 8.579x faster | 2.742x |

Interpretation of `native_vs_ctypes`:
- values `< 1.0` mean native is slower than ctypes-c.
- values `> 1.0` mean native is faster than ctypes-c.

### 9.2 Native Delta vs Prior Optimized Native Baseline

Reference baseline artifact:
- `AGENTS/benchmark.native.optimized.json` (report dated `2026-02-16`)

| case | native_prev_s | native_now_s | improvement |
| --- | ---: | ---: | ---: |
| cavity | 0.104879 | 0.099065 | 5.54% |
| pocket | 0.099812 | 0.089402 | 10.43% |
| pore | 0.097166 | 0.085016 | 12.50% |
| hub | 0.104784 | 0.095459 | 8.90% |
| 4jpn | 2.888686 | 2.902769 | -0.49% |
| 4jpp_assembly | 9.264783 | 8.778116 | 5.25% |

Notes:
- Small-case native runtimes improved by ~`5.5%` to `12.5%`.
- Large case (`4jpp_assembly`) improved by `5.25%`.
- Medium case (`4jpn`) shows a small regression (`0.49%`), within normal benchmark variance observed across repeated runs.

### 9.3 Memory Snapshot

| case | python_mb | ctypes_mb | native_mb |
| --- | ---: | ---: | ---: |
| cavity | 177.148 | 176.945 | 176.785 |
| pocket | 175.383 | 175.742 | 175.367 |
| pore | 172.875 | 172.496 | 172.652 |
| hub | 171.902 | 171.551 | 171.953 |
| 4jpn | 698.867 | 698.816 | 698.965 |
| 4jpp_assembly | 1714.551 | 1714.582 | 1714.410 |

No meaningful memory regression observed from this pass.

### 9.4 Output Consistency Check

Across `python`, `ctypes-c`, and `native` runs for this pass, each case matched on:
- `num_detected_volumes`
- `largest_type`
- `largest_volume`

This indicates no high-level classification drift from the orchestration-overhead changes.

## 10. Stage-Level Profile Snapshot (2026-02-19)

Scope:
- Goal: identify remaining runtime concentration after Python orchestration-overhead reductions.
- Command used:
  - `VOLUMIZER_BACKEND=native uv run --python 3.11 python scripts/benchmark.py --group all --repeats 1 --resolution 2.0 --print-stage-breakdown --stage-top-n 6 --output-json AGENTS/benchmark.native.stage-profiling.all.sample.json`
- Artifact:
  - `AGENTS/benchmark.native.stage-profiling.all.sample.json`

Top stage contributors by case (native backend):

- `4jpp_assembly` (`8.799919s` total)
  - `classify_components`: `6.920783s` (`78.646%`)
  - `prepare_structure`: `1.565332s` (`17.788%`)
  - `load_structure`: `0.098293s` (`1.117%`)
- `4jpn` (`2.882467s` total)
  - `classify_components`: `2.137557s` (`74.157%`)
  - `prepare_structure`: `0.541521s` (`18.787%`)
  - `load_structure`: `0.069114s` (`2.398%`)

Interpretation:
- The dominant optimization target is still component classification (`classify_components`) on medium/large cases.
- The second largest bucket is structure preparation (`prepare_structure`, mostly cleaning + principal-axis alignment), especially visible once classification is fast.
- Remaining stages are comparatively small (<~2.5% each on medium/large), so expected ROI is lower there.

Recommended next optimization focus:
1. Reduce `classify_components` cost further (native algorithmic/data-structure tuning and reduced Python mapping overhead after native return).
2. Investigate `prepare_structure` cost (avoid redundant work, evaluate optional alignment strategies, or selectively skip alignment where safe).

## 11. Native Classification Split Profile (Kernel vs Mapping) (2026-02-19)

Scope:
- Goal: split `classify_components` into native kernel time vs Python-side mapping time.
- Command used:
  - `VOLUMIZER_BACKEND=native uv run --python 3.11 python scripts/benchmark.py --group all --repeats 1 --resolution 2.0 --print-stage-breakdown --stage-top-n 8 --output-json AGENTS/benchmark.native.stage-profiling.split.all.sample.json`
- Artifact:
  - `AGENTS/benchmark.native.stage-profiling.split.all.sample.json`

Key medium/large results:
- `4jpp_assembly` (`8.882311s` total)
  - `classify_components_native_kernel`: `3.976869s` (`44.773%`)
  - `classify_components_native_mapping`: `2.996221s` (`33.732%`)
  - `prepare_structure`: `1.589959s` (`17.900%`)
- `4jpn` (`3.020501s` total)
  - `classify_components_native_kernel`: `1.667231s` (`55.197%`)
  - `classify_components_native_mapping`: `0.573792s` (`18.997%`)
  - `prepare_structure`: `0.573989s` (`19.003%`)

Interpretation:
- `classify_components` remains the dominant bucket, but the split shows meaningful Python-side mapping cost after native return, especially on large structures.
- For `4jpp_assembly`, mapping is ~`75%` of kernel time (`2.996s` vs `3.977s`), so optimizing only kernel code will leave substantial time on the table.
- Next highest non-classification bucket is still `prepare_structure` (~`18-19%` on medium/large).

Immediate optimization targets informed by split profile:
1. Reduce Python mapping overhead in `classify_buried_components_native` path (object creation, set conversions, component assembly).
2. Continue native kernel tuning for classifier core loops.
3. Evaluate cost/necessity of full `prepare_structure` alignment in workflows where strict alignment is not required.

## 12. Mapping Overhead Reduction via Fast Axial-Length Kernel (2026-02-20)

Scope:
- Goal: reduce `classify_components_native_mapping` overhead identified in Section 11.
- Change summary:
  - Replaced per-component axial-length alignment call with direct NumPy PCA/eigendecomposition in `volumizer/voxel.py`.
  - Kept classification semantics and output schema unchanged (validated by full test suite).

Artifacts used for before/after comparison (both `repeats=3`, native backend):
- Before:
  - `AGENTS/benchmark.native.stage-profiling.split.medium.r3.v4.json`
  - `AGENTS/benchmark.native.stage-profiling.split.large.r3.v4.json`
- After:
  - `AGENTS/benchmark.native.stage-profiling.split.medium.r3.v5.json`
  - `AGENTS/benchmark.native.stage-profiling.split.large.r3.v5.json`
  - `AGENTS/benchmark.native.stage-profiling.split.all.r3.v5.json`

### 12.1 Medium/Large Delta

| case | mean_s_before | mean_s_after | total_improvement | mapping_before_s | mapping_after_s | mapping_improvement |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 4jpn | 3.741052 | 2.359652 | 36.93% | 1.357268 | 0.011684 | 99.14% |
| 4jpp_assembly | 7.153828 | 5.874129 | 17.89% | 1.296058 | 0.018259 | 98.59% |

### 12.2 Current Stage Distribution (Native, repeats=3)

From `AGENTS/benchmark.native.stage-profiling.split.all.r3.v5.json`:

- `4jpn` (`2.542518s` mean)
  - `classify_components_native_kernel`: `1.698309s` (`66.796%`)
  - `prepare_structure`: `0.603458s` (`23.735%`)
  - `classify_components_native_mapping`: `0.013243s` (`0.521%`)
- `4jpp_assembly` (`5.939490s` mean)
  - `classify_components_native_kernel`: `3.947005s` (`66.454%`)
  - `prepare_structure`: `1.643937s` (`27.678%`)
  - `classify_components_native_mapping`: `0.017856s` (`0.301%`)

### 12.3 Updated Optimization Priority

With mapping overhead now near-zero on medium/large cases, next highest-ROI work is:
1. Native classifier kernel optimization (`classify_components_native_kernel`).
2. Structure preparation cost reduction (`prepare_structure`: cleaning + principal-axis alignment).

## 13. Prepare-Structure Overhead Reduction via Thin-SVD Principal Alignment (2026-02-20)

Scope:
- Goal: reduce `prepare_structure` overhead identified in Section 12.3.
- Change summary:
  - Replaced `align.align_structure()` passthrough to `biotite.structure.orient_principal_components()` with an equivalent local implementation that uses thin SVD (`full_matrices=False`).
  - Preserved orientation semantics (same rotation loop/tolerance/order/sign handling) while avoiding expensive allocation of an `N x N` left-singular matrix for `N x 3` coordinates.
  - Added alignment parity regression test (`tests/test_align.py`) to assert coordinate-level equality against Biotite on representative small/medium/large structures.

Artifacts used for before/after comparison (both `repeats=3`, native backend):
- Before:
  - `AGENTS/benchmark.native.stage-profiling.split.medium.r3.v5.json`
  - `AGENTS/benchmark.native.stage-profiling.split.large.r3.v5.json`
- After:
  - `AGENTS/benchmark.native.stage-profiling.split.medium.r3.v6.json`
  - `AGENTS/benchmark.native.stage-profiling.split.large.r3.v6.json`
  - `AGENTS/benchmark.native.stage-profiling.split.all.r3.v6.json`

### 13.1 Medium/Large Delta

| case | mean_s_before | mean_s_after | total_improvement | prepare_before_s | prepare_after_s | prepare_improvement |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 4jpn | 2.359652 | 1.846656 | 21.74% | 0.533345 | 0.002122 | 99.60% |
| 4jpp_assembly | 5.874129 | 4.206156 | 28.40% | 1.570291 | 0.003052 | 99.81% |

### 13.2 Current Stage Distribution (Native, repeats=3)

From `AGENTS/benchmark.native.stage-profiling.split.all.r3.v6.json`:

- `4jpn` (`1.798947s` mean)
  - `classify_components_native_kernel`: `1.583837s` (`88.042%`)
  - `classify_components_native_mapping`: `0.011675s` (`0.649%`)
  - `prepare_structure`: `0.002063s` (`0.115%`)
- `4jpp_assembly` (`4.274768s` mean)
  - `classify_components_native_kernel`: `3.933186s` (`92.009%`)
  - `classify_components_native_mapping`: `0.017876s` (`0.418%`)
  - `prepare_structure`: `0.002920s` (`0.068%`)

### 13.3 Memory Snapshot

Compared to the prior split-profile repeat artifacts (`v5`), peak RSS dropped substantially:

| case | max_rss_before_mb | max_rss_after_mb | change |
| --- | ---: | ---: | ---: |
| 4jpn | 699.000 | 181.434 | -74.04% |
| 4jpp_assembly | 1715.316 | 193.277 | -88.73% |

This matches the thin-SVD change: avoiding allocation of large dense `U` matrices in principal-axis alignment removes a major transient memory spike.

### 13.4 Updated Optimization Priority

With both Python mapping and prepare-structure overhead now near-zero on medium/large cases, next highest-ROI work is:
1. Native classifier kernel optimization (`classify_components_native_kernel`).
2. Secondary pipeline stages now visible after kernel work (`load_structure`, `get_first_shell_exposed_voxels`, `volumes_to_structure`).
