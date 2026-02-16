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
