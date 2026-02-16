# Volumizer Backend Performance Report (2026-02-16)

## 1. Scope

- Benchmark harness: `scripts/benchmark.py`
- Resolution: `2.0 A`
- Repeats: `3`
- Group: `all` (`small`, `medium`, `large`)
- Backend runs:
- `python`: `VOLUMIZER_BACKEND=python` with `src/voxel.so` and `src/fib_sphere.so` temporarily moved out
- `ctypes-c`: `VOLUMIZER_BACKEND=python` with compiled `.so` accelerators present
- `native`: `VOLUMIZER_BACKEND=native` with `volumizer_native` built via `maturin`

Artifacts:
- `AGENTS/benchmark.python.json`
- `AGENTS/benchmark.ctypes-c.json`
- `AGENTS/benchmark.native.json`

## 2. Environment Snapshot

- Date generated (UTC): `2026-02-16 04:32:27 UTC`
- OS: `Fedora Linux 43 (Xfce)`
- Kernel: `6.18.8-200.fc43.x86_64`
- CPU: `Intel(R) Core(TM) i9-10900X CPU @ 3.70GHz` (`20` logical CPUs)
- RAM: `131563892 kB` (~`125.47 GiB`)
- Python: `3.11.14`

## 3. Runtime Comparison

Mean runtime in seconds from each backend summary:

| case | python_s | ctypes_s | native_s | native_vs_python | ctypes_vs_python | native_vs_ctypes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cavity | 0.394818 | 0.169685 | 0.196238 | 2.012x faster | 2.327x faster | 0.865x |
| pocket | 0.539714 | 0.205055 | 0.217411 | 2.482x faster | 2.632x faster | 0.943x |
| pore | 0.760463 | 0.189863 | 0.259198 | 2.934x faster | 4.005x faster | 0.733x |
| hub | 0.717067 | 0.189071 | 0.260075 | 2.757x faster | 3.793x faster | 0.727x |
| 4jpn | 83.107042 | 8.548185 | 20.007940 | 4.154x faster | 9.722x faster | 0.427x |
| 4jpp_assembly | 207.212250 | 24.393836 | 57.477134 | 3.605x faster | 8.494x faster | 0.424x |

Interpretation of `native_vs_ctypes`:
- values `< 1.0` mean native is slower than the current ctypes-C path.

## 4. Memory Snapshot

Peak RSS from summary rows is broadly similar across backends:
- medium (`4jpn`): ~`698-699 MB`
- large (`4jpp_assembly`): ~`1713-1714 MB`

No meaningful memory regression was observed for native in these runs.

## 5. Output Consistency Check

Across all three backend runs, each case matched on:
- `num_detected_volumes`
- `largest_type`
- `largest_volume`

This indicates no obvious high-level classification drift in benchmark outputs.

## 6. Success Criteria Check (Migration Plan)

Against the defined performance targets:
- Medium proteins (`4jpn`) native speedup target `>= 3x` vs python: **met** (`4.154x`)
- Large proteins (`4jpp_assembly`) native speedup target `>= 5x` vs python: **not met** (`3.605x`)

## 7. Current Conclusion

- Native backend is materially faster than pure Python.
- Native backend is currently significantly slower than the existing ctypes-C path (roughly `2.3x` slower on medium/large fixtures).
- Remaining optimization work should focus on reducing Python orchestration overhead around native kernels, especially in atom shell expansion and classification-adjacent loops.

