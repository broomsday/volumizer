# Volumizer Backend Performance Report

## 1. Summary

The native Rust backend (`volumizer_native`) delivers major speedups and memory reductions over both the pure-Python and ctypes-C backends across all benchmark cases. On medium and large structures the native backend is **392x** and **674x** faster than pure Python, and **40x** and **77x** faster than ctypes-C. Peak memory dropped **74%** on medium and **89%** on large structures.

All backends produce identical classification output (volume counts, types, and sizes) on every benchmark case.

## 2. Environment

All measurements were taken on a single host under consistent conditions.

- OS: Fedora Linux 43 (Xfce), kernel `6.18.8-200.fc43.x86_64`
- CPU: Intel Core i9-10900X @ 3.70 GHz (20 logical CPUs)
- RAM: 125.5 GiB
- Python: 3.11.14
- Benchmark harness: `scripts/benchmark.py`
- Resolution: 2.0 A, 3 repeats per case

## 3. Benchmark Cases

| category | case | input | detected volumes | largest type | largest volume (A^3) |
| --- | --- | --- | ---: | --- | ---: |
| small | cavity | `tests/pdbs/cavity.pdb` | 4 | cavity | 184 |
| small | pocket | `tests/pdbs/pocket.pdb` | 2 | pocket | 840 |
| small | pore | `tests/pdbs/pore.pdb` | 2 | pore | 1,904 |
| small | hub | `tests/pdbs/hub.pdb` | 2 | hub | 2,176 |
| medium | 4jpn | `tests/pdbs/4jpn.pdb` | 40 | pore | 45,048 |
| large | 4jpp_assembly | `tests/pdbs/4jpp.cif` | 66 | pore | 71,160 |

## 4. Runtime Comparison

Mean wall-clock seconds (3 repeats, 2.0 A resolution):

| case | python (s) | ctypes-c (s) | native (s) | vs python | vs ctypes-c |
| --- | ---: | ---: | ---: | ---: | ---: |
| cavity | 0.395 | 0.169 | 0.030 | 13x | 5.7x |
| pocket | 0.540 | 0.169 | 0.030 | 18x | 5.7x |
| pore | 0.760 | 0.171 | 0.030 | 25x | 5.6x |
| hub | 0.717 | 0.183 | 0.029 | 25x | 6.3x |
| 4jpn | 83.107 | 8.464 | 0.212 | 392x | 40x |
| 4jpp_assembly | 207.212 | 23.776 | 0.308 | 674x | 77x |

## 5. Memory Comparison

Peak resident set size (MB):

| case | python | ctypes-c | native | reduction vs python |
| --- | ---: | ---: | ---: | ---: |
| cavity | 176.6 | 177.2 | 163.6 | 7% |
| pocket | 174.6 | 175.3 | 163.6 | 6% |
| pore | 172.1 | 172.6 | 164.0 | 5% |
| hub | 171.6 | 171.9 | 163.6 | 5% |
| 4jpn | 698.4 | 698.3 | 182.0 | 74% |
| 4jpp_assembly | 1,714.3 | 1,714.8 | 191.9 | 89% |

The large memory reductions on medium/large structures came primarily from replacing full-matrix SVD with thin-SVD for principal-axis alignment, which eliminated a transient `N x N` allocation.

## 6. Optimization Journey

The native backend went through multiple optimization passes between Feb 13 and Feb 27, 2026. The initial native port was already faster than pure Python but only modestly faster than ctypes-C on small cases and slower on some medium/large stages. Successive passes targeted the dominant runtime stage at each step:

1. **Batch fib-sphere expansion + hash-set neighbor detection** reduced initial native runtime by 46-63% on small cases and 84-86% on medium/large. This eliminated per-center call overhead and replaced O(N*M) neighbor scans with O(1) hash lookups.

2. **Orchestration overhead reduction** removed unnecessary Python-side dataframe operations around native kernel calls, yielding 5-12% gains on small cases and ~5% on large.

3. **Fast axial-length kernel** replaced per-component alignment calls with direct NumPy PCA in the classification mapping path, reducing Python mapping overhead by 99% and total runtime by 17-37%.

4. **Thin-SVD principal alignment** replaced Biotite's full-matrix SVD with an equivalent thin-SVD implementation, reducing `prepare_structure` to near-zero and cutting peak memory by 74-89%.

5. **Lookup-based BFS + remaining-set elimination** replaced set-scan BFS with coordinate-lookup BFS and removed per-iteration set rebuilds in the Rust classifier kernel. This reduced kernel time by 76-82% and total runtime by 65-72%.

6. **Hash-lookup classify-type acceleration** replaced repeated full-array scans in component-type classification with hash-based adjacency checks and subset BFS. Kernel time dropped another 84-92%.

7. **Specialized native kernels** for first-shell exposed selection, volumes-to-structure vectorization, add-extra-points path simplification, and exposed/buried selection API each contributed incremental end-to-end gains of 2-17%.

Net improvement within the native backend itself across all passes:

| case | initial native (s) | final native (s) | improvement |
| --- | ---: | ---: | ---: |
| cavity | 0.196 | 0.030 | 85% |
| pocket | 0.217 | 0.030 | 86% |
| pore | 0.259 | 0.030 | 88% |
| hub | 0.260 | 0.029 | 89% |
| 4jpn | 20.008 | 0.212 | 99% |
| 4jpp_assembly | 57.477 | 0.308 | 99% |

## 7. Final Stage Distribution

With all optimizations applied, runtime is broadly distributed across stages on medium and large cases. No single stage dominates by more than ~32%.

### 4jpn (medium, 0.212s total)

| stage | time (s) | share |
| --- | ---: | ---: |
| load_structure | 0.069 | 32% |
| classify_components (kernel) | 0.045 | 21% |
| get_first_shell_exposed_voxels | 0.037 | 17% |
| add_extra_points | 0.017 | 8% |
| classify_components (mapping) | 0.011 | 5% |
| get_exposed_and_buried_voxels | 0.011 | 5% |
| add_voxel_grid | 0.008 | 4% |
| other | 0.014 | 7% |

### 4jpp_assembly (large, 0.308s total)

| stage | time (s) | share |
| --- | ---: | ---: |
| load_structure | 0.097 | 31% |
| classify_components (kernel) | 0.069 | 23% |
| get_first_shell_exposed_voxels | 0.046 | 15% |
| add_extra_points | 0.026 | 8% |
| classify_components (mapping) | 0.018 | 6% |
| get_exposed_and_buried_voxels | 0.018 | 6% |
| add_voxel_grid | 0.013 | 4% |
| other | 0.021 | 7% |

## 8. Success Criteria

Against the targets defined in `AGENTS/MIGRATION_PLAN.md`:

| criterion | target | result | status |
| --- | --- | --- | --- |
| Medium (4jpn) speedup vs Python | >= 3x | 392x | met |
| Large (4jpp) speedup vs Python | >= 5x | 674x | met |
| Peak memory increase | <= 1.5x | 0.11x (89% reduction) | met |
| Volume classification parity | exact match | all cases match | met |
| Axial length parity | abs delta <= 0.01 | within tolerance | met |

## 9. Output Consistency

Across all three backends (`python`, `ctypes-c`, `native`), every benchmark case produces identical results for:
- `num_detected_volumes`
- `largest_type`
- `largest_volume`

No classification drift was observed at any point during the optimization work.
