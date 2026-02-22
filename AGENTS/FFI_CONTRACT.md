# Volumizer Native FFI Contract

## 1. Purpose

Define stable boundaries between Python orchestration and the future Rust core so we can migrate kernels incrementally without changing user-facing API behavior.

This contract targets:
- `fib_sphere` bulk point generation
- voxel neighbor detection
- BFS component expansion
- exposed/buried solvent split
- first-shell exposed voxel selection
- full volume classification

Current implementation status:
- Implemented in `volumizer_native`:
- `fibonacci_sphere_points`
- `fibonacci_sphere_points_batch`
- `get_neighbor_voxel_indices`
- `get_first_shell_exposed_indices`
- `get_first_shell_exposed_selection`
- `bfs_component_indices`
- `get_exposed_and_buried_voxel_indices`
- `get_exposed_and_buried_selection`
- `classify_buried_components`
- Still pending:
- Stable packaging/import contract for installed wheels across platforms

## 2. Conventions

- Coordinate system: Angstrom units, right-handed XYZ.
- Index dtype: `int32` for `(N,3)` voxel/index arrays unless otherwise noted.
- Split-axis selection APIs use `int64` axis arrays and return `int64` flat indices.
- Coordinate dtype: `float32` in native calls unless otherwise noted.
- Contiguous arrays required at boundary (`C` order).
- All return buffers are deterministic and sorted where relevant.
- Errors: native raises Python `ValueError` for invalid inputs.

## 3. Backend Selection

Python runtime should support:
- `VOLUMIZER_BACKEND=python` (default) -> Python orchestration path (may still use `ctypes` C helpers when present)
- `VOLUMIZER_BACKEND=native` -> require native module, fail loudly if unavailable
- `VOLUMIZER_BACKEND=auto` -> use native when import succeeds, otherwise Python

## 4. Function Contracts

## 4.1 `fibonacci_sphere_points`

Intent:
- Generate `samples` points on a sphere centered at `(x,y,z)` with radius `radius`.

Input:
- `radius: float`
- `x: float`
- `y: float`
- `z: float`
- `samples: int`

Output:
- `np.ndarray[float32]` with shape `(samples, 3)` (columns `x,y,z`)

Notes:
- Semantics must match current `volumizer.fib_sphere.fibonacci_sphere`.

## 4.2 `fibonacci_sphere_points_batch`

Intent:
- Generate sphere points for many centers in one call to reduce Python/native call overhead.

Input:
- `radius: float`
- `centers: np.ndarray[float32]` shape `(N, 3)`
- `samples: int`

Output:
- `np.ndarray[float32]` with shape `(N * samples, 3)`

Notes:
- Output ordering is deterministic:
- all `samples` points for `centers[0]`, then `centers[1]`, etc.
- Semantics per center must match `fibonacci_sphere_points`.

## 4.3 `get_neighbor_voxel_indices`

Intent:
- Return indices of query voxels that are 6-neighbors of any reference voxel.

Input:
- `query_voxels: np.ndarray[int32]` shape `(N, 3)`
- `reference_voxels: np.ndarray[int32]` shape `(M, 3)`

Output:
- `np.ndarray[int32]` shape `(K,)`, sorted ascending, containing indices into query array.

Notes:
- Neighbor rule is ordinal adjacency (Manhattan distance exactly 1).
- Must preserve parity with current `get_neighbor_voxels_*` behavior.

## 4.4 `bfs_component_indices`

Intent:
- Perform BFS over candidate voxel indices and return one connected component.

Input:
- `voxels: np.ndarray[int32]` shape `(N, 3)`
- `searchable_indices: np.ndarray[int32]` shape `(S,)` unique candidate indices

Output:
- `np.ndarray[int32]` shape `(C,)`, indices that belong to discovered component.

Notes:
- Connectivity rule is same 6-neighbor adjacency.
- Component seed semantics must match current Python/C implementations.

## 4.5 `get_exposed_and_buried_voxel_indices`

Intent:
- Classify solvent voxels into exposed vs buried based on axis-occlusion heuristics.

Input:
- `solvent_voxels: np.ndarray[int32]` shape `(S, 3)`
- `protein_voxels: np.ndarray[int32]` shape `(P, 3)`
- `grid_dims: np.ndarray[int32]` shape `(3,)`

Output:
- `dict` with:
- `exposed_indices: np.ndarray[int32]` shape `(E,)`, indices into `solvent_voxels`
- `buried_indices: np.ndarray[int32]` shape `(B,)`, indices into `solvent_voxels`

Notes:
- Must preserve parity with current `get_exposed_and_buried_voxels` semantics.
- Returned indices are deterministic and ordered by input traversal.
- This is the legacy native split API; Python prefers `get_exposed_and_buried_selection` when available.

## 4.6 `classify_buried_components`

Intent:
- Full native kernel for agglomerating buried voxels and classifying cavity/pocket/pore/hub.

Input:
- `buried_voxels: np.ndarray[int32]` shape `(B, 3)`
- `exposed_voxels: np.ndarray[int32]` shape `(E, 3)`
- `grid_dims: np.ndarray[int32]` shape `(3,)`
- `min_num_voxels: int`
- `voxel_size: float`

Output:
- Structured dict/tuple with:
- `component_type_codes: np.ndarray[int8]` (`0=occluded,1=cavity,2=pocket,3=pore,4=hub`)
- `component_offsets: np.ndarray[int32]` prefix offsets into `component_voxel_indices_flat`
- `component_voxel_indices_flat: np.ndarray[int32]`
- `surface_offsets: np.ndarray[int32]` prefix offsets into `component_surface_indices_flat`
- `component_surface_indices_flat: np.ndarray[int32]`
- optional `kernel_stage_timings_seconds: dict[str, float]` for benchmark instrumentation

Notes:
- Flattened buffers avoid Python object allocation overhead.
- Python layer maps type codes to labels, constructs `VoxelGroup`, and computes centers/axial lengths.

## 4.7 `get_first_shell_exposed_indices`

Intent:
- Return indices of exposed voxels that are 6-neighbors of any buried voxel.

Input:
- `exposed_voxels: np.ndarray[int32]` shape `(E, 3)`
- `buried_voxels: np.ndarray[int32]` shape `(B, 3)`
- `grid_dims: np.ndarray[int32]` shape `(3,)`

Output:
- `np.ndarray[int32]` shape `(K,)`, indices into `exposed_voxels`

Notes:
- Maintains first-shell semantics used by `voxel.get_first_shell_exposed_voxels()`.

## 4.8 `get_first_shell_exposed_selection`

Intent:
- Return first-shell selection indices and precomputed flat voxel indices to reduce Python-side mapping overhead.

Input:
- `exposed_x: np.ndarray[int64]` shape `(E,)`
- `exposed_y: np.ndarray[int64]` shape `(E,)`
- `exposed_z: np.ndarray[int64]` shape `(E,)`
- `buried_x: np.ndarray[int64]` shape `(B,)`
- `buried_y: np.ndarray[int64]` shape `(B,)`
- `buried_z: np.ndarray[int64]` shape `(B,)`
- `grid_dims: np.ndarray[int32]` shape `(3,)`

Output:
- `dict` with:
- `neighbor_indices: np.ndarray[int32]` shape `(K,)`, indices into exposed axis arrays
- `flat_indices: np.ndarray[int64]` shape `(K,)`, flat voxel-grid indices for selected voxels

Notes:
- Preferred fast path when available; Python falls back to `get_first_shell_exposed_indices`.

## 4.9 `get_exposed_and_buried_selection`

Intent:
- Return exposed/buried selections with query indices and flat voxel-grid indices in one call.

Input:
- `solvent_x: np.ndarray[int64]` shape `(S,)`
- `solvent_y: np.ndarray[int64]` shape `(S,)`
- `solvent_z: np.ndarray[int64]` shape `(S,)`
- `protein_x: np.ndarray[int64]` shape `(P,)`
- `protein_y: np.ndarray[int64]` shape `(P,)`
- `protein_z: np.ndarray[int64]` shape `(P,)`
- `grid_dims: np.ndarray[int32]` shape `(3,)`

Output:
- `dict` with:
- `exposed_query_indices: np.ndarray[int32]` shape `(E,)`
- `buried_query_indices: np.ndarray[int32]` shape `(B,)`
- `exposed_flat_indices: np.ndarray[int64]` shape `(E,)`
- `buried_flat_indices: np.ndarray[int64]` shape `(B,)`

Notes:
- Preferred low-copy native split API for Python mapping path.
- Python falls back to `get_exposed_and_buried_voxel_indices` when this API is unavailable.

## 5. Validation Requirements

- Unit parity tests for each function vs Python baseline.
- Deterministic output ordering verified in tests.
- Cross-check tolerances:
- volume exact parity on fixtures
- dimensions abs delta <= `0.01`

## 6. Versioning

- Track native ABI/API as `NATIVE_CONTRACT_VERSION`.
- Phase 1 target value: `1`.
- Any breaking signature change requires increment and wrapper update.
