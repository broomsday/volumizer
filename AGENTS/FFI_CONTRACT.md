# Volumizer Native FFI Contract

## 1. Purpose

Define stable boundaries between Python orchestration and the future Rust core so we can migrate kernels incrementally without changing user-facing API behavior.

This contract targets:
- `fib_sphere` bulk point generation
- voxel neighbor detection
- BFS component expansion
- exposed/buried solvent split
- (later) full volume classification

Current implementation status:
- Implemented in `volumizer_native`:
- `fibonacci_sphere_points`
- `get_neighbor_voxel_indices`
- `bfs_component_indices`
- `get_exposed_and_buried_voxel_indices`
- `classify_buried_components`
- Still pending:
- Stable packaging/import contract for installed wheels across platforms

## 2. Conventions

- Coordinate system: Angstrom units, right-handed XYZ.
- Index dtype: `int32`.
- Coordinate dtype: `float32` in native calls unless otherwise noted.
- Contiguous arrays required at boundary (`C` order).
- All return buffers are deterministic and sorted where relevant.
- Errors: native raises Python `ValueError` for invalid inputs.

## 3. Backend Selection

Python runtime should support:
- `VOLUMIZER_BACKEND=python` (default) -> always Python implementation
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

## 4.2 `get_neighbor_voxel_indices`

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

## 4.3 `bfs_component_indices`

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

## 4.4 `get_exposed_and_buried_voxel_indices`

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

## 4.5 `classify_buried_components` (Phase 2+)

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
- `component_centers: np.ndarray[float32]` shape `(C,3)`
- `component_axial_lengths: np.ndarray[float32]` shape `(C,3)`

Notes:
- Flattened buffers avoid Python object allocation overhead.
- Python layer will map type codes to labels and construct `VoxelGroup`.

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
