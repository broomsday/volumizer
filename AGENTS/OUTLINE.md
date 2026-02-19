# Volumizer Project Outline

## 1. Project Purpose

`volumizer` analyzes 3D biomolecular structures (primarily from the PDB/RCSB ecosystem) to find solvent-occluded internal volumes and classify them by how they connect to the external solvent surface.

Primary volume classes:
- `cavity`: no surface connection
- `pocket`: one surface connection
- `pore`: two surface connections
- `hub`: three or more surface connections

The package also tracks small/internal connected solvent groups as `occluded` volumes.

## 2. End-to-End Workflow

Main orchestration is in `volumizer/volumizer.py`.

1. Load a structure from `.pdb`, `.cif`, `.mmtf`, or other formats (`volumizer/pdb.py`).
2. Clean structure (optionally remove non-protein residues, hydrogens; handle multi-model input).
3. Align coordinates to principal axes (`volumizer/align.py`).
4. Convert atoms to coordinate dataframe (`x`, `y`, `z`, `element`).
5. Expand each atom with extra shell points (Fibonacci sphere sampling; `volumizer/fib_sphere.py`) to better represent atom volume at chosen voxel size.
6. Build point cloud and voxel grid (`volumizer/voxel.py`, PyntCloud backend).
7. Split voxels into protein-filled vs solvent voxels.
8. Classify solvent voxels as exposed or buried using axis-occlusion heuristics.
9. Keep exposed voxels adjacent to buried voxels as the relevant external “first shell”.
10. Agglomerate buried voxels into connected components (BFS).
11. Classify each component by number of distinct solvent-contact surfaces:
    - 0 => cavity
    - 1 => pocket
    - 2 => pore
    - 3+ => hub
12. Compute per-volume metrics (volume, center, principal-axis dimensions).
13. Emit:
    - annotation dataframe (`id`, `type`, `volume`, `x`, `y`, `z`)
    - structure containing volume voxels encoded as pseudo-atoms for visualization.

## 3. Classification Logic (Current Implementation)

`volumizer/voxel.py` is the computational core.

- Voxel connectivity uses 6-neighbor (axis-adjacent) logic.
- Buried/exposed determination is based on whether a solvent voxel is blocked on opposing axis directions by protein voxels.
- Volume type is inferred from how many disconnected solvent-contact surface regions touch each buried component.
- Very small connected groups are labeled `occluded` before final reporting.

## 4. Output Artifacts

### Annotation DataFrame
Built in `volumizer/utils.py::make_annotation_dataframe()`:
- one row per detected hub/pore/pocket/cavity
- sorted by descending volume
- stores principal-axis lengths as `x`, `y`, `z`

### Annotated PDB-style Structure
Built in `volumizer/pdb.py`:
- each voxel is emitted as a pseudo-atom
- residue name encodes volume type (`HUB`, `POR`, `POK`, `CAV`, `OCC`)
- B-factor marks direct surface-contact voxels (`50.0` vs `0.0`)

## 5. Codebase Map

- `volumizer/volumizer.py`: high-level API (`volumize_pdb`, `volumize_structure`, save helpers)
- `volumizer/cli.py`: CLI parser/orchestration (`analyze`, `cluster`, `cache` subcommands, checkpoint/resume, progress logging (JSONL + periodic ETA), metadata cache controls, manifest write/replay, summary subset replay, deterministic cluster sharding, failures-manifest export)
- `volumizer/voxel.py`: voxel generation, burial/exposure split, BFS agglomeration, type assignment, geometric metrics
- `volumizer/fib_sphere.py`: atom shell point generation
- `volumizer/pdb.py`: structure I/O, cleanup, PDB formatting, annotation structure assembly
- `volumizer/utils.py`: runtime configuration and annotation helpers
- `volumizer/constants.py`: radii, thresholds, voxel/type mappings
- `src/voxel.c`, `src/fib_sphere.c`: optional accelerated routines called via `ctypes`
- `tests/`: behavior and regression tests for classification and C/Python parity

## 6. Performance Architecture (Python + C)

The project uses Python by default and conditionally enables native acceleration if both shared libraries exist:
- `src/voxel.so`
- `src/fib_sphere.so`

These are manually built via `src/compile_c_libs.sh`. Runtime selection is controlled by `utils.using_performant()`.

Accelerated/native paths currently cover:
- Fibonacci sphere coordinate generation
- Neighbor voxel detection
- BFS component growth

Most orchestration, data reshaping, and classification logic remains in Python.

## 6a. Migration Status (Current)

Implemented so far:
- Optional Rust extension scaffold exists in `native/` (`pyo3` + `maturin`).
- Native backend selection exists via `VOLUMIZER_BACKEND` in `volumizer/native_backend.py`.
- Native kernels currently implemented:
- Fibonacci sphere point generation.
- Fibonacci sphere batch point generation.
- Neighbor voxel index detection.
- BFS component expansion.
- Exposed/buried solvent split.
- Buried-component classification output in flattened buffer format.
- Python dispatch/fallback paths are wired in `volumizer/fib_sphere.py` and `volumizer/voxel.py`.
- Native parity and integration tests are present and passing when the module is built.

Still pending for full switch-over:
- Further reduce Python-side dataframe/orchestration overhead around native paths.
- Package native module for standard install path (prebuilt wheels/CI matrix).
- Maintain and refresh backend-vs-backend performance reporting.

## 7. Important Configuration Knobs

In `volumizer/utils.py`:
- `set_resolution(resolution)`: controls voxel size and voxel volume globally
- `set_non_protein(bool)`: include/exclude non-protein residues

Default behavior emphasizes robust protein-focused analysis:
- voxel size: `3.0 Å`
- non-protein excluded
- hydrogens excluded

## 8. Recode Direction (C/Rust-Oriented)

For a C/Rust migration, the highest-impact targets are:

1. `volumizer/voxel.py` component labeling and surface classification loops.
2. `volumizer/fib_sphere.py` atom expansion loops and dataframe-heavy concatenation.
3. Data representation conversions between pandas, numpy, and PyntCloud.

Likely architecture for next iteration:
- keep Python API layer and output formatting initially
- move core volumization kernel to a native library (Rust or C)
- expose minimal FFI boundary with contiguous arrays (avoid per-voxel Python calls)
- preserve current classification semantics to keep test fixtures valid

## 9. UI / Usability Direction

Current user interface includes both Python function calls and a subcommand CLI (`volumizer analyze|cluster|cache`) for file/PDB-ID/cluster inputs with CIF + JSON outputs plus metadata-cache inspection/maintenance. Legacy flag-only invocation is still auto-routed for compatibility. For broader adoption:

1. Continue CLI UX hardening with validation and safety features (for example manifest/summary validation subcommands, failure-threshold guards, and runtime limits) on top of the current baseline (`analyze`/`cluster`/`cache`, resume, method/resolution filters, jobs/retries, metadata cache + negative cache, dry-run, checkpoint, progress JSONL + progress interval, manifest write/replay, summary subset replay, deterministic sharding, failures-manifest export).
2. Improve install ergonomics for native acceleration (prebuilt wheels or optional Rust extension build).
3. Add user-facing docs focused on:
   - quick start from raw PDB/CIF
   - interpretation of volume classes
   - visualization workflow (PyMOL/ChimeraX)
