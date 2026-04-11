# Volumizer Agent Guide

This is the single current-state guide for agents working in this repository.
It replaces the old migration/status/performance/gallery planning notes.
Keep this file concise and current. Do not reintroduce a large doc split unless
there is a strong reason.

## What This Project Does

`volumizer` analyzes 3D biomolecular structures and identifies internal solvent
volumes. The main reported classes are:

- `cavity`: no solvent-surface connection
- `pocket`: one solvent-surface connection
- `pore`: two solvent-surface connections
- `hub`: three or more solvent-surface connections

Small buried connected solvent groups are tracked as `occluded`.

Current CLI behavior:

- hub classification still exists in the core pipeline
- CLI-written annotation JSON/CIF outputs omit hubs by default
- `--include-hubs` restores hub emission for CLI outputs

Typical inputs:

- local `.pdb`, `.cif`, `.mmtf`
- downloaded RCSB `.cif` entries
- RCSB sequence-cluster representative sets via the CLI

Typical outputs:

- annotated structure file (`.annotated.cif` or PDB-formatted output)
- annotation JSON / dataframe-style records
- run summary JSON for batch jobs
- gallery SQLite index and cached render PNGs

## Current State

- The Rust native backend is the normal fast path.
- Python remains the API, CLI, I/O, and output-shaping layer.
- `VOLUMIZER_BACKEND=auto` is the default. If `volumizer_native` imports, the
  native path is used; otherwise the Python path runs.
- The gallery web service is working in MVP form:
  - SQLite index build
  - cached thumbnail rendering
  - FastAPI API + static frontend
  - Mol* viewer data endpoint and detail flow
- Gallery thumbnail rendering supports parallel structure-level jobs, local
  Mol* assets from `node_modules/molstar`, configurable
  `software|hardware|auto` backends, `compatibility|fast` axis modes, and
  optional `--timing-jsonl` telemetry.
- Scientific parity and deterministic output ordering matter more than local
  refactor elegance.

## Core Pipeline

The main orchestration entrypoint is `volumizer/volumizer.py`.

End-to-end flow:

1. Load structure: `pdb.load_structure()`
2. Clean structure: `pdb.clean_structure()`
3. Align to principal axes: `align.align_structure()`
4. Convert atoms to coordinates: `pdb.get_structure_coords()`
5. Expand atoms with shell points: `fib_sphere.add_extra_points()`
6. Build point cloud and voxel grid: `voxel.coords_to_point_cloud()`,
   `voxel.add_voxel_grid()`
7. Split protein vs solvent voxels:
   `voxel.get_protein_and_solvent_voxels()`
8. Split solvent into exposed vs buried:
   `voxel.get_exposed_and_buried_voxels()`
9. Keep only exposed voxels adjacent to buried voxels:
   `voxel.get_first_shell_exposed_voxels()`
10. Agglomerate buried components and classify them:
    `voxel.get_pores_pockets_cavities_occluded()`
11. Build annotation dataframe:
    `utils.make_annotation_dataframe()`
12. Build pseudo-atom visualization structure:
    `pdb.volumes_to_structure()`

## Where To Look

| Area | Start here | Important functions / notes |
| --- | --- | --- |
| Public pipeline API | `volumizer/volumizer.py` | `annotate_structure_volumes()`, `prepare_pdb_structure()`, `volumize_structure()`, `volumize_pdb()`, `volumize_pdb_and_save()` |
| Structure loading and output | `volumizer/pdb.py` | `load_structure()`, `clean_structure()`, `get_structure_coords()`, `volumes_to_structure()`, `make_volumized_pdb_lines()`, `save_structure()` |
| Alignment | `volumizer/align.py` | `align_structure()` uses the current principal-axis alignment logic |
| Atom shell expansion | `volumizer/fib_sphere.py` | `fibonacci_sphere()`, `estimate_fibonacci_sphere_samples()`, `add_extra_points()`, `add_extra_points_native()`, `add_extra_points_c()` |
| Voxelization and classification | `volumizer/voxel.py` | `get_exposed_and_buried_voxels()`, `get_first_shell_exposed_voxels()`, `breadth_first_search()`, `classify_buried_components_native()`, `get_pores_pockets_cavities_occluded()`, cross-section metric helpers |
| Runtime backend selection | `volumizer/native_backend.py` | `get_requested_backend()`, `get_native_module()`, `using_native()`, `active_backend()` |
| Global config / summaries | `volumizer/utils.py` | `set_resolution()`, `set_non_protein()`, `make_annotation_dataframe()`, backend helper wrappers |
| Shared types | `volumizer/types.py` | `VoxelGroup`, `Annotation`, `ComponentData` |
| Constants / thresholds | `volumizer/constants.py` | radii, voxel/type mappings, default thresholds |
| RCSB download + cluster helpers | `volumizer/rcsb.py` | `download_structure_cif()`, `fetch_entry_metadata()`, `entry_passes_filters()`, `fetch_cluster_representative_entry_ids()` |
| CLI entrypoint | `volumizer/cli.py` | `build_parser()`, `resolve_input_structures()`, `analyze_structure_file()`, `_analyze_structures()`, `run_cli()` |
| Gallery index build | `volumizer/gallery_index.py` | `build_gallery_index()` parses run outputs and writes `data/gallery.db` |
| Gallery querying | `volumizer/gallery_query.py` | `query_gallery_index()` builds SQL-backed filtered result sets |
| Gallery thumbnail rendering | `volumizer/gallery_render.py` | `render_gallery_thumbnails()` coordinates Node/Mol* render jobs |
| Gallery web service | `volumizer/web/app.py` | `create_app()` wires FastAPI routes and static UI |
| Gallery DB access | `volumizer/web/db.py` | `list_runs()`, `get_hit_detail()`, `get_hit_artifacts()` |
| Native Rust module | `native/src/lib.rs` | Python exports include `fibonacci_sphere_points`, `fibonacci_sphere_points_batch`, `get_neighbor_voxel_indices`, `get_first_shell_exposed_selection`, `get_exposed_and_buried_selection`, `classify_buried_components` |
| Legacy C helpers | `src/fib_sphere.c`, `src/voxel.c` | still supported as fallback/compatibility path, but not the main optimization target |

## Backend Model

Backend control lives in `volumizer/native_backend.py`.

Supported modes:

- `VOLUMIZER_BACKEND=auto`: default, prefer Rust if importable
- `VOLUMIZER_BACKEND=native`: require Rust and fail loudly if unavailable
- `VOLUMIZER_BACKEND=python`: force Python path

Practical guidance:

- Treat the Rust path as canonical for performance work.
- Keep Python and native semantics aligned.
- Avoid adding per-item Python/native calls; batch contiguous-array boundaries
  are the intended design.
- If behavior changes in classification, ordering, or geometric metrics, expect
  tests and downstream tooling to break.

## Structure Loading Notes

`pdb.load_structure()` is a high-leverage function.

- It handles local structure parsing and assembly policy.
- Assembly policy options are `biological`, `asymmetric`, and `auto`.
- This is still one of the most important places to inspect when runtime moves.
- Loading timings are already instrumented in sub-stages, so performance work
  here should preserve timing visibility.

## CLI Surface

Main subcommands:

- `volumizer analyze`
- `volumizer cluster`
- `volumizer cache`

Capabilities already in place:

- local-file analysis
- PDB ID download flow
- cluster representative workflows
- metadata cache and negative cache
- checkpoint/resume
- progress JSONL and periodic ETA output
- manifest write/replay
- replay from prior summaries
- deterministic sharding
- failures-manifest export
- assembly-policy control

If a user asks about batch execution, start in `volumizer/cli.py`.

## Output Artifacts

Primary analysis outputs:

- annotated structure:
  - pseudo-atoms encode detected volume voxels
  - residue names encode type (`POR`, `POK`, `CAV`, `OCC`)
  - `HUB` is only written when CLI hub emission is enabled
- annotation records:
  - per-type counts
  - per-volume sizes
  - axial lengths
  - cross-section circularity / uniformity
- run summary JSON for CLI jobs

Primary gallery outputs:

- `data/gallery.db`
- `data/runs/<run_id>/...`
- `data/renders/<run_id>/<source_label>/x.png|y.png|z.png`

## Gallery MVP

The gallery stack is local-first and intentionally simple.

Data flow:

1. Run volumizer analysis (`analyze` or `cluster`)
2. Build SQLite index from run outputs:
   `volumizer/gallery_index.py`
3. Generate cached renders:
   `volumizer/gallery_render.py`
   `scripts/molstar_render_single.mjs`
4. Serve local API and frontend:
   `volumizer/web/app.py`

Current responsibilities:

- `gallery_index.py`:
  - ingests summary + annotation outputs
  - computes structure-level metrics
  - stores per-volume and aggregate values in SQLite
- `gallery_query.py`:
  - filter + sort + pagination logic
- `gallery_render.py`:
  - dispatches Mol* thumbnail rendering through Node
- `web/app.py`:
  - serves `/api/runs`
  - serves `/api/hits`
  - serves `/api/hits/{structure_id}`
  - serves `/api/hits/{structure_id}/viewer-data`
  - serves structure, annotation, and thumbnail files

The frontend is intentionally thin. Avoid turning it into a large SPA unless
the product direction materially changes.

## Native Rust Surface

The Rust extension is in `native/` and exposed as `volumizer_native`.

Key exported functions:

- `contract_version()`
- `backend_info()`
- `fibonacci_sphere_points()`
- `fibonacci_sphere_points_batch()`
- `get_neighbor_voxel_indices()`
- `get_first_shell_exposed_indices()`
- `get_first_shell_exposed_selection()`
- `bfs_component_indices()`
- `get_exposed_and_buried_voxel_indices()`
- `get_exposed_and_buried_selection()`
- `classify_buried_components()`

When working on native code:

- preserve deterministic ordering
- preserve Python-visible shapes and dtypes
- prefer flat arrays and offsets over Python object construction
- keep parity tests close to any behavior change

## Tests And Benchmarks

Useful test clusters:

- core behavior:
  - `tests/test_voxel.py`
  - `tests/test_fib_sphere.py`
  - `tests/test_pdb.py`
- native parity / backend:
  - `tests/test_native_backend.py`
  - `tests/test_pipeline_native_parity.py`
  - `tests/test_voxel_native_parity.py`
  - `tests/test_fib_sphere_native_parity.py`
- CLI:
  - `tests/test_cli.py`
- gallery:
  - `tests/test_gallery_index.py`
  - `tests/test_gallery_query.py`
  - `tests/test_gallery_render.py`
  - `tests/test_gallery_web.py`

Performance harness:

- `scripts/benchmark.py`
- keep comparisons on the same machine
- keep resolution fixed when comparing before/after
- benchmark changes that touch `load_structure`, voxel classification, or
  gallery indexing/render hot paths

## Working Rules For Agents

- Start from this file, then open the owning module directly.
- Prefer reading code over relying on old assumptions.
- Preserve classification semantics unless the task explicitly changes them.
- Preserve output stability where IDs, ordering, and saved artifact names matter.
- Treat `load_structure()` and `voxel.py` as high-risk areas for silent drift.
- Treat the gallery DB schema and API response shapes as compatibility surfaces.
- Keep documentation centralized here; if the architecture changes, update this
  file instead of reviving long status/migration narratives.
- Remove temporary planning notes once the work lands; fold durable guidance
  back into this file.
