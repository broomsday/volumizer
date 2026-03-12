# RCSB Hit Gallery + Mol* Viewer Plan (2026-03-02)

## 1. Goal
Build a local web workflow on top of volumizer outputs that can:
1. Filter hits by pore metrics (size, length, min/max diameter).
2. Filter hits by structure metrics (chain count, residue count, sequence-unique chains).
3. Show a gallery of cached quick renders (3 axis-aligned images per hit).
4. Open any hit in an interactive 3D molecular viewer using Mol*.

## 2. Key Decision: Tech for Steps 3 and 4

- Java is not required.
- Browser-side JavaScript is required for Mol* interactive viewing.
- TypeScript is optional:
  - For a minimal implementation, plain JS is enough.
  - For a larger maintainable frontend, TypeScript is recommended.

For step 3 pre-rendering, using Mol* is a good choice so gallery thumbnails match interactive viewer styling.

## 3. Proposed Architecture

### 3.1 Data Pipeline (Python + existing volumizer)
1. Use `volumizer analyze` / `volumizer cluster` to generate:
   - `*.annotated.cif`
   - `*.annotation.json`
   - `run.summary.json`
2. Run an index/enrichment step that reads those outputs and writes a local search index (SQLite preferred).

### 3.2 Render Pipeline (Node + Mol* + headless browser)
1. A renderer worker script loads each hit in Mol*.
2. Captures 3 PNGs per hit (`x`, `y`, `z` viewpoints).
3. Stores PNG paths/status in SQLite.
4. Re-runs only for missing/stale renders.

### 3.3 Local Web App
- Backend: Python API server (FastAPI recommended) + static file serving.
- Frontend: lightweight HTML/CSS/JS app.
- Interactive view: embed Mol* on detail page/modal.

### 3.4 Recommended Local-First Server Shape
- Keep the web layer thin and Python-native.
- Use FastAPI only as a read-mostly API and static file server.
- Reuse existing Python gallery helpers for indexing, SQL querying, and thumbnail generation.
- Keep SQLite as the only datastore for local/modest datasets.
- Avoid an ORM initially; direct SQL or thin helper functions are enough for the current schema.
- Avoid a separate SPA frontend initially; use static HTML/CSS/JS with a small amount of browser-side fetch logic.
- Keep thumbnail generation offline/asynchronous from browsing; the web app should read cached PNGs, not render them on demand.

## 4. Filesystem + Data Layout

## 4.1 Suggested Layout
- `data/runs/<run_id>/`:
  - volumizer output files (`.annotated.cif`, `.annotation.json`, `run.summary.json`)
- `data/gallery.db`:
  - SQLite index for filters and render metadata
- `data/renders/<run_id>/<source_label>/`:
  - `x.png`, `y.png`, `z.png`

## 4.2 SQLite Schema (initial)

### `runs`
- `run_id` (pk)
- `created_at`
- `source_summary_path`
- `resolution`
- `assembly_policy`

### `structures`
- `structure_id` (pk)
- `run_id` (fk)
- `source_label`
- `pdb_id` (nullable)
- `input_path`
- `annotated_cif_path`
- `annotation_json_path`
- `num_chains`
- `num_residues`
- `num_sequence_unique_chains`

### `volumes`
- `volume_id` (pk)
- `structure_id` (fk)
- `kind` (`pore|pocket|cavity|hub`)
- `rank_in_kind`
- `volume_a3`
- `length_a`
- `min_diameter_a`
- `max_diameter_a`
- `centroid_x`, `centroid_y`, `centroid_z`

### `structure_aggregates`
- `structure_id` (pk/fk)
- `num_pores`
- `largest_pore_volume_a3`
- `largest_pore_length_a`
- `largest_pore_min_diameter_a`
- `largest_pore_max_diameter_a`

### `renders`
- `structure_id` (pk/fk)
- `x_png_path`, `y_png_path`, `z_png_path`
- `render_style_hash`
- `render_status` (`pending|done|failed`)
- `render_error` (nullable)
- `updated_at`

## 5. Filter Strategy

## 5.1 Pore Filters
- Structure-level and pore-level views:
  - `min/max largest_pore_volume_a3`
  - `min/max largest_pore_length_a`
  - `min/max largest_pore_min_diameter_a`
  - `min/max largest_pore_max_diameter_a`
  - `min_pore_count`

## 5.2 Structural Filters
- `min/max num_chains`
- `min/max num_residues`
- `min/max num_sequence_unique_chains`

## 5.3 Query Execution
- Use SQL WHERE clauses + indexes:
  - index on `structures(num_chains, num_residues, num_sequence_unique_chains)`
  - index on `structure_aggregates` pore metric columns
- API supports pagination + sorting (`largest_pore_volume`, `num_residues`, etc.).

## 6. Render Pipeline Details (Step 3)

## 6.1 Why Mol* for renders
- Same visual language as interactive viewer.
- No mismatch between gallery and detail view.
- Single representation config to maintain.

## 6.2 Render job flow
1. Determine stale/missing entries from `renders` table.
2. For each structure:
   - Load annotated CIF in Mol*.
   - Apply consistent style (cartoon/surface + highlighted pore pseudo-atoms).
   - Set camera to +X, +Y, +Z canonical views.
   - Capture PNGs at fixed size (e.g. `320x240`).
3. Write files and update `renders` row.

## 6.3 Determinism
- Keep fixed:
  - background color
  - representation presets
  - camera distance and clipping
  - image dimensions
- Maintain `render_style_hash` to invalidate old thumbnails when style changes.

## 7. Interactive Viewer Details (Step 4)

1. Detail page receives `structure_id`.
2. Loads corresponding annotated CIF path from API.
3. Initializes Mol* viewer in browser.
4. Optional UX:
   - toggle pseudo-atoms on/off
   - focus on largest pore centroid
   - switch between `x/y/z` presets and free camera

## 8. API Plan (FastAPI)

## 8.1 Endpoints
- `GET /api/runs`
- `GET /api/hits` with filter params and pagination
- `GET /api/hits/{structure_id}`
- `GET /api/hits/{structure_id}/viewer-data`
- `POST /api/renders/rebuild` (optional admin/local-only)

## 8.2 Example `/api/hits` params
- `pore_volume_min`, `pore_volume_max`
- `pore_length_min`, `pore_length_max`
- `pore_dmin_min`, `pore_dmin_max`
- `pore_dmax_min`, `pore_dmax_max`
- `chains_min`, `chains_max`
- `residues_min`, `residues_max`
- `seq_unique_chains_min`, `seq_unique_chains_max`
- `limit`, `offset`, `sort_by`, `sort_dir`

## 9. Frontend Plan

## 9.1 Gallery page
- Left filter panel.
- Hit grid cards (3 thumbnails + summary stats).
- Pagination and sort controls.
- Card click opens detail modal/page.

## 9.2 Detail view
- Mol* canvas.
- Metadata panel with pore + structure stats.
- Links back to source files (`annotated.cif`, `annotation.json`).

## 10. Implementation Phases

## Phase A: Data + Filters (Python/Rust-heavy)
- Add enrichment/index script:
  - ingest summary + annotation JSON + CIF
  - compute structure metrics + pore metrics
  - write SQLite
- Deliverable: queryable DB + CLI query smoke test.

## Phase B: Thumbnail Rendering
- Add Node renderer (Mol* + Playwright/Puppeteer).
- Add render queue logic and cached outputs.
- Current implementation: `scripts/render_gallery_thumbnails.py` (queue) + `scripts/molstar_render_single.mjs` (Mol* worker).
- Runtime prerequisites: `playwright` package and Chromium browser (`npm install --save-dev playwright` + `npx playwright install chromium`).
- Deliverable: 3 PNGs per hit with deterministic style.

## Phase C: Local API + Gallery UI
- FastAPI endpoints + static frontend.
- Filter controls + paged gallery.
- Initial implementation now exists in `volumizer/web/app.py` with static assets under `volumizer/web/static/`.
- Deliverable: local website showing filtered hits + thumbnails.

## Phase D: Mol* Interactive View
- Integrate Mol* viewer on detail page.
- Add basic selection/toggle controls.
- Deliverable: click hit -> interactive molecular view.

## 11. Testing Plan

- Unit tests:
  - metric extraction from annotation JSON/CIF
  - SQL filter query correctness
- Integration tests:
  - index build from a sample run
  - API filter endpoints return expected subsets
- Render smoke tests:
  - generated files exist and non-zero size
- UI smoke tests:
  - filter -> card count updates
  - card click opens Mol* detail view

## 12. Risks and Mitigations

1. Render throughput for very large runs.
- Mitigation: background queue, incremental rendering, cache by style hash.

2. Metric ambiguity for multi-pore structures.
- Mitigation: store both aggregate metrics and pore-level table.

3. Mol* memory usage in headless batch mode.
- Mitigation: process in bounded batches, restart browser context periodically.

4. Data drift when rerunning analysis.
- Mitigation: run/version IDs and immutable run folders.

## 13. Acceptance Criteria (MVP)

- User can point the app at a `run.summary.json` and build index.
- User can filter by all requested pore + structural parameters.
- Gallery shows 3 cached thumbnails per hit.
- Clicking a hit opens Mol* interactive viewer for that structure.
- End-to-end works locally with documented commands.

## 14. Immediate Next Actions

1. Install Node render dependencies (`npm install --save-dev playwright` + `npx playwright install chromium`).
2. Run `scripts/render_gallery_thumbnails.py` on a real indexed run and validate generated `x/y/z` PNG quality.
3. Add minimal FastAPI endpoints (`/api/runs`, `/api/hits`, `/api/hits/{structure_id}`) for Phase C.
4. Scaffold gallery page with filter controls and thumbnail grid wired to `/api/hits`.
