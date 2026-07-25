# Implementation Plan: Upload-and-Analyze Landing Page

## Goal

Replace the gallery-only entry point with a **landing page** offering two paths:

1. **Upload a PDB** → run the volumizer analysis locally → show the annotated
   result directly in the interactive viewer → structure is added to the gallery.
2. **Browse gallery** → the existing gallery experience.

## Confirmed design decisions

- **Result UX**: show the live interactive Mol* 3D view immediately on completion;
  render the PNG gallery thumbnails lazily in the background.
- **Params**: the upload form exposes **resolution** + **assembly policy**; everything
  else uses existing defaults.
- **Grouping**: uploaded structures share a single `uploads` run in the gallery DB.

---

## Current architecture (reference)

- **Backend**: `volumizer/web/app.py` — FastAPI app, currently **read-only**. Serves the
  SPA (`static/index.html`, `static/app.js`, `static/styles.css`) and JSON APIs over the
  SQLite gallery DB (`data/gallery.db`).
- **DB helpers**: `volumizer/web/db.py` (read-only queries).
- **Indexing**: `volumizer/gallery_index.py` — `build_gallery_index(summary_path, db_path,
  run_id, replace_run, ...)` is **incremental and per-run** (`CREATE TABLE IF NOT EXISTS`,
  inserts one run into an existing DB). Schema of interest: `runs`, `structures`,
  `structure_aggregates`, `volumes`, `renders`. The `renders` row carries
  `render_status` (`pending` / `done` / `failed`).
- **Rendering**: `volumizer/gallery_render.py` — `render_gallery_thumbnails(...)` selects
  DB rows whose `render_status = 'pending'` (optionally `failed`) and produces PNGs via
  Playwright + Mol*. **Running it against the DB naturally picks up any new pending row.**
- **Analysis**: `volumizer.cli.analyze_structure_file(source_label, input_path, output_dir,
  min_voxels, min_volume, overwrite, assembly_policy=..., max_residues=None,
  include_hubs=False, enable_necked_pocket_cavity=True)` produces
  `<label>.annotated.cif` + `<label>.annotation.json` in `output_dir`. The isolated
  subprocess wrapper is `volumizer/_analysis_worker.py`.
- **Viewer**: the detail view already renders live 3D from the annotated CIF + annotation
  JSON via `GET /api/hits/{id}/viewer-data` — **independent of the PNG thumbnails**.

### Key enabling facts

1. `build_gallery_index` can append to an existing DB; the `gallery` bash script only wipes
   the DB (`rm -f`) for a *full* rebuild. The library call does not require a wipe.
2. `gallery_render` renders exactly the rows marked `pending`, so lazy background
   thumbnailing is just "insert a pending render row, then invoke the renderer."
3. The live viewer needs no PNGs, so the user sees the result the moment analysis +
   indexing finish.

### The one gotcha

The shared `uploads` run does **not** map onto `build_gallery_index(..., replace_run=True)`:
`replace_run` **deletes the run (cascading to `renders`)** before re-inserting from a
summary. Per-upload use would wipe prior uploads' rows and force re-rendering all their
thumbnails. **Resolution:** add an incremental "index one structure into an existing run"
helper that upserts a single structure without deleting the run.

---

## Work breakdown

Three independently testable stages. Build and verify in order.

Progress:
- **Stage 1**: complete. `index_single_structure(...)` now shares the batch indexing
  insert path through `_index_one(...)`; focused gallery index and web regression tests pass.
- **Stage 2**: not started.
- **Stage 3**: not started.

### Stage 1 — Incremental single-structure indexing

**File**: `volumizer/gallery_index.py`

Add a helper that upserts one analyzed structure into an existing `uploads` run:

```
def index_single_structure(
    db_path: Path,
    run_id: str,               # "uploads"
    source_label: str,
    input_path: Path,
    annotated_cif_path: Path,
    annotation_json_path: Path,
    resolution: float,
    assembly_policy: str,
) -> int:                      # returns structure_id
```

Behavior:
- `connection.executescript(_SCHEMA_SQL)` (idempotent) and `PRAGMA foreign_keys = ON`.
- Ensure the `uploads` run row exists (`INSERT OR IGNORE` into `runs` with `created_at`,
  a synthetic `source_summary_path` such as `data/runs/uploads/run.summary.json`,
  `resolution`, `assembly_policy`). Note per-run `resolution` is stored once at run
  creation; document that uploads share the first upload's run-level resolution (the
  per-structure viewer resolution comes from `run_meta.resolution`). If uploads must
  support mixed resolutions later, revisit run-level storage.
- Upsert the structure: delete any existing `(run_id, source_label)` structure row (cascade
  clears its aggregates/volumes/renders), then insert fresh. Reuse the existing metric
  computation (`_compute_structure_metrics`, `_normalize_volume_rows`, `_kind_aggregates`,
  `_infer_pdb_id`) already used by `build_gallery_index` — factor the per-structure insert
  body of `build_gallery_index` into a shared internal function to avoid duplication.
- Insert the `renders` row with `render_status = 'pending'`, null PNG paths,
  `updated_at = now`.
- Return the new `structure_id`.

Refactor note: extract the loop body of `build_gallery_index` (structure + aliases +
volumes + aggregates + renders insert) into `_index_one(connection, entry, ...)` so both
the batch path and `index_single_structure` call it. For upload indexing, construct a
small `result`-like dict (`source`, `input_path`, `structure_output`,
`annotation_output`) and pass that through the same internal helper rather than creating a
second structure-insert data path.

**Tests** (`tests/`): new test that builds a temp DB, calls `index_single_structure` twice
with different labels, and asserts: both structures present under run `uploads`, aggregates
+ volumes populated, render rows `pending`, and re-indexing the same label replaces rather
than duplicates.

### Stage 2 — Backend upload + analysis job API

**Files**: `volumizer/web/app.py`, new `volumizer/web/jobs.py`, `volumizer/web/db.py`.

**Job registry** (`volumizer/web/jobs.py`):
- In-process registry (`dict[str, Job]`) guarded by a lock. Single-user local app, so a
  single background worker thread + queue serializes jobs (avoids SQLite write contention).
- `Job` fields: `job_id`, `status` (`queued` / `running` / `indexing` / `done` /
  `error`), `message`, `structure_id` (on success), `error` (on failure),
  `created_at`, `updated_at`, and `thumbnail_status` (`pending` / `rendering` / `done` /
  `failed` / `skipped`) so the viewer can be shown as soon as indexing completes.
- `submit_analysis(...)` enqueues; a worker runs the pipeline. Design the registry with
  injectable analysis and thumbnail-render functions so FastAPI tests can avoid slow
  real analysis and Playwright rendering.

**Pipeline orchestration** (in the worker):
1. Save the uploaded bytes to `data/runs/uploads/incoming/<job_id>.<ext>` (accept `.pdb`,
   `.cif`, `.mmcif`; validate extension + non-empty + size cap).
2. `status = running`: run analysis via the `_analysis_worker` subprocess (mirrors the CLI's
   isolation) or a direct `analyze_structure_file` call in the worker thread. Prefer the
   subprocess for crash isolation and to reuse existing arg plumbing. Output dir:
   `data/runs/uploads/<source_label>/`. `source_label` derived from the filename (sanitized,
   deduped against existing labels in the run).
3. On analysis success: `status = indexing`; call
   `gallery_index.index_single_structure(...)` → `structure_id`.
4. `status = done`, set `structure_id` immediately after indexing so the frontend can route
   to the live detail viewer without waiting for PNG thumbnails.
5. Separately set `thumbnail_status = rendering` and invoke `gallery_render` against the DB
   in the background (subprocess `python -m volumizer.gallery_render --db <db>` or the
   library function). This picks up the new `pending` row. Thumbnail failure is **non-fatal**
   and updates only `thumbnail_status`; the indexed structure remains usable.
   Any exception → `status = error` with a user-safe message; keep the traceback server-side.

**Endpoints** (`app.py`):
- `POST /api/analyze` — `multipart/form-data`: `file`, `resolution` (float), `assembly_policy`
  (enum from `pdb.VALID_ASSEMBLY_POLICIES`). Validates inputs, enqueues a job, returns
  `{ "job_id": ... , "status": "queued" }`. Reject if analysis is disabled (see config note).
- `GET /api/analyze/{job_id}` — returns the job status payload; `structure_id` present when
  `done` so the frontend can route to `/api/hits/{id}`.
- Optionally `GET /api/analyze` — current/queued job for basic single-flight UX.

**Config / safety notes**:
- Add an env flag (e.g. `VOLUMIZER_ENABLE_UPLOAD`, default on for local) so the read-only
  deployment can disable writes. Surface it in `/api/health` (e.g. `upload_enabled`,
  `renderer_available`) so the frontend can hide the upload UI when unavailable.
- Enforce an upload size limit and a per-analysis `max_residues` guard (reuse
  `PostAssemblyResidueLimitExceeded` handling) to bound runtime.
- Confirm the FastAPI multipart dependency (`python-multipart`) is available; add it to the
  project dependencies if the first endpoint test reports it missing.
- `data/runs/uploads/` should be created on demand.

### Stage 3 — Frontend landing view + upload flow

**Files**: `volumizer/web/static/index.html`, `static/app.js`, `static/styles.css`.

- **View switching**: introduce a lightweight client-side view state
  (`landing` / `gallery` / `detail`) in `app.js`. Landing is the default. "Browse gallery"
  reveals the existing gallery workspace; the detail view is the existing viewer. No new
  HTML pages — toggle top-level sections in the existing shell. Optionally reflect state in
  the URL hash (`#/`, `#/gallery`, `#/hit/{id}`) for reload/back support.
- **Landing markup** (`index.html`): a hero with two primary actions (Upload / Browse) and
  an upload form: file input, `resolution` number input (default from health/config or a
  sensible constant), `assembly_policy` select (options from a small static list mirroring
  `VALID_ASSEMBLY_POLICIES`), submit button. Hide the form if `upload_enabled` is false.
- **Upload flow** (`app.js`):
  1. On submit, `POST /api/analyze` with `FormData`; disable the form (single-flight).
  2. Show a progress panel driven by polling `GET /api/analyze/{job_id}` (~1s interval),
     surfacing `status` transitions (`queued → running → indexing → done`) and, after
     `done`, showing thumbnail progress via `thumbnail_status` without blocking the viewer.
  3. On `done`, route to the detail view for `structure_id` (reuse the existing hit-detail
     loader that calls `/api/hits/{id}` + `/api/hits/{id}/viewer-data`) so the live 3D shows
     immediately. The gallery card thumbnail fills in once background rendering completes.
  4. On `error`, show the message and re-enable the form.
- **Styling** (`styles.css`): landing/hero + upload form + progress panel styles consistent
  with the existing masthead/pill visual language.

---

## Files touched (summary)

| File | Change |
|------|--------|
| `volumizer/gallery_index.py` | Add `index_single_structure`; extract shared `_index_one` insert helper |
| `volumizer/web/jobs.py` | **New** — in-process job registry + serialized worker + pipeline orchestration |
| `volumizer/web/app.py` | Add `POST /api/analyze`, `GET /api/analyze/{job_id}`; extend `/api/health`; upload-enabled guard |
| `volumizer/web/db.py` | (If needed) helper to fetch existing `uploads` source_labels for dedupe |
| `volumizer/web/static/index.html` | Landing view + upload form markup |
| `volumizer/web/static/app.js` | View switching, upload submit, job polling, route-to-detail |
| `volumizer/web/static/styles.css` | Landing/hero/upload/progress styles |
| `tests/` | Unit test for `index_single_structure`; API test for the analyze endpoints (mock the analysis) |

---

## Testing & verification

- **Stage 1**: pytest for `index_single_structure` (insert, re-insert idempotency, aggregates,
  pending render rows).
- **Stage 2**: FastAPI `TestClient` test that posts a small fixture PDB (or mocks the analysis
  worker) and polls the job to `done`, asserting a `structure_id` and a queryable
  `/api/hits/{id}`. Mock analysis and thumbnail rendering through the job registry's
  injectable test seams. Assert `error` handling for a bad upload.
- **Stage 3**: manual run via the `gallery` script / `uvicorn`, upload a known small PDB
  (e.g. a `tests/pdbs` fixture), confirm the live viewer appears on completion and the
  gallery card thumbnail fills in shortly after.
- **Regression**: existing gallery browsing, filters, and detail viewer remain unchanged;
  read-only deployments with `upload_enabled = false` hide the upload UI cleanly.

## Open questions / follow-ups (non-blocking)

- Run-level vs per-structure resolution for the shared `uploads` run (see Stage 1 note) —
  fine for now since uploads will typically share a default; revisit if mixed resolutions
  are needed.
- Cleanup policy for `data/runs/uploads/incoming/` temp files.
- Whether to also wire the upload path into the `gallery` bash orchestration or keep it
  purely server-driven (recommended: server-driven, the script stays a batch tool).
