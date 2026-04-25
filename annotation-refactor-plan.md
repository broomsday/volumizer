## Annotation Refactor Plan

### Goal

Remove the old "display-only" annotation concept and make the relevant code paths
change the true final annotation written to JSON, CIF, gallery index data, and
viewer inputs.

### Confirmed Current State

- `--merge-mouth-gap-voxels` is already applied during raw component
  classification, before JSON/CIF annotation output is built.
- The current display-only relabel is the necked-pocket override applied after
  the raw annotation dataframe is built.
- The viewer does not persist an explicit "mouth count". It renders "mouth"
  voxels from `surface_indices` via `B_iso_or_equiv > 0`.
- `7ODW` at `3.0 A` currently classifies incorrectly for the intended outcome:
  it yields 4 pores, 10 hubs, 123 pockets, and 0 cavities, with the largest row
  being a `pore` of volume `3043440.0`.
- Two of the smaller/larger `7ODW` pore rows are currently called `pore` despite
  having a single connected direct-mouth component. Those are produced by the
  directional fallback in the raw classifier rather than by a true two-mouth
  split, so `7ODW` is an appropriate regression target for the refactor.

### Concrete Edit Sequence

1. Add `7ODW` as a regression fixture and test case.
   - Add `tests/pdbs/7ODW.cif`.
   - Add an end-to-end test at `3.0 A` asserting:
     - largest annotated type is `cavity`
     - largest annotated volume is `3043440.0`
     - number of `pore` rows is `0`

2. Move neck-based relabeling out of the display layer and into the true
   annotation pipeline.
   - Remove the display-only post-processing helpers from
     `volumizer/volumizer.py`.
   - Replace them with a backend-neutral refinement pass over classified
     `VoxelGroup` objects before `Annotation` is constructed.

3. Generalize the existing neck heuristic in `volumizer/voxel.py`.
   - Reuse the existing neck-width/core-shell logic as the basis for final
     relabeling.
   - Add helpers for:
     - computing direct mouth components for one volume
     - checking whether a specific mouth is thin-necked
     - refining a raw `pocket` or `pore` into its final annotation type

4. Apply final refinement rules before annotation materialization.
   - Raw `pocket` with one thin neck -> final `cavity`
   - Raw `pore` with two thin necks -> final `cavity`
   - Raw `pore` with one thin neck -> final `pocket`
   - Otherwise unchanged

5. Rebucket, sort, and renumber after refinement.
   - After final type refinement, rebuild the `hubs`, `pores`, `pockets`, and
     `cavities` dictionaries.
   - Then apply existing filter/sort logic so ids, counts, largest-type
     summaries, dataframe rows, and output structure residue ids all reflect the
     final annotation, not the raw topology bucket.

6. Remove `display_type` plumbing from output generation.
   - Stop writing or preferring `display_type` in CLI payloads.
   - Stop passing display overrides into structure serialization.
   - Make CIF residue names come directly from final type.
   - Make gallery indexing use canonical `type` only for new outputs.
   - Keep read-side fallback in gallery indexing only if needed for old saved
     outputs.

7. Keep CLI compatibility, but change semantics and help text.
   - Keep `--no-necked-pocket-cavity` for now so existing invocations still
     parse.
   - Change its help text and internal meaning so it disables final neck-based
     relabeling, not a display-only override.
   - Update any related cluster/worker help text and payload fields to match.

8. Update tests to match canonical final types.
   - Replace display-type tests with final-type tests in:
     - `tests/test_voxel.py`
     - `tests/test_cli.py`
     - `tests/test_pdb.py`
     - `tests/test_gallery_index.py`
   - Add direct unit tests for:
     - necked pocket -> cavity
     - pore with two thin mouths -> cavity
     - pore with one thin mouth -> pocket
     - raw pore from the single-mouth directional fallback being corrected by
       the final refinement pass

9. Preserve backend parity.
   - Because the new refinement pass happens after raw component classification,
     it should be shared by both Python and native backends.
   - Re-run existing parity coverage and add targeted assertions if needed so
     final annotation outputs remain backend-consistent.

### Important Clarification About `--merge-mouth-gap-voxels`

`--merge-mouth-gap-voxels` already changes raw classification before final
annotation is written. It can change whether a raw component is typed as a
`pore` or `pocket`.

What it does not currently provide is a persisted "final mouth group count" data
model. The only serialized mouth/surface representation is the set of shell
voxels stored in `surface_indices` and rendered through B-factors. That means
the refactor should not assume that viewer mouth rendering is already tracking
merged mouth groups as a first-class concept.

### Implementation Focus

The first implementation target should be a single backend-neutral refinement
pass that corrects final annotation types and makes `7ODW` pass before the
display-only cleanup is fully scrubbed from every downstream test.
