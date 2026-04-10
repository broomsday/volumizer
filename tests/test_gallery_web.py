import json
from pathlib import Path
import sqlite3
import subprocess
import textwrap

from fastapi.testclient import TestClient

from volumizer import gallery_index, gallery_render
from volumizer.paths import TEST_DIR
from volumizer.web.app import create_app


TEST_INPUT_PDB = TEST_DIR / "pdbs" / "cavity.pdb"


def _build_molstar_assets(tmp_path: Path) -> Path:
    asset_root = tmp_path / "molstar"
    asset_root.mkdir(parents=True, exist_ok=True)
    (asset_root / "molstar.js").write_text("window.molstar = { Viewer: {} };", encoding="utf-8")
    (asset_root / "molstar.css").write_text("body { background: #fff; }\n", encoding="utf-8")
    return asset_root


def _write_annotation(path: Path, volume: float) -> None:
    path.write_text(
        json.dumps(
            {
                "source": path.stem,
                "num_volumes": 1,
                "frac_alpha": 0.40,
                "frac_beta": 0.30,
                "frac_coil": 0.30,
                "volumes": [
                    {
                        "id": 0,
                        "type": "pore",
                        "volume": volume,
                        "x": 12.0,
                        "y": 5.0,
                        "z": 3.0,
                        "cross_section_circularity": 0.85,
                        "cross_section_uniformity": 0.72,
                    }
                ],
            }
        ),
        encoding="utf-8",
    )


def _build_web_fixture(tmp_path: Path) -> tuple[Path, str, dict[str, int]]:
    run_id = "web-fixture"
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    pdb_ids = {
        "hit-a": "4JPN",
        "hit-b": "1ABC",
    }

    results = []
    for source_label, volume in (("hit-a", 200.0), ("hit-b", 80.0)):
        annotation_path = run_dir / f"{source_label}.annotation.json"
        structure_output_path = run_dir / f"{source_label}.annotated.cif"
        _write_annotation(annotation_path, volume)
        structure_output_path.write_text("data_test\n#\n", encoding="utf-8")

        result_entry = {
            "source": source_label,
            "pdb_id": pdb_ids[source_label],
            "input_path": str(TEST_INPUT_PDB),
            "structure_output": str(structure_output_path),
            "annotation_output": str(annotation_path),
        }
        if source_label == "hit-a":
            result_entry["cluster_member_pdb_ids"] = ["4JPN", "4JPP"]
        results.append(result_entry)

    summary_path = run_dir / "run.summary.json"
    summary_path.write_text(
        json.dumps(
            {
                "config": {
                    "assembly_policy": "biological",
                    "resolution": 3.0,
                    "keep_non_protein": False,
                    "output_dir": str(run_dir),
                },
                "results": results,
                "errors": [],
                "skipped": [],
                "planned": [],
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    db_path = tmp_path / "gallery.db"
    gallery_index.build_gallery_index(
        summary_path=summary_path,
        db_path=db_path,
        run_id=run_id,
        replace_run=False,
        strict=True,
    )

    render_root = tmp_path / "renders"

    def fake_render(
        structure_path: Path,
        output_dir: Path,
        width: int,
        height: int,
        style: dict,
    ) -> None:
        assert structure_path.is_file()
        assert width > 0
        assert height > 0
        assert style["width"] == width
        assert style["height"] == height
        output_dir.mkdir(parents=True, exist_ok=True)
        for axis in ("x", "y", "z"):
            (output_dir / f"{axis}.png").write_bytes(f"png-{axis}".encode("utf-8"))

    gallery_render.render_gallery_thumbnails(
        db_path=db_path,
        render_root=render_root,
        run_id=run_id,
        render_fn=fake_render,
    )

    with sqlite3.connect(db_path) as connection:
        rows = connection.execute(
            "SELECT structure_id, source_label FROM structures ORDER BY source_label ASC"
        ).fetchall()

    structure_ids = {str(source_label): int(structure_id) for structure_id, source_label in rows}
    return db_path, run_id, structure_ids


def test_gallery_web_root_and_health(tmp_path: Path):
    db_path, _, _ = _build_web_fixture(tmp_path)
    asset_root = _build_molstar_assets(tmp_path)
    client = TestClient(create_app(db_path, molstar_asset_root=asset_root))

    root_response = client.get("/")
    assert root_response.status_code == 200
    assert "Volumizer Gallery" in root_response.text
    assert 'name="pdb_id_query"' in root_response.text

    health_response = client.get("/api/health")
    assert health_response.status_code == 200
    payload = health_response.json()
    assert payload["status"] == "ok"
    assert payload["db_exists"] is True
    assert payload["db"] == str(db_path.resolve())
    assert payload["molstar_assets_available"] is True
    assert payload["molstar_asset_root"] == str(asset_root.resolve())


def test_gallery_web_serves_local_molstar_assets(tmp_path: Path):
    db_path, _, _ = _build_web_fixture(tmp_path)
    asset_root = _build_molstar_assets(tmp_path)
    client = TestClient(create_app(db_path, molstar_asset_root=asset_root))

    js_response = client.get("/assets/molstar.js")
    assert js_response.status_code == 200
    assert "window.molstar" in js_response.text

    css_response = client.get("/assets/molstar.css")
    assert css_response.status_code == 200
    assert "background" in css_response.text


def test_gallery_web_lists_runs_and_hits(tmp_path: Path):
    db_path, run_id, structure_ids = _build_web_fixture(tmp_path)
    client = TestClient(create_app(db_path))

    runs_response = client.get("/api/runs")
    assert runs_response.status_code == 200
    runs_payload = runs_response.json()
    assert runs_payload["runs"][0]["run_id"] == run_id
    assert runs_payload["runs"][0]["structure_count"] == 2

    hits_response = client.get("/api/hits")
    assert hits_response.status_code == 200
    hits_payload = hits_response.json()
    assert hits_payload["total_count"] == 2
    assert len(hits_payload["rows"]) == 2

    hit_a_id = structure_ids["hit-a"]
    hit_a_row = next(row for row in hits_payload["rows"] if row["structure_id"] == hit_a_id)
    assert hit_a_row["urls"]["detail"] == f"/api/hits/{hit_a_id}"
    assert hit_a_row["urls"]["structure"] == f"/files/structure/{hit_a_id}"
    assert hit_a_row["urls"]["thumbnail_x"] == f"/files/render/{hit_a_id}/x.png"

    filtered_response = client.get("/api/hits?pore_volume_min=150")
    assert filtered_response.status_code == 200
    filtered_payload = filtered_response.json()
    assert filtered_payload["total_count"] == 1
    assert filtered_payload["rows"][0]["source_label"] == "hit-a"

    filtered_max_response = client.get("/api/hits?pore_volume_max=150")
    assert filtered_max_response.status_code == 200
    filtered_max_payload = filtered_max_response.json()
    assert filtered_max_payload["total_count"] == 1
    assert filtered_max_payload["rows"][0]["source_label"] == "hit-b"

    pdb_filtered_response = client.get("/api/hits?pdb_id_query=jp")
    assert pdb_filtered_response.status_code == 200
    pdb_filtered_payload = pdb_filtered_response.json()
    assert pdb_filtered_payload["total_count"] == 1
    assert pdb_filtered_payload["rows"][0]["source_label"] == "hit-a"

    alias_filtered_response = client.get("/api/hits?pdb_id_query=4jpp")
    assert alias_filtered_response.status_code == 200
    alias_filtered_payload = alias_filtered_response.json()
    assert alias_filtered_payload["total_count"] == 1
    assert alias_filtered_payload["rows"][0]["source_label"] == "hit-a"


def test_gallery_web_api_filters_by_protein_ranges(tmp_path: Path):
    db_path, _, _ = _build_web_fixture(tmp_path)
    with sqlite3.connect(db_path) as connection:
        connection.execute(
            "UPDATE structures SET num_chains = 1, num_residues = 100, num_sequence_unique_chains = 1 WHERE source_label = 'hit-a'"
        )
        connection.execute(
            "UPDATE structures SET num_chains = 4, num_residues = 300, num_sequence_unique_chains = 4 WHERE source_label = 'hit-b'"
        )
        connection.commit()

    client = TestClient(create_app(db_path))

    filtered_response = client.get(
        "/api/hits?chains_max=2&residues_max=150&seq_unique_chains_max=1"
    )
    assert filtered_response.status_code == 200
    filtered_payload = filtered_response.json()
    assert filtered_payload["total_count"] == 1
    assert filtered_payload["rows"][0]["source_label"] == "hit-a"


def test_gallery_web_api_filters_by_length_maximum(tmp_path: Path):
    db_path, _, _ = _build_web_fixture(tmp_path)
    with sqlite3.connect(db_path) as connection:
        connection.execute(
            """
            UPDATE structure_aggregates
            SET largest_pore_length_a = 40
            WHERE structure_id = (
                SELECT structure_id FROM structures WHERE source_label = 'hit-a'
            )
            """
        )
        connection.execute(
            """
            UPDATE structure_aggregates
            SET largest_pore_length_a = 12
            WHERE structure_id = (
                SELECT structure_id FROM structures WHERE source_label = 'hit-b'
            )
            """
        )
        connection.commit()

    client = TestClient(create_app(db_path))

    filtered_response = client.get("/api/hits?pore_length_max=20")
    assert filtered_response.status_code == 200
    filtered_payload = filtered_response.json()
    assert filtered_payload["total_count"] == 1
    assert filtered_payload["rows"][0]["source_label"] == "hit-b"


def test_gallery_web_static_search_serializes_form_fields_generically():
    static_dir = Path(__file__).resolve().parents[1] / "volumizer" / "web" / "static"
    app_js = (static_dir / "app.js").read_text(encoding="utf-8")
    index_html = (static_dir / "index.html").read_text(encoding="utf-8")

    assert "for (const field of elements.filtersForm.elements)" in app_js
    assert "const numericNames = [" not in app_js
    assert 'name="seq_unique_chains_max"' in index_html
    assert 'name="chains_max"' in index_html
    assert 'name="residues_max"' in index_html
    assert 'name="pore_volume_max"' in index_html
    assert 'name="pore_length_max"' in index_html
    assert 'name="pocket_volume_max"' in index_html
    assert 'name="pocket_length_max"' in index_html
    assert 'name="cavity_volume_max"' in index_html
    assert 'name="cavity_length_max"' in index_html
    assert 'name="hub_volume_max"' in index_html
    assert 'name="hub_length_max"' in index_html
    assert '<option value="">None</option>' in index_html
    assert 'name="seq_unique_chains_max" value="3"' not in index_html
    assert 'name="frac_coil_max" value="0.6"' not in index_html
    assert 'name="pore_circularity_min" value="0.5"' not in index_html
    assert 'name="pore_uniformity_min" value="0.5"' not in index_html
    assert '<option value="">Default</option>' in index_html
    assert '<input type="checkbox" data-metric="num_chains" checked /> Chains' in index_html
    assert '<input type="checkbox" data-metric="num_residues" checked /> Residues' in index_html
    assert '<input type="checkbox" data-metric="num_sequence_unique_chains" checked /> Unique chains' in index_html
    assert '<input type="checkbox" data-metric="frac_coil" checked /> Coil' not in index_html
    assert 'href="/assets/molstar.css"' in index_html
    assert 'src="/assets/molstar.js"' in index_html
    assert "unpkg.com/molstar" not in index_html


def test_gallery_web_static_filter_presets_restore_active_selection():
    static_dir = Path(__file__).resolve().parents[1] / "volumizer" / "web" / "static"
    app_js = (static_dir / "app.js").read_text(encoding="utf-8")

    assert "const ACTIVE_FILTER_PRESET_KEY = 'volumizer_active_filter_preset';" in app_js
    assert "function restoreActiveFilterPreset()" in app_js
    assert "setActiveFilterPresetName(name);" in app_js
    assert "restoreActiveFilterPreset();" in app_js
    assert "event.preventDefault();" in app_js
    assert "const restoredActivePreset = restoreActiveFilterPreset();" in app_js


def test_gallery_web_static_sparse_filter_presets_reset_previous_values():
    static_dir = Path(__file__).resolve().parents[1] / "volumizer" / "web" / "static"
    app_js_path = static_dir / "app.js"

    node_script = textwrap.dedent(
        f"""
        const fs = require('node:fs');
        const vm = require('node:vm');

        const source = fs.readFileSync({str(app_js_path)!r}, 'utf8')
          + '\\n;globalThis.__galleryTest = {{ applyFilterState, state }};';

        function createField(name, value, defaultValue = '') {{
          return {{
            name,
            value,
            defaultValue,
            type: 'text',
            checked: false,
            addEventListener() {{}},
          }};
        }}

        const chainsMax = createField('chains_max', '99', '');
        const seqUniqueChainsMax = createField('seq_unique_chains_max', '1', '3');
        const sortBy = createField('sort_by', 'num_chains', 'largest_pore_volume_a3');
        const sortDir = createField('sort_dir', 'asc', 'desc');
        const limit = createField('limit', '12', '24');

        const filtersForm = {{
          elements: [chainsMax, seqUniqueChainsMax, sortBy, sortDir, limit],
          reset() {{
            for (const el of this.elements) {{
              el.value = el.defaultValue;
            }}
          }},
          addEventListener() {{}},
        }};

        const pdbIdQueryInput = {{ value: '1ABC', addEventListener() {{}} }};
        const displayPropsPopup = {{ querySelectorAll() {{ return []; }}, hidden: true }};
        const detailDialog = {{
          open: false,
          addEventListener() {{}},
          getBoundingClientRect() {{
            return {{ top: 0, left: 0, width: 0, height: 0 }};
          }},
          showModal() {{}},
          close() {{}},
        }};

        function stubElement() {{
          return {{
            value: '',
            innerHTML: '',
            textContent: '',
            hidden: false,
            disabled: false,
            open: false,
            options: [],
            append() {{}},
            appendChild() {{}},
            addEventListener() {{}},
            querySelectorAll() {{ return []; }},
            getBoundingClientRect() {{
              return {{ top: 0, left: 0, width: 0, height: 0 }};
            }},
            showModal() {{}},
            close() {{}},
          }};
        }}

        const elementIds = [
          'filters-form',
          'pdb-id-query',
          'run-id',
          'limit',
          'sort-by',
          'sort-dir',
          'prev-page',
          'next-page',
          'reset-filters',
          'loading-state',
          'empty-state',
          'results-grid',
          'db-pill',
          'result-pill',
          'page-status',
          'detail-dialog',
          'close-detail',
          'detail-title',
          'detail-links',
          'detail-metrics',
          'detail-volumes',
          'viewer-host',
          'display-props-toggle',
          'display-props-popup',
          'filter-preset-select',
          'filter-preset-load',
          'filter-preset-save',
          'filter-preset-delete',
          'display-preset-select',
          'display-preset-load',
          'display-preset-save',
          'display-preset-delete',
        ];

        const elementsById = Object.fromEntries(
          elementIds.map((id) => [id, stubElement()]),
        );
        elementsById['filters-form'] = filtersForm;
        elementsById['pdb-id-query'] = pdbIdQueryInput;
        elementsById['display-props-popup'] = displayPropsPopup;
        elementsById['detail-dialog'] = detailDialog;
        elementsById['limit'] = limit;
        elementsById['sort-by'] = sortBy;
        elementsById['sort-dir'] = sortDir;

        const document = {{
          getElementById(id) {{
            return elementsById[id] || stubElement();
          }},
          createElement() {{
            return stubElement();
          }},
        }};

        const context = {{
          console,
          document,
          window: {{ addEventListener() {{}}, molstar: null }},
          localStorage: {{ getItem() {{ return null; }}, setItem() {{}}, removeItem() {{}} }},
          URLSearchParams,
          fetch: async () => ({{ ok: true, json: async () => ({{}}) }}),
          confirm: () => true,
          prompt: () => '',
          Set,
          Number,
          String,
          JSON,
          Promise,
          Object,
          Array,
          Math,
        }};

        vm.createContext(context);
        vm.runInContext(source, context);

        context.__galleryTest.applyFilterState({{
          seq_unique_chains_max: '3',
          sort_by: 'largest_pore_volume_a3',
          sort_dir: 'desc',
          limit: '24',
        }});

        console.log(JSON.stringify({{
          chains_max: chainsMax.value,
          seq_unique_chains_max: seqUniqueChainsMax.value,
          sort_by: sortBy.value,
          sort_dir: sortDir.value,
          limit: limit.value,
          pdb_id_query: pdbIdQueryInput.value,
          current_limit: context.__galleryTest.state.currentLimit,
          current_offset: context.__galleryTest.state.currentOffset,
        }}));
        """
    )

    result = subprocess.run(
        ["node", "-e", node_script],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(result.stdout)

    assert payload["chains_max"] == ""
    assert payload["seq_unique_chains_max"] == "3"
    assert payload["sort_by"] == "largest_pore_volume_a3"
    assert payload["sort_dir"] == "desc"
    assert payload["limit"] == "24"
    assert payload["pdb_id_query"] == ""
    assert payload["current_limit"] == 24
    assert payload["current_offset"] == 0


def test_gallery_web_static_none_filter_preset_clears_restored_values():
    static_dir = Path(__file__).resolve().parents[1] / "volumizer" / "web" / "static"
    app_js_path = static_dir / "app.js"

    node_script = textwrap.dedent(
        f"""
        const fs = require('node:fs');
        const vm = require('node:vm');

        const source = fs.readFileSync({str(app_js_path)!r}, 'utf8')
          + '\\n;globalThis.__galleryTest = {{ applyFilterState, state }};';

        function createField(name, value, defaultValue = '') {{
          return {{
            name,
            value,
            defaultValue,
            type: 'text',
            checked: false,
            addEventListener() {{}},
          }};
        }}

        const runId = createField('run_id', 'legacy-run', '');
        const seqUniqueChainsMax = createField('seq_unique_chains_max', '3', '');
        const fracCoilMax = createField('frac_coil_max', '0.6', '');
        const poreVolumeMin = createField('pore_volume_min', '10000', '');
        const poreDminMin = createField('pore_dmin_min', '35', '');
        const poreCircularityMin = createField('pore_circularity_min', '0.8', '');
        const poreUniformityMin = createField('pore_uniformity_min', '0.5', '');
        const sortBy = createField('sort_by', 'largest_pore_volume_a3', 'largest_pore_volume_a3');
        const sortDir = createField('sort_dir', 'desc', 'desc');
        const limit = createField('limit', '24', '24');

        const filtersForm = {{
          elements: [
            runId,
            seqUniqueChainsMax,
            fracCoilMax,
            poreVolumeMin,
            poreDminMin,
            poreCircularityMin,
            poreUniformityMin,
            sortBy,
            sortDir,
            limit,
          ],
          reset() {{
            runId.value = 'legacy-run';
            seqUniqueChainsMax.value = '3';
            fracCoilMax.value = '0.6';
            poreVolumeMin.value = '10000';
            poreDminMin.value = '35';
            poreCircularityMin.value = '0.8';
            poreUniformityMin.value = '0.5';
            sortBy.value = 'largest_pore_volume_a3';
            sortDir.value = 'desc';
            limit.value = '24';
          }},
          addEventListener() {{}},
        }};

        const pdbIdQueryInput = {{ value: '1ABC', addEventListener() {{}} }};
        const displayPropsPopup = {{ querySelectorAll() {{ return []; }}, hidden: true }};
        const detailDialog = {{
          open: false,
          addEventListener() {{}},
          getBoundingClientRect() {{
            return {{ top: 0, left: 0, width: 0, height: 0 }};
          }},
          showModal() {{}},
          close() {{}},
        }};

        function stubElement() {{
          return {{
            value: '',
            innerHTML: '',
            textContent: '',
            hidden: false,
            disabled: false,
            open: false,
            options: [],
            append() {{}},
            appendChild() {{}},
            addEventListener() {{}},
            querySelectorAll() {{ return []; }},
            getBoundingClientRect() {{
              return {{ top: 0, left: 0, width: 0, height: 0 }};
            }},
            showModal() {{}},
            close() {{}},
          }};
        }}

        const elementIds = [
          'filters-form',
          'pdb-id-query',
          'run-id',
          'limit',
          'sort-by',
          'sort-dir',
          'prev-page',
          'next-page',
          'reset-filters',
          'loading-state',
          'empty-state',
          'results-grid',
          'db-pill',
          'result-pill',
          'page-status',
          'detail-dialog',
          'close-detail',
          'detail-title',
          'detail-links',
          'detail-metrics',
          'detail-volumes',
          'viewer-host',
          'display-props-toggle',
          'display-props-popup',
          'filter-preset-select',
          'filter-preset-load',
          'filter-preset-save',
          'filter-preset-delete',
          'display-preset-select',
          'display-preset-load',
          'display-preset-save',
          'display-preset-delete',
        ];

        const elementsById = Object.fromEntries(
          elementIds.map((id) => [id, stubElement()]),
        );
        elementsById['filters-form'] = filtersForm;
        elementsById['pdb-id-query'] = pdbIdQueryInput;
        elementsById['display-props-popup'] = displayPropsPopup;
        elementsById['detail-dialog'] = detailDialog;
        elementsById['run-id'] = runId;
        elementsById['limit'] = limit;
        elementsById['sort-by'] = sortBy;
        elementsById['sort-dir'] = sortDir;

        const document = {{
          getElementById(id) {{
            return elementsById[id] || stubElement();
          }},
          createElement() {{
            return stubElement();
          }},
        }};

        const context = {{
          console,
          document,
          window: {{ addEventListener() {{}}, molstar: null }},
          localStorage: {{ getItem() {{ return null; }}, setItem() {{}}, removeItem() {{}} }},
          URLSearchParams,
          fetch: async () => ({{ ok: true, json: async () => ({{}}) }}),
          confirm: () => true,
          prompt: () => '',
          alert: () => {{}},
          Set,
          Number,
          String,
          JSON,
          Promise,
          Object,
          Array,
          Math,
        }};

        vm.createContext(context);
        vm.runInContext(source, context);
        context.__galleryTest.applyFilterState({{}});

        console.log(JSON.stringify({{
          run_id: runId.value,
          seq_unique_chains_max: seqUniqueChainsMax.value,
          frac_coil_max: fracCoilMax.value,
          pore_volume_min: poreVolumeMin.value,
          pore_dmin_min: poreDminMin.value,
          pore_circularity_min: poreCircularityMin.value,
          pore_uniformity_min: poreUniformityMin.value,
          sort_by: sortBy.value,
          sort_dir: sortDir.value,
          limit: limit.value,
          pdb_id_query: pdbIdQueryInput.value,
          current_limit: context.__galleryTest.state.currentLimit,
          current_offset: context.__galleryTest.state.currentOffset,
        }}));
        """
    )

    result = subprocess.run(
        ["node", "-e", node_script],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(result.stdout)

    assert payload["run_id"] == ""
    assert payload["seq_unique_chains_max"] == ""
    assert payload["frac_coil_max"] == ""
    assert payload["pore_volume_min"] == ""
    assert payload["pore_dmin_min"] == ""
    assert payload["pore_circularity_min"] == ""
    assert payload["pore_uniformity_min"] == ""
    assert payload["sort_by"] == "largest_pore_volume_a3"
    assert payload["sort_dir"] == "desc"
    assert payload["limit"] == "24"
    assert payload["pdb_id_query"] == ""
    assert payload["current_limit"] == 24
    assert payload["current_offset"] == 0


def test_gallery_web_static_builtin_filter_presets_exist():
    static_dir = Path(__file__).resolve().parents[1] / "volumizer" / "web" / "static"
    app_js_path = static_dir / "app.js"

    node_script = textwrap.dedent(
        f"""
        const fs = require('node:fs');
        const vm = require('node:vm');

        const source = fs.readFileSync({str(app_js_path)!r}, 'utf8')
          + '\\n;globalThis.__galleryTest = {{'
          + ' BUILTIN_FILTER_PRESET_PREFIX,'
          + ' BUILTIN_FILTER_PRESETS,'
          + ' getFilterPresetData,'
          + ' populateFilterPresetSelect,'
          + ' elements'
          + ' }};';

        function createField(name, value = '', defaultValue = '') {{
          return {{
            name,
            value,
            defaultValue,
            type: 'text',
            checked: false,
            addEventListener() {{}},
          }};
        }}

        function createSelect() {{
          return {{
            value: '',
            innerHTML: '',
            textContent: '',
            hidden: false,
            disabled: false,
            open: false,
            options: [],
            addEventListener() {{}},
            append(option) {{
              this.options.push(option);
            }},
            appendChild(option) {{
              this.options.push(option);
            }},
            querySelectorAll() {{ return []; }},
            getBoundingClientRect() {{
              return {{ top: 0, left: 0, width: 0, height: 0 }};
            }},
            showModal() {{}},
            close() {{}},
          }};
        }}

        function stubElement() {{
          return {{
            value: '',
            innerHTML: '',
            textContent: '',
            hidden: false,
            disabled: false,
            open: false,
            options: [],
            append() {{}},
            appendChild() {{}},
            addEventListener() {{}},
            querySelectorAll() {{ return []; }},
            getBoundingClientRect() {{
              return {{ top: 0, left: 0, width: 0, height: 0 }};
            }},
            showModal() {{}},
            close() {{}},
          }};
        }}

        const filtersForm = {{
          elements: [
            createField('seq_unique_chains_max'),
            createField('frac_coil_max'),
            createField('pore_volume_min'),
            createField('pore_dmin_min'),
            createField('pore_circularity_min'),
            createField('pore_uniformity_min'),
            createField('limit', '24', '24'),
          ],
          reset() {{}},
          addEventListener() {{}},
        }};

        const filterPresetSelect = createSelect();
        const elementIds = [
          'filters-form',
          'pdb-id-query',
          'run-id',
          'limit',
          'sort-by',
          'sort-dir',
          'prev-page',
          'next-page',
          'reset-filters',
          'loading-state',
          'empty-state',
          'results-grid',
          'db-pill',
          'result-pill',
          'page-status',
          'detail-dialog',
          'close-detail',
          'detail-title',
          'detail-links',
          'detail-metrics',
          'detail-volumes',
          'viewer-host',
          'display-props-toggle',
          'display-props-popup',
          'filter-preset-select',
          'filter-preset-load',
          'filter-preset-save',
          'filter-preset-delete',
          'display-preset-select',
          'display-preset-load',
          'display-preset-save',
          'display-preset-delete',
        ];

        const elementsById = Object.fromEntries(
          elementIds.map((id) => [id, stubElement()]),
        );
        elementsById['filters-form'] = filtersForm;
        elementsById['pdb-id-query'] = createField('pdb_id_query');
        elementsById['limit'] = filtersForm.elements[6];
        elementsById['sort-by'] = createField('sort_by', 'largest_pore_volume');
        elementsById['sort-dir'] = createField('sort_dir', 'desc');
        elementsById['detail-dialog'] = stubElement();
        elementsById['display-props-popup'] = stubElement();
        elementsById['filter-preset-select'] = filterPresetSelect;

        const document = {{
          getElementById(id) {{
            return elementsById[id] || stubElement();
          }},
          createElement() {{
            return {{
              value: '',
              textContent: '',
            }};
          }},
        }};

        const context = {{
          console,
          document,
          window: {{ addEventListener() {{}}, molstar: null }},
          localStorage: {{
            getItem(key) {{
              if (key === 'volumizer_filter_presets') {{
                return JSON.stringify({{ Custom: {{ seq_unique_chains_max: '9' }} }});
              }}
              return null;
            }},
            setItem() {{}},
            removeItem() {{}},
          }},
          URLSearchParams,
          fetch: async () => ({{ ok: true, json: async () => ({{}}) }}),
          confirm: () => true,
          prompt: () => '',
          alert: () => {{}},
          Set,
          Number,
          String,
          JSON,
          Promise,
          Object,
          Array,
          Math,
        }};

        vm.createContext(context);
        vm.runInContext(source, context);
        context.__galleryTest.populateFilterPresetSelect();

        const poresValue = (
          context.__galleryTest.BUILTIN_FILTER_PRESET_PREFIX + 'Pores'
        );
        const pocketsValue = (
          context.__galleryTest.BUILTIN_FILTER_PRESET_PREFIX + 'Pockets'
        );
        const cavitiesValue = (
          context.__galleryTest.BUILTIN_FILTER_PRESET_PREFIX + 'Cavities'
        );
        const porePreset = context.__galleryTest.getFilterPresetData(poresValue);
        const pocketPreset = context.__galleryTest.getFilterPresetData(pocketsValue);
        const cavityPreset = context.__galleryTest.getFilterPresetData(cavitiesValue);

        console.log(JSON.stringify({{
          option_labels: context.__galleryTest.elements.filterPresetSelect.options.map(
            (option) => option.textContent,
          ),
          pore_preset: porePreset,
          pocket_preset: pocketPreset,
          cavity_preset: cavityPreset,
        }}));
        """
    )

    result = subprocess.run(
        ["node", "-e", node_script],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(result.stdout)

    assert payload["option_labels"] == ["None", "Pores", "Pockets", "Cavities", "Custom"]
    assert payload["pore_preset"] == {
        "seq_unique_chains_max": "3",
        "frac_coil_max": "0.6",
        "pore_volume_min": "10000",
        "pore_dmin_min": "35",
        "pore_circularity_min": "0.8",
        "pore_uniformity_min": "0.5",
    }
    assert payload["pocket_preset"] == {
        "seq_unique_chains_max": "3",
        "frac_coil_max": "0.6",
    }
    assert payload["cavity_preset"] == {
        "seq_unique_chains_max": "3",
        "frac_coil_max": "0.6",
    }


def test_gallery_web_static_builtin_display_presets_exist():
    static_dir = Path(__file__).resolve().parents[1] / "volumizer" / "web" / "static"
    app_js_path = static_dir / "app.js"

    node_script = textwrap.dedent(
        f"""
        const fs = require('node:fs');
        const vm = require('node:vm');

        const source = fs.readFileSync({str(app_js_path)!r}, 'utf8')
          + '\\n;globalThis.__galleryTest = {{'
          + ' BUILTIN_DISPLAY_PRESET_PREFIX,'
          + ' getDisplayPresetData,'
          + ' populateDisplayPresetSelect,'
          + ' elements'
          + ' }};';

        function createField(name, value = '', defaultValue = '') {{
          return {{
            name,
            value,
            defaultValue,
            type: 'text',
            checked: false,
            addEventListener() {{}},
          }};
        }}

        function createSelect() {{
          return {{
            value: '',
            innerHTML: '',
            textContent: '',
            hidden: false,
            disabled: false,
            open: false,
            options: [],
            addEventListener() {{}},
            append(option) {{
              this.options.push(option);
            }},
            appendChild(option) {{
              this.options.push(option);
            }},
            querySelectorAll() {{ return []; }},
            getBoundingClientRect() {{
              return {{ top: 0, left: 0, width: 0, height: 0 }};
            }},
            showModal() {{}},
            close() {{}},
          }};
        }}

        function stubElement() {{
          return {{
            value: '',
            innerHTML: '',
            textContent: '',
            hidden: false,
            disabled: false,
            open: false,
            options: [],
            append() {{}},
            appendChild() {{}},
            addEventListener() {{}},
            querySelectorAll() {{ return []; }},
            getBoundingClientRect() {{
              return {{ top: 0, left: 0, width: 0, height: 0 }};
            }},
            showModal() {{}},
            close() {{}},
          }};
        }}

        const filtersForm = {{
          elements: [
            createField('limit', '24', '24'),
          ],
          reset() {{}},
          addEventListener() {{}},
        }};

        const displayPresetSelect = createSelect();
        const elementIds = [
          'filters-form',
          'pdb-id-query',
          'run-id',
          'limit',
          'sort-by',
          'sort-dir',
          'prev-page',
          'next-page',
          'reset-filters',
          'loading-state',
          'empty-state',
          'results-grid',
          'db-pill',
          'result-pill',
          'page-status',
          'detail-dialog',
          'close-detail',
          'detail-title',
          'detail-links',
          'detail-metrics',
          'detail-volumes',
          'viewer-host',
          'display-props-toggle',
          'display-props-popup',
          'filter-preset-select',
          'filter-preset-load',
          'filter-preset-save',
          'filter-preset-delete',
          'display-preset-select',
          'display-preset-load',
          'display-preset-save',
          'display-preset-delete',
        ];

        const elementsById = Object.fromEntries(
          elementIds.map((id) => [id, stubElement()]),
        );
        elementsById['filters-form'] = filtersForm;
        elementsById['pdb-id-query'] = createField('pdb_id_query');
        elementsById['limit'] = filtersForm.elements[0];
        elementsById['sort-by'] = createField('sort_by', 'largest_pore_volume');
        elementsById['sort-dir'] = createField('sort_dir', 'desc');
        elementsById['detail-dialog'] = stubElement();
        elementsById['display-props-popup'] = stubElement();
        elementsById['display-preset-select'] = displayPresetSelect;

        const document = {{
          getElementById(id) {{
            return elementsById[id] || stubElement();
          }},
          createElement() {{
            return {{
              value: '',
              textContent: '',
            }};
          }},
        }};

        const context = {{
          console,
          document,
          window: {{ addEventListener() {{}}, molstar: null }},
          localStorage: {{
            getItem(key) {{
              if (key === 'volumizer_display_presets') {{
                return JSON.stringify({{ Custom: ['num_hubs'] }});
              }}
              return null;
            }},
            setItem() {{}},
            removeItem() {{}},
          }},
          URLSearchParams,
          fetch: async () => ({{ ok: true, json: async () => ({{}}) }}),
          confirm: () => true,
          prompt: () => '',
          alert: () => {{}},
          Set,
          Number,
          String,
          JSON,
          Promise,
          Object,
          Array,
          Math,
        }};

        vm.createContext(context);
        vm.runInContext(source, context);
        context.__galleryTest.populateDisplayPresetSelect();

        const poresValue = (
          context.__galleryTest.BUILTIN_DISPLAY_PRESET_PREFIX + 'Pores'
        );
        const pocketsValue = (
          context.__galleryTest.BUILTIN_DISPLAY_PRESET_PREFIX + 'Pockets'
        );
        const cavitiesValue = (
          context.__galleryTest.BUILTIN_DISPLAY_PRESET_PREFIX + 'Cavities'
        );

        console.log(JSON.stringify({{
          option_labels: context.__galleryTest.elements.displayPresetSelect.options.map(
            (option) => option.textContent,
          ),
          default_preset: context.__galleryTest.getDisplayPresetData(''),
          pore_preset: context.__galleryTest.getDisplayPresetData(poresValue),
          pocket_preset: context.__galleryTest.getDisplayPresetData(pocketsValue),
          cavity_preset: context.__galleryTest.getDisplayPresetData(cavitiesValue),
        }}));
        """
    )

    result = subprocess.run(
        ["node", "-e", node_script],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(result.stdout)

    assert payload["option_labels"] == ["Default", "Pores", "Pockets", "Cavities", "Custom"]
    assert payload["default_preset"] == [
        "num_chains",
        "num_residues",
        "num_sequence_unique_chains",
    ]
    assert payload["pore_preset"] == [
        "num_chains",
        "num_residues",
        "num_sequence_unique_chains",
        "largest_pore_volume_a3",
        "largest_pore_min_diameter_a",
        "largest_pore_max_diameter_a",
        "largest_pore_uniformity",
        "largest_pore_circularity",
    ]
    assert payload["pocket_preset"] == [
        "num_chains",
        "num_residues",
        "num_sequence_unique_chains",
        "largest_pocket_volume_a3",
        "largest_pocket_max_diameter_a",
        "largest_pocket_length_a",
        "largest_pocket_circularity",
    ]
    assert payload["cavity_preset"] == [
        "num_chains",
        "num_residues",
        "num_sequence_unique_chains",
        "largest_cavity_volume_a3",
        "largest_cavity_max_diameter_a",
        "largest_cavity_length_a",
        "largest_cavity_circularity",
        "largest_cavity_uniformity",
    ]


def test_gallery_web_detail_and_viewer_data(tmp_path: Path):
    db_path, _, structure_ids = _build_web_fixture(tmp_path)
    structure_id = structure_ids["hit-a"]
    client = TestClient(create_app(db_path))

    detail_response = client.get(f"/api/hits/{structure_id}")
    assert detail_response.status_code == 200
    detail_payload = detail_response.json()
    assert detail_payload["source_label"] == "hit-a"
    assert len(detail_payload["volumes"]) == 1
    assert detail_payload["volumes"][0]["kind"] == "pore"
    assert detail_payload["volumes"][0]["cross_section_circularity"] == 0.85
    assert detail_payload["volumes"][0]["cross_section_uniformity"] == 0.72

    viewer_response = client.get(f"/api/hits/{structure_id}/viewer-data")
    assert viewer_response.status_code == 200
    viewer_payload = viewer_response.json()
    assert viewer_payload["structure_format"] == "mmcif"
    assert viewer_payload["structure_url"] == f"/files/structure/{structure_id}"
    assert viewer_payload["annotation_url"] == f"/files/annotation/{structure_id}"
    assert viewer_payload["thumbnails"]["x"] == f"/files/render/{structure_id}/x.png"


def test_gallery_web_serves_artifact_files(tmp_path: Path):
    db_path, _, structure_ids = _build_web_fixture(tmp_path)
    structure_id = structure_ids["hit-a"]
    client = TestClient(create_app(db_path))

    structure_response = client.get(f"/files/structure/{structure_id}")
    assert structure_response.status_code == 200
    assert "data_test" in structure_response.text

    annotation_response = client.get(f"/files/annotation/{structure_id}")
    assert annotation_response.status_code == 200
    assert annotation_response.json()["volumes"][0]["type"] == "pore"

    render_response = client.get(f"/files/render/{structure_id}/x.png")
    assert render_response.status_code == 200
    assert render_response.content == b"png-x"
