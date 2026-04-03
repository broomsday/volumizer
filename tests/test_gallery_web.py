import json
from pathlib import Path
import sqlite3

from fastapi.testclient import TestClient

from volumizer import gallery_index, gallery_render
from volumizer.paths import TEST_DIR
from volumizer.web.app import create_app


TEST_INPUT_PDB = TEST_DIR / "pdbs" / "cavity.pdb"


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
    client = TestClient(create_app(db_path))

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


def test_gallery_web_static_search_serializes_form_fields_generically():
    static_dir = Path(__file__).resolve().parents[1] / "volumizer" / "web" / "static"
    app_js = (static_dir / "app.js").read_text(encoding="utf-8")
    index_html = (static_dir / "index.html").read_text(encoding="utf-8")

    assert "for (const field of elements.filtersForm.elements)" in app_js
    assert "const numericNames = [" not in app_js
    assert 'name="seq_unique_chains_max"' in index_html
    assert 'name="chains_max"' in index_html
    assert 'name="residues_max"' in index_html


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
