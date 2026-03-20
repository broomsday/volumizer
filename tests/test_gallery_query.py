import json
from pathlib import Path
import sqlite3

import pytest

from volumizer import gallery_index, gallery_query
from volumizer.paths import TEST_DIR


TEST_INPUT_PDB = TEST_DIR / "pdbs" / "cavity.pdb"


def _write_annotation(
    path: Path,
    volumes: list[dict],
    frac_alpha: float | None = None,
    frac_beta: float | None = None,
    frac_coil: float | None = None,
) -> None:
    payload: dict = {
        "source": path.stem,
        "num_volumes": len(volumes),
        "volumes": volumes,
    }
    if frac_alpha is not None:
        payload["frac_alpha"] = frac_alpha
    if frac_beta is not None:
        payload["frac_beta"] = frac_beta
    if frac_coil is not None:
        payload["frac_coil"] = frac_coil
    path.write_text(json.dumps(payload), encoding="utf-8")


def _build_query_fixture_db(tmp_path: Path) -> tuple[Path, str]:
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)

    annotations = {
        "hit-a": {
            "volumes": [
                {"id": 0, "type": "pore", "volume": 200.0, "x": 12.0, "y": 5.0, "z": 3.0,
                 "cross_section_circularity": 0.85, "cross_section_uniformity": 0.72},
                {"id": 0, "type": "pocket", "volume": 40.0, "x": 6.0, "y": 4.0, "z": 2.0,
                 "cross_section_circularity": 0.60, "cross_section_uniformity": 0.55},
            ],
            "frac_alpha": 0.50,
            "frac_beta": 0.20,
            "frac_coil": 0.30,
        },
        "hit-b": {
            "volumes": [
                {"id": 0, "type": "pore", "volume": 80.0, "x": 7.0, "y": 3.0, "z": 2.0,
                 "cross_section_circularity": 0.90, "cross_section_uniformity": 0.80},
            ],
            "frac_alpha": 0.10,
            "frac_beta": 0.60,
            "frac_coil": 0.30,
        },
        "hit-c": {
            "volumes": [
                {"id": 0, "type": "pocket", "volume": 30.0, "x": 5.0, "y": 2.0, "z": 1.0,
                 "cross_section_circularity": 0.45, "cross_section_uniformity": 0.35},
            ],
            "frac_alpha": 0.80,
            "frac_beta": 0.05,
            "frac_coil": 0.15,
        },
    }

    results = []
    for source_label, entry in annotations.items():
        annotation_path = run_dir / f"{source_label}.annotation.json"
        structure_output_path = run_dir / f"{source_label}.annotated.cif"
        _write_annotation(
            annotation_path,
            entry["volumes"],
            frac_alpha=entry.get("frac_alpha"),
            frac_beta=entry.get("frac_beta"),
            frac_coil=entry.get("frac_coil"),
        )
        structure_output_path.write_text("data_test\n#\n", encoding="utf-8")

        results.append(
            {
                "source": source_label,
                "input_path": str(TEST_INPUT_PDB),
                "structure_output": str(structure_output_path),
                "annotation_output": str(annotation_path),
            }
        )

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
    run_id = "query-fixture"
    gallery_index.build_gallery_index(
        summary_path=summary_path,
        db_path=db_path,
        run_id=run_id,
        replace_run=False,
        strict=True,
    )

    with sqlite3.connect(db_path) as connection:
        connection.execute(
            "UPDATE structures SET num_chains = 1, num_residues = 100, num_sequence_unique_chains = 1 WHERE source_label = 'hit-a'"
        )
        connection.execute(
            "UPDATE structures SET num_chains = 3, num_residues = 250, num_sequence_unique_chains = 2 WHERE source_label = 'hit-b'"
        )
        connection.execute(
            "UPDATE structures SET num_chains = 2, num_residues = 180, num_sequence_unique_chains = 2 WHERE source_label = 'hit-c'"
        )
        connection.commit()

    return db_path, run_id


def test_query_gallery_index_default_sort_and_pagination(tmp_path: Path):
    db_path, run_id = _build_query_fixture_db(tmp_path)

    page_1 = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        limit=2,
        offset=0,
    )

    assert page_1["total_count"] == 3
    assert len(page_1["rows"]) == 2
    assert [row["source_label"] for row in page_1["rows"]] == ["hit-a", "hit-b"]

    page_2 = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        limit=1,
        offset=1,
    )
    assert page_2["total_count"] == 3
    assert [row["source_label"] for row in page_2["rows"]] == ["hit-b"]


def test_query_gallery_index_filters_pore_and_structure_ranges(tmp_path: Path):
    db_path, run_id = _build_query_fixture_db(tmp_path)

    pore_filtered = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        pore_volume_min=100.0,
    )
    assert [row["source_label"] for row in pore_filtered["rows"]] == ["hit-a"]

    diameter_filtered = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        pore_dmax_max=3.5,
    )
    assert [row["source_label"] for row in diameter_filtered["rows"]] == ["hit-b"]

    structure_filtered = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        chains_min=2,
        chains_max=2,
        residues_min=150,
        residues_max=200,
        seq_unique_chains_min=2,
    )
    assert [row["source_label"] for row in structure_filtered["rows"]] == ["hit-c"]


def test_query_gallery_index_filters_pocket_and_cavity_ranges(tmp_path: Path):
    db_path, run_id = _build_query_fixture_db(tmp_path)

    pocket_filtered = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        pocket_volume_min=25.0,
    )
    labels = [row["source_label"] for row in pocket_filtered["rows"]]
    assert "hit-a" in labels
    assert "hit-c" in labels
    assert "hit-b" not in labels

    pocket_only = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        num_pockets_min=1,
        num_pores_max=0,
    )
    assert [row["source_label"] for row in pocket_only["rows"]] == ["hit-c"]


def test_query_gallery_index_returns_and_filters_sse_fractions(tmp_path: Path):
    db_path, run_id = _build_query_fixture_db(tmp_path)

    result = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        sort_by="frac_alpha",
        sort_dir="desc",
    )
    labels = [row["source_label"] for row in result["rows"]]
    assert labels[0] == "hit-c"
    assert result["rows"][0]["frac_alpha"] == 0.80

    high_alpha = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        frac_alpha_min=0.40,
    )
    labels = [row["source_label"] for row in high_alpha["rows"]]
    assert "hit-a" in labels
    assert "hit-c" in labels
    assert "hit-b" not in labels


def test_query_gallery_index_sorts_by_pocket_volume(tmp_path: Path):
    db_path, run_id = _build_query_fixture_db(tmp_path)

    result = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        sort_by="largest_pocket_volume",
        sort_dir="desc",
    )
    labels = [row["source_label"] for row in result["rows"]]
    assert labels[0] == "hit-a"
    assert result["rows"][0]["largest_pocket_volume_a3"] == 40.0


def test_query_gallery_index_filters_and_sorts_by_cross_section_metrics(tmp_path: Path):
    db_path, run_id = _build_query_fixture_db(tmp_path)

    # hit-b has pore circularity 0.90, hit-a has 0.85
    sorted_by_circ = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        sort_by="largest_pore_circularity",
        sort_dir="desc",
    )
    labels = [row["source_label"] for row in sorted_by_circ["rows"]]
    assert labels[0] == "hit-b"
    assert sorted_by_circ["rows"][0]["largest_pore_circularity"] == 0.90

    # Filter pore uniformity >= 0.75 should only match hit-b (0.80)
    filtered = gallery_query.query_gallery_index(
        db_path=db_path,
        run_id=run_id,
        pore_uniformity_min=0.75,
    )
    assert [row["source_label"] for row in filtered["rows"]] == ["hit-b"]


def test_query_gallery_index_validates_sort_inputs(tmp_path: Path):
    db_path, run_id = _build_query_fixture_db(tmp_path)

    with pytest.raises(ValueError, match="Unsupported sort_by"):
        gallery_query.query_gallery_index(
            db_path=db_path,
            run_id=run_id,
            sort_by="unknown-column",
        )

    with pytest.raises(ValueError, match="sort_dir"):
        gallery_query.query_gallery_index(
            db_path=db_path,
            run_id=run_id,
            sort_dir="sideways",
        )
