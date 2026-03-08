import json
from pathlib import Path
import sqlite3

import pytest

from volumizer import gallery_index
from volumizer.paths import TEST_DIR


TEST_INPUT_PDB = TEST_DIR / "pdbs" / "cavity.pdb"


def _write_summary(
    summary_path: Path,
    annotation_path: Path,
    structure_output_path: Path,
) -> None:
    payload = {
        "config": {
            "assembly_policy": "biological",
            "resolution": 3.0,
            "keep_non_protein": False,
            "output_dir": str(summary_path.parent),
        },
        "results": [
            {
                "source": "hit-a",
                "input_path": str(TEST_INPUT_PDB),
                "structure_output": str(structure_output_path),
                "annotation_output": str(annotation_path),
            }
        ],
        "errors": [],
        "skipped": [],
        "planned": [],
    }
    summary_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def test_build_gallery_index_from_annotation_payload_records(tmp_path: Path):
    summary_path = tmp_path / "run.summary.json"
    annotation_path = tmp_path / "hit-a.annotation.json"
    structure_output_path = tmp_path / "hit-a.annotated.cif"
    db_path = tmp_path / "gallery.db"

    structure_output_path.write_text("data_test\n#\n", encoding="utf-8")
    annotation_path.write_text(
        json.dumps(
            {
                "source": "hit-a",
                "num_volumes": 2,
                "volumes": [
                    {
                        "id": 0,
                        "type": "pore",
                        "volume": 200.0,
                        "x": 12.0,
                        "y": 5.0,
                        "z": 3.0,
                    },
                    {
                        "id": 0,
                        "type": "pocket",
                        "volume": 40.0,
                        "x": 6.0,
                        "y": 4.0,
                        "z": 2.0,
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    _write_summary(summary_path, annotation_path, structure_output_path)

    result = gallery_index.build_gallery_index(
        summary_path=summary_path,
        db_path=db_path,
        run_id="phase-a-test",
        replace_run=False,
        strict=True,
    )

    assert result["run_id"] == "phase-a-test"
    assert result["indexed_structures"] == 1
    assert result["indexed_volumes"] == 2
    assert result["skipped_structures"] == 0

    with sqlite3.connect(db_path) as connection:
        run_row = connection.execute(
            "SELECT run_id, assembly_policy, resolution FROM runs"
        ).fetchone()
        assert run_row == ("phase-a-test", "biological", 3.0)

        structure_row = connection.execute(
            """
            SELECT source_label, num_chains, num_residues, num_sequence_unique_chains
            FROM structures
            """
        ).fetchone()
        assert structure_row is not None
        assert structure_row[0] == "hit-a"
        assert structure_row[1] is not None and structure_row[1] > 0
        assert structure_row[2] is not None and structure_row[2] > 0
        assert structure_row[3] is not None and structure_row[3] > 0

        pore_row = connection.execute(
            """
            SELECT kind, rank_in_kind, volume_a3, length_a, min_diameter_a, max_diameter_a
            FROM volumes
            WHERE kind = 'pore'
            """
        ).fetchone()
        assert pore_row == ("pore", 1, 200.0, 12.0, 3.0, 5.0)

        aggregate_row = connection.execute(
            """
            SELECT num_pores, largest_pore_volume_a3, largest_pore_length_a,
                   largest_pore_min_diameter_a, largest_pore_max_diameter_a
            FROM structure_aggregates
            """
        ).fetchone()
        assert aggregate_row == (1, 200.0, 12.0, 3.0, 5.0)

        render_row = connection.execute(
            "SELECT render_status FROM renders"
        ).fetchone()
        assert render_row == ("pending",)


def test_build_gallery_index_supports_dataframe_oriented_json_and_replace(tmp_path: Path):
    summary_path = tmp_path / "run.summary.json"
    annotation_path = tmp_path / "hit-a.annotation.json"
    structure_output_path = tmp_path / "hit-a.annotated.cif"
    db_path = tmp_path / "gallery.db"

    structure_output_path.write_text("data_test\n#\n", encoding="utf-8")

    def _write_df_payload(volume: float) -> None:
        annotation_path.write_text(
            json.dumps(
                {
                    "id": {"0": 0},
                    "type": {"0": "pore"},
                    "volume": {"0": volume},
                    "x": {"0": 9.0},
                    "y": {"0": 4.0},
                    "z": {"0": 2.0},
                }
            ),
            encoding="utf-8",
        )

    _write_df_payload(120.0)
    _write_summary(summary_path, annotation_path, structure_output_path)

    gallery_index.build_gallery_index(
        summary_path=summary_path,
        db_path=db_path,
        run_id="replace-test",
        replace_run=False,
        strict=True,
    )

    with pytest.raises(ValueError, match="Run already exists"):
        gallery_index.build_gallery_index(
            summary_path=summary_path,
            db_path=db_path,
            run_id="replace-test",
            replace_run=False,
            strict=True,
        )

    _write_df_payload(222.0)
    gallery_index.build_gallery_index(
        summary_path=summary_path,
        db_path=db_path,
        run_id="replace-test",
        replace_run=True,
        strict=True,
    )

    with sqlite3.connect(db_path) as connection:
        rows = connection.execute(
            """
            SELECT volume_a3, length_a, min_diameter_a, max_diameter_a
            FROM volumes
            """
        ).fetchall()

    assert rows == [(222.0, 9.0, 2.0, 4.0)]
