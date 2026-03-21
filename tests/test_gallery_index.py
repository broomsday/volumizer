import json
from pathlib import Path
import sqlite3

import pytest

from volumizer import gallery_index, pdb
from volumizer.gallery_index import _sequence_overlap, _count_sequence_unique_chains
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
                "frac_alpha": 0.45,
                "frac_beta": 0.20,
                "frac_coil": 0.35,
                "volumes": [
                    {
                        "id": 0,
                        "type": "pore",
                        "volume": 200.0,
                        "x": 12.0,
                        "y": 5.0,
                        "z": 3.0,
                        "cross_section_circularity": 0.85,
                        "cross_section_uniformity": 0.72,
                    },
                    {
                        "id": 0,
                        "type": "pocket",
                        "volume": 40.0,
                        "x": 6.0,
                        "y": 4.0,
                        "z": 2.0,
                        "cross_section_circularity": 0.60,
                        "cross_section_uniformity": 0.55,
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
            SELECT source_label, num_chains, num_residues, num_sequence_unique_chains,
                   frac_alpha, frac_beta, frac_coil
            FROM structures
            """
        ).fetchone()
        assert structure_row is not None
        assert structure_row[0] == "hit-a"
        assert structure_row[1] is not None and structure_row[1] > 0
        assert structure_row[2] is not None and structure_row[2] > 0
        assert structure_row[3] is not None and structure_row[3] > 0
        assert structure_row[4] == 0.45
        assert structure_row[5] == 0.20
        assert structure_row[6] == 0.35

        pore_row = connection.execute(
            """
            SELECT kind, rank_in_kind, volume_a3, length_a, min_diameter_a, max_diameter_a,
                   cross_section_circularity, cross_section_uniformity
            FROM volumes
            WHERE kind = 'pore'
            """
        ).fetchone()
        assert pore_row == ("pore", 1, 200.0, 12.0, 3.0, 5.0, 0.85, 0.72)

        aggregate_row = connection.execute(
            """
            SELECT num_pores, largest_pore_volume_a3, largest_pore_length_a,
                   largest_pore_min_diameter_a, largest_pore_max_diameter_a,
                   largest_pore_circularity, largest_pore_uniformity,
                   num_pockets, largest_pocket_volume_a3, largest_pocket_length_a,
                   largest_pocket_min_diameter_a, largest_pocket_max_diameter_a,
                   largest_pocket_circularity, largest_pocket_uniformity,
                   num_cavities, num_hubs
            FROM structure_aggregates
            """
        ).fetchone()
        assert aggregate_row == (
            1, 200.0, 12.0, 3.0, 5.0, 0.85, 0.72,
            1, 40.0, 6.0, 2.0, 4.0, 0.60, 0.55,
            0, 0,
        )

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


def test_build_gallery_index_infers_pdb_id_from_source_when_missing(tmp_path: Path):
    summary_path = tmp_path / "run.summary.json"
    annotation_path = tmp_path / "1dzf.annotation.json"
    structure_output_path = tmp_path / "1dzf.annotated.cif"
    db_path = tmp_path / "gallery.db"

    structure_output_path.write_text("data_test\n#\n", encoding="utf-8")
    annotation_path.write_text(
        json.dumps(
            {
                "source": "1dzf",
                "num_volumes": 1,
                "volumes": [
                    {
                        "id": 0,
                        "type": "pore",
                        "volume": 111.0,
                        "x": 9.0,
                        "y": 4.0,
                        "z": 2.0,
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    summary_path.write_text(
        json.dumps(
            {
                "config": {
                    "assembly_policy": "biological",
                    "resolution": 3.0,
                    "keep_non_protein": False,
                    "output_dir": str(tmp_path),
                },
                "results": [
                    {
                        "source": "1dzf",
                        "input_path": str(TEST_INPUT_PDB),
                        "structure_output": str(structure_output_path),
                        "annotation_output": str(annotation_path),
                    }
                ],
                "errors": [],
                "skipped": [],
                "planned": [],
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    gallery_index.build_gallery_index(
        summary_path=summary_path,
        db_path=db_path,
        run_id="pdb-id-infer-test",
        replace_run=False,
        strict=True,
    )

    with sqlite3.connect(db_path) as connection:
        row = connection.execute(
            "SELECT source_label, pdb_id FROM structures"
        ).fetchone()

    assert row == ("1dzf", "1DZF")


class TestSequenceOverlap:
    def test_identical_sequences(self):
        seq = ("ALA", "GLY", "VAL")
        assert _sequence_overlap(seq, seq) == 1.0

    def test_empty_sequences(self):
        assert _sequence_overlap((), ()) == 1.0

    def test_one_empty(self):
        assert _sequence_overlap(("ALA",), ()) == 0.0

    def test_prefix_truncation(self):
        full = ("GLY", "GLY", "MET", "GLU", "LYS", "ALA", "VAL", "THR", "PRO", "LEU")
        truncated = ("MET", "GLU", "LYS", "ALA", "VAL", "THR", "PRO", "LEU")
        overlap = _sequence_overlap(full, truncated)
        # 8 matches out of 10 max length = 0.8
        assert overlap == pytest.approx(0.8)

    def test_suffix_truncation(self):
        full = ("ALA", "GLY", "VAL", "LEU", "ILE")
        truncated = ("ALA", "GLY", "VAL")
        overlap = _sequence_overlap(full, truncated)
        assert overlap == pytest.approx(0.6)

    def test_completely_different(self):
        seq1 = ("ALA", "ALA", "ALA")
        seq2 = ("GLY", "GLY", "GLY")
        assert _sequence_overlap(seq1, seq2) == 0.0


class TestCountSequenceUniqueChains:
    def test_all_identical(self):
        seqs = {
            "A": ("ALA", "GLY", "VAL"),
            "B": ("ALA", "GLY", "VAL"),
            "C": ("ALA", "GLY", "VAL"),
        }
        assert _count_sequence_unique_chains(seqs) == 1

    def test_all_different(self):
        seqs = {
            "A": ("ALA",) * 100,
            "B": ("GLY",) * 100,
        }
        assert _count_sequence_unique_chains(seqs) == 2

    def test_fuzzy_match_terminal_truncation(self):
        # 100-residue chain, one copy missing 3 N-terminal residues
        base = tuple(["ALA"] * 50 + ["GLY"] * 50)
        truncated = base[3:]
        seqs = {"A": base, "B": truncated}
        # 97/100 = 0.97 > 0.95 threshold
        assert _count_sequence_unique_chains(seqs) == 1

    def test_below_threshold_counted_as_different(self):
        base = tuple(["ALA"] * 100)
        # Differ by 6 residues at the start -> 94/100 = 0.94 < 0.95
        different = tuple(["GLY"] * 6 + ["ALA"] * 94)
        seqs = {"A": base, "B": different}
        assert _count_sequence_unique_chains(seqs) == 2

    def test_empty(self):
        assert _count_sequence_unique_chains({}) == 0

    def test_single_chain(self):
        assert _count_sequence_unique_chains({"A": ("ALA",)}) == 1


RCSB_2ZBT = Path("data/runs/rcsb70/downloads/2ZBT.cif")


@pytest.mark.skipif(not RCSB_2ZBT.exists(), reason="2ZBT test data not available")
class TestBiologicalAssemblyMetrics2ZBT:
    def test_2zbt_has_12_chains_1_unique(self):
        chains, residues, unique = gallery_index._compute_structure_metrics(
            RCSB_2ZBT, "biological"
        )
        assert chains == 12
        assert unique == 1
        assert residues is not None and residues > 3000
