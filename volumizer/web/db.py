"""
Database helpers for the local gallery web app.
"""

from __future__ import annotations

from pathlib import Path
import sqlite3
from typing import Any


_STRUCTURE_FORMAT_BY_SUFFIX = {
    ".bcif": "bcif",
    ".cif": "mmcif",
    ".mmcif": "mmcif",
    ".pdb": "pdb",
}

VALID_RENDER_AXES = frozenset({"x", "y", "z"})


def _connect(db_path: Path) -> sqlite3.Connection:
    resolved = Path(db_path).resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"Gallery DB does not exist: {resolved}")

    connection = sqlite3.connect(resolved)
    connection.row_factory = sqlite3.Row
    return connection


def resolve_artifact_path(raw_path: str | None) -> Path | None:
    if raw_path is None:
        return None

    path = Path(str(raw_path)).expanduser()
    if not path.is_absolute():
        path = path.resolve()
    return path


def detect_structure_format(structure_path: Path | None) -> str | None:
    if structure_path is None:
        return None
    return _STRUCTURE_FORMAT_BY_SUFFIX.get(structure_path.suffix.lower())


def list_runs(db_path: Path) -> list[dict[str, Any]]:
    with _connect(db_path) as connection:
        rows = connection.execute(
            """
            SELECT
                r.run_id,
                r.created_at,
                r.source_summary_path,
                r.resolution,
                r.assembly_policy,
                COUNT(s.structure_id) AS structure_count
            FROM runs r
            LEFT JOIN structures s ON s.run_id = r.run_id
            GROUP BY
                r.run_id,
                r.created_at,
                r.source_summary_path,
                r.resolution,
                r.assembly_policy
            ORDER BY r.created_at DESC, r.run_id ASC
            """
        ).fetchall()

    return [dict(row) for row in rows]


def get_hit_detail(db_path: Path, structure_id: int) -> dict[str, Any] | None:
    with _connect(db_path) as connection:
        row = connection.execute(
            """
            SELECT
                s.structure_id,
                s.run_id,
                s.source_label,
                s.pdb_id,
                s.input_path,
                s.annotated_cif_path,
                s.annotation_json_path,
                s.num_chains,
                s.num_residues,
                s.num_sequence_unique_chains,
                s.frac_alpha,
                s.frac_beta,
                s.frac_coil,
                a.num_pores,
                a.largest_pore_volume_a3,
                a.largest_pore_length_a,
                a.largest_pore_min_diameter_a,
                a.largest_pore_max_diameter_a,
                a.num_pockets,
                a.largest_pocket_volume_a3,
                a.largest_pocket_length_a,
                a.largest_pocket_min_diameter_a,
                a.largest_pocket_max_diameter_a,
                a.num_cavities,
                a.largest_cavity_volume_a3,
                a.largest_cavity_length_a,
                a.largest_cavity_min_diameter_a,
                a.largest_cavity_max_diameter_a,
                a.num_hubs,
                a.largest_hub_volume_a3,
                a.largest_hub_length_a,
                a.largest_hub_min_diameter_a,
                a.largest_hub_max_diameter_a,
                a.largest_pore_circularity,
                a.largest_pore_uniformity,
                a.largest_pocket_circularity,
                a.largest_pocket_uniformity,
                a.largest_cavity_circularity,
                a.largest_cavity_uniformity,
                a.largest_hub_circularity,
                a.largest_hub_uniformity,
                r.x_png_path,
                r.y_png_path,
                r.z_png_path,
                r.render_style_hash,
                r.render_status,
                r.render_error,
                r.updated_at
            FROM structures s
            LEFT JOIN structure_aggregates a ON a.structure_id = s.structure_id
            LEFT JOIN renders r ON r.structure_id = s.structure_id
            WHERE s.structure_id = ?
            """,
            (int(structure_id),),
        ).fetchone()
        if row is None:
            return None

        volume_rows = connection.execute(
            """
            SELECT
                volume_id,
                structure_id,
                kind,
                rank_in_kind,
                volume_a3,
                length_a,
                min_diameter_a,
                max_diameter_a,
                centroid_x,
                centroid_y,
                centroid_z,
                cross_section_circularity,
                cross_section_uniformity
            FROM volumes
            WHERE structure_id = ?
            ORDER BY kind ASC, rank_in_kind ASC
            """,
            (int(structure_id),),
        ).fetchall()

    detail = dict(row)
    detail["volumes"] = [dict(volume_row) for volume_row in volume_rows]
    return detail


def get_hit_artifacts(db_path: Path, structure_id: int) -> dict[str, Any] | None:
    with _connect(db_path) as connection:
        row = connection.execute(
            """
            SELECT
                s.structure_id,
                s.source_label,
                s.annotated_cif_path,
                s.annotation_json_path,
                r.x_png_path,
                r.y_png_path,
                r.z_png_path
            FROM structures s
            LEFT JOIN renders r ON r.structure_id = s.structure_id
            WHERE s.structure_id = ?
            """,
            (int(structure_id),),
        ).fetchone()

    return None if row is None else dict(row)
