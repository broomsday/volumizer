"""
Query helpers for local gallery SQLite index.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sqlite3
from typing import Any


_SORT_COLUMN_MAP = {
    "structure_id": "s.structure_id",
    "source_label": "s.source_label",
    "num_chains": "s.num_chains",
    "num_residues": "s.num_residues",
    "num_sequence_unique_chains": "s.num_sequence_unique_chains",
    "frac_alpha": "s.frac_alpha",
    "frac_beta": "s.frac_beta",
    "frac_coil": "s.frac_coil",
    "num_pores": "a.num_pores",
    "largest_pore_volume": "a.largest_pore_volume_a3",
    "largest_pore_length": "a.largest_pore_length_a",
    "largest_pore_min_diameter": "a.largest_pore_min_diameter_a",
    "largest_pore_max_diameter": "a.largest_pore_max_diameter_a",
    "num_pockets": "a.num_pockets",
    "largest_pocket_volume": "a.largest_pocket_volume_a3",
    "largest_pocket_length": "a.largest_pocket_length_a",
    "largest_pocket_min_diameter": "a.largest_pocket_min_diameter_a",
    "largest_pocket_max_diameter": "a.largest_pocket_max_diameter_a",
    "num_cavities": "a.num_cavities",
    "largest_cavity_volume": "a.largest_cavity_volume_a3",
    "largest_cavity_length": "a.largest_cavity_length_a",
    "largest_cavity_min_diameter": "a.largest_cavity_min_diameter_a",
    "largest_cavity_max_diameter": "a.largest_cavity_max_diameter_a",
    "num_hubs": "a.num_hubs",
    "largest_hub_volume": "a.largest_hub_volume_a3",
    "largest_hub_length": "a.largest_hub_length_a",
    "largest_hub_min_diameter": "a.largest_hub_min_diameter_a",
    "largest_hub_max_diameter": "a.largest_hub_max_diameter_a",
    "largest_pore_circularity": "a.largest_pore_circularity",
    "largest_pore_uniformity": "a.largest_pore_uniformity",
    "largest_pocket_circularity": "a.largest_pocket_circularity",
    "largest_pocket_uniformity": "a.largest_pocket_uniformity",
    "largest_cavity_circularity": "a.largest_cavity_circularity",
    "largest_cavity_uniformity": "a.largest_cavity_uniformity",
    "largest_hub_circularity": "a.largest_hub_circularity",
    "largest_hub_uniformity": "a.largest_hub_uniformity",
}
_DEFAULT_SORT_BY = "largest_pore_volume"
_DEFAULT_SORT_DIR = "desc"


def _safe_non_negative_int(value: int, name: str) -> int:
    normalized = int(value)
    if normalized < 0:
        raise ValueError(f"{name} must be >= 0")
    return normalized


def _normalize_sort(sort_by: str, sort_dir: str) -> tuple[str, str]:
    normalized_sort_by = str(sort_by).strip()
    if normalized_sort_by not in _SORT_COLUMN_MAP:
        raise ValueError(
            f"Unsupported sort_by: {sort_by}. "
            f"Expected one of: {', '.join(sorted(_SORT_COLUMN_MAP.keys()))}"
        )

    normalized_sort_dir = str(sort_dir).strip().lower()
    if normalized_sort_dir not in {"asc", "desc"}:
        raise ValueError("sort_dir must be one of: asc, desc")

    return normalized_sort_by, normalized_sort_dir


def _append_range_filter(
    where_clauses: list[str],
    params: list[Any],
    column: str,
    min_value: float | int | None,
    max_value: float | int | None,
) -> None:
    if min_value is not None:
        where_clauses.append(f"{column} >= ?")
        params.append(min_value)
    if max_value is not None:
        where_clauses.append(f"{column} <= ?")
        params.append(max_value)


def query_gallery_index(
    db_path: Path,
    run_id: str | None = None,
    pore_volume_min: float | None = None,
    pore_volume_max: float | None = None,
    pore_length_min: float | None = None,
    pore_length_max: float | None = None,
    pore_dmin_min: float | None = None,
    pore_dmin_max: float | None = None,
    pore_dmax_min: float | None = None,
    pore_dmax_max: float | None = None,
    num_pores_min: int | None = None,
    num_pores_max: int | None = None,
    pocket_volume_min: float | None = None,
    pocket_volume_max: float | None = None,
    pocket_length_min: float | None = None,
    pocket_length_max: float | None = None,
    pocket_dmin_min: float | None = None,
    pocket_dmin_max: float | None = None,
    pocket_dmax_min: float | None = None,
    pocket_dmax_max: float | None = None,
    num_pockets_min: int | None = None,
    num_pockets_max: int | None = None,
    cavity_volume_min: float | None = None,
    cavity_volume_max: float | None = None,
    cavity_length_min: float | None = None,
    cavity_length_max: float | None = None,
    cavity_dmin_min: float | None = None,
    cavity_dmin_max: float | None = None,
    cavity_dmax_min: float | None = None,
    cavity_dmax_max: float | None = None,
    num_cavities_min: int | None = None,
    num_cavities_max: int | None = None,
    hub_volume_min: float | None = None,
    hub_volume_max: float | None = None,
    hub_length_min: float | None = None,
    hub_length_max: float | None = None,
    hub_dmin_min: float | None = None,
    hub_dmin_max: float | None = None,
    hub_dmax_min: float | None = None,
    hub_dmax_max: float | None = None,
    num_hubs_min: int | None = None,
    num_hubs_max: int | None = None,
    chains_min: int | None = None,
    chains_max: int | None = None,
    residues_min: int | None = None,
    residues_max: int | None = None,
    seq_unique_chains_min: int | None = None,
    seq_unique_chains_max: int | None = None,
    frac_alpha_min: float | None = None,
    frac_alpha_max: float | None = None,
    frac_beta_min: float | None = None,
    frac_beta_max: float | None = None,
    frac_coil_min: float | None = None,
    frac_coil_max: float | None = None,
    pore_circularity_min: float | None = None,
    pore_circularity_max: float | None = None,
    pore_uniformity_min: float | None = None,
    pore_uniformity_max: float | None = None,
    pocket_circularity_min: float | None = None,
    pocket_circularity_max: float | None = None,
    pocket_uniformity_min: float | None = None,
    pocket_uniformity_max: float | None = None,
    cavity_circularity_min: float | None = None,
    cavity_circularity_max: float | None = None,
    cavity_uniformity_min: float | None = None,
    cavity_uniformity_max: float | None = None,
    hub_circularity_min: float | None = None,
    hub_circularity_max: float | None = None,
    hub_uniformity_min: float | None = None,
    hub_uniformity_max: float | None = None,
    limit: int = 50,
    offset: int = 0,
    sort_by: str = _DEFAULT_SORT_BY,
    sort_dir: str = _DEFAULT_SORT_DIR,
) -> dict[str, Any]:
    db_path = Path(db_path).resolve()
    if not db_path.is_file():
        raise FileNotFoundError(f"Gallery DB does not exist: {db_path}")

    normalized_limit = _safe_non_negative_int(limit, "limit")
    normalized_offset = _safe_non_negative_int(offset, "offset")
    normalized_sort_by, normalized_sort_dir = _normalize_sort(sort_by, sort_dir)

    where_clauses = ["1=1"]
    params: list[Any] = []

    if run_id is not None:
        where_clauses.append("s.run_id = ?")
        params.append(str(run_id))

    _append_range_filter(
        where_clauses,
        params,
        "a.largest_pore_volume_a3",
        pore_volume_min,
        pore_volume_max,
    )
    _append_range_filter(
        where_clauses,
        params,
        "a.largest_pore_length_a",
        pore_length_min,
        pore_length_max,
    )
    _append_range_filter(
        where_clauses,
        params,
        "a.largest_pore_min_diameter_a",
        pore_dmin_min,
        pore_dmin_max,
    )
    _append_range_filter(
        where_clauses,
        params,
        "a.largest_pore_max_diameter_a",
        pore_dmax_min,
        pore_dmax_max,
    )
    _append_range_filter(where_clauses, params, "a.num_pores", num_pores_min, num_pores_max)

    _append_range_filter(
        where_clauses, params, "a.largest_pocket_volume_a3", pocket_volume_min, pocket_volume_max,
    )
    _append_range_filter(
        where_clauses, params, "a.largest_pocket_length_a", pocket_length_min, pocket_length_max,
    )
    _append_range_filter(
        where_clauses, params, "a.largest_pocket_min_diameter_a", pocket_dmin_min, pocket_dmin_max,
    )
    _append_range_filter(
        where_clauses, params, "a.largest_pocket_max_diameter_a", pocket_dmax_min, pocket_dmax_max,
    )
    _append_range_filter(where_clauses, params, "a.num_pockets", num_pockets_min, num_pockets_max)

    _append_range_filter(
        where_clauses, params, "a.largest_cavity_volume_a3", cavity_volume_min, cavity_volume_max,
    )
    _append_range_filter(
        where_clauses, params, "a.largest_cavity_length_a", cavity_length_min, cavity_length_max,
    )
    _append_range_filter(
        where_clauses, params, "a.largest_cavity_min_diameter_a", cavity_dmin_min, cavity_dmin_max,
    )
    _append_range_filter(
        where_clauses, params, "a.largest_cavity_max_diameter_a", cavity_dmax_min, cavity_dmax_max,
    )
    _append_range_filter(where_clauses, params, "a.num_cavities", num_cavities_min, num_cavities_max)

    _append_range_filter(
        where_clauses, params, "a.largest_hub_volume_a3", hub_volume_min, hub_volume_max,
    )
    _append_range_filter(
        where_clauses, params, "a.largest_hub_length_a", hub_length_min, hub_length_max,
    )
    _append_range_filter(
        where_clauses, params, "a.largest_hub_min_diameter_a", hub_dmin_min, hub_dmin_max,
    )
    _append_range_filter(
        where_clauses, params, "a.largest_hub_max_diameter_a", hub_dmax_min, hub_dmax_max,
    )
    _append_range_filter(where_clauses, params, "a.num_hubs", num_hubs_min, num_hubs_max)

    _append_range_filter(where_clauses, params, "s.num_chains", chains_min, chains_max)
    _append_range_filter(where_clauses, params, "s.num_residues", residues_min, residues_max)
    _append_range_filter(
        where_clauses,
        params,
        "s.num_sequence_unique_chains",
        seq_unique_chains_min,
        seq_unique_chains_max,
    )
    _append_range_filter(where_clauses, params, "s.frac_alpha", frac_alpha_min, frac_alpha_max)
    _append_range_filter(where_clauses, params, "s.frac_beta", frac_beta_min, frac_beta_max)
    _append_range_filter(where_clauses, params, "s.frac_coil", frac_coil_min, frac_coil_max)

    _append_range_filter(where_clauses, params, "a.largest_pore_circularity", pore_circularity_min, pore_circularity_max)
    _append_range_filter(where_clauses, params, "a.largest_pore_uniformity", pore_uniformity_min, pore_uniformity_max)
    _append_range_filter(where_clauses, params, "a.largest_pocket_circularity", pocket_circularity_min, pocket_circularity_max)
    _append_range_filter(where_clauses, params, "a.largest_pocket_uniformity", pocket_uniformity_min, pocket_uniformity_max)
    _append_range_filter(where_clauses, params, "a.largest_cavity_circularity", cavity_circularity_min, cavity_circularity_max)
    _append_range_filter(where_clauses, params, "a.largest_cavity_uniformity", cavity_uniformity_min, cavity_uniformity_max)
    _append_range_filter(where_clauses, params, "a.largest_hub_circularity", hub_circularity_min, hub_circularity_max)
    _append_range_filter(where_clauses, params, "a.largest_hub_uniformity", hub_uniformity_min, hub_uniformity_max)

    where_sql = " AND ".join(where_clauses)
    sort_sql = _SORT_COLUMN_MAP[normalized_sort_by]
    order_sql = (
        f"ORDER BY ({sort_sql} IS NULL) ASC, {sort_sql} {normalized_sort_dir.upper()}, "
        "s.structure_id ASC"
    )

    base_from_sql = (
        "FROM structures s "
        "INNER JOIN structure_aggregates a ON a.structure_id = s.structure_id "
        "LEFT JOIN renders r ON r.structure_id = s.structure_id "
        f"WHERE {where_sql}"
    )

    count_sql = f"SELECT COUNT(*) {base_from_sql}"
    query_sql = (
        "SELECT "
        "s.structure_id, s.run_id, s.source_label, s.pdb_id, s.input_path, "
        "s.annotated_cif_path, s.annotation_json_path, "
        "s.num_chains, s.num_residues, s.num_sequence_unique_chains, "
        "s.frac_alpha, s.frac_beta, s.frac_coil, "
        "a.num_pores, a.largest_pore_volume_a3, a.largest_pore_length_a, "
        "a.largest_pore_min_diameter_a, a.largest_pore_max_diameter_a, "
        "a.num_pockets, a.largest_pocket_volume_a3, a.largest_pocket_length_a, "
        "a.largest_pocket_min_diameter_a, a.largest_pocket_max_diameter_a, "
        "a.num_cavities, a.largest_cavity_volume_a3, a.largest_cavity_length_a, "
        "a.largest_cavity_min_diameter_a, a.largest_cavity_max_diameter_a, "
        "a.num_hubs, a.largest_hub_volume_a3, a.largest_hub_length_a, "
        "a.largest_hub_min_diameter_a, a.largest_hub_max_diameter_a, "
        "a.largest_pore_circularity, a.largest_pore_uniformity, "
        "a.largest_pocket_circularity, a.largest_pocket_uniformity, "
        "a.largest_cavity_circularity, a.largest_cavity_uniformity, "
        "a.largest_hub_circularity, a.largest_hub_uniformity, "
        "r.x_png_path, r.y_png_path, r.z_png_path, r.render_style_hash, r.render_status "
        f"{base_from_sql} "
        f"{order_sql} "
        "LIMIT ? OFFSET ?"
    )

    with sqlite3.connect(db_path) as connection:
        connection.row_factory = sqlite3.Row

        total_count = int(connection.execute(count_sql, params).fetchone()[0])
        row_params = [*params, normalized_limit, normalized_offset]
        rows = connection.execute(query_sql, row_params).fetchall()

    row_dicts = [dict(row) for row in rows]
    return {
        "db": str(db_path),
        "run_id": run_id,
        "filters": {
            "pore_volume_min": pore_volume_min,
            "pore_volume_max": pore_volume_max,
            "pore_length_min": pore_length_min,
            "pore_length_max": pore_length_max,
            "pore_dmin_min": pore_dmin_min,
            "pore_dmin_max": pore_dmin_max,
            "pore_dmax_min": pore_dmax_min,
            "pore_dmax_max": pore_dmax_max,
            "num_pores_min": num_pores_min,
            "num_pores_max": num_pores_max,
            "pocket_volume_min": pocket_volume_min,
            "pocket_volume_max": pocket_volume_max,
            "pocket_length_min": pocket_length_min,
            "pocket_length_max": pocket_length_max,
            "pocket_dmin_min": pocket_dmin_min,
            "pocket_dmin_max": pocket_dmin_max,
            "pocket_dmax_min": pocket_dmax_min,
            "pocket_dmax_max": pocket_dmax_max,
            "num_pockets_min": num_pockets_min,
            "num_pockets_max": num_pockets_max,
            "cavity_volume_min": cavity_volume_min,
            "cavity_volume_max": cavity_volume_max,
            "cavity_length_min": cavity_length_min,
            "cavity_length_max": cavity_length_max,
            "cavity_dmin_min": cavity_dmin_min,
            "cavity_dmin_max": cavity_dmin_max,
            "cavity_dmax_min": cavity_dmax_min,
            "cavity_dmax_max": cavity_dmax_max,
            "num_cavities_min": num_cavities_min,
            "num_cavities_max": num_cavities_max,
            "hub_volume_min": hub_volume_min,
            "hub_volume_max": hub_volume_max,
            "hub_length_min": hub_length_min,
            "hub_length_max": hub_length_max,
            "hub_dmin_min": hub_dmin_min,
            "hub_dmin_max": hub_dmin_max,
            "hub_dmax_min": hub_dmax_min,
            "hub_dmax_max": hub_dmax_max,
            "num_hubs_min": num_hubs_min,
            "num_hubs_max": num_hubs_max,
            "chains_min": chains_min,
            "chains_max": chains_max,
            "residues_min": residues_min,
            "residues_max": residues_max,
            "seq_unique_chains_min": seq_unique_chains_min,
            "seq_unique_chains_max": seq_unique_chains_max,
            "frac_alpha_min": frac_alpha_min,
            "frac_alpha_max": frac_alpha_max,
            "frac_beta_min": frac_beta_min,
            "frac_beta_max": frac_beta_max,
            "frac_coil_min": frac_coil_min,
            "frac_coil_max": frac_coil_max,
            "pore_circularity_min": pore_circularity_min,
            "pore_circularity_max": pore_circularity_max,
            "pore_uniformity_min": pore_uniformity_min,
            "pore_uniformity_max": pore_uniformity_max,
            "pocket_circularity_min": pocket_circularity_min,
            "pocket_circularity_max": pocket_circularity_max,
            "pocket_uniformity_min": pocket_uniformity_min,
            "pocket_uniformity_max": pocket_uniformity_max,
            "cavity_circularity_min": cavity_circularity_min,
            "cavity_circularity_max": cavity_circularity_max,
            "cavity_uniformity_min": cavity_uniformity_min,
            "cavity_uniformity_max": cavity_uniformity_max,
            "hub_circularity_min": hub_circularity_min,
            "hub_circularity_max": hub_circularity_max,
            "hub_uniformity_min": hub_uniformity_min,
            "hub_uniformity_max": hub_uniformity_max,
        },
        "sort_by": normalized_sort_by,
        "sort_dir": normalized_sort_dir,
        "limit": normalized_limit,
        "offset": normalized_offset,
        "total_count": total_count,
        "rows": row_dicts,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Query local gallery SQLite index with pore/structure filters."
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=Path("data") / "gallery.db",
        help="SQLite database path (default: data/gallery.db).",
    )
    parser.add_argument("--run-id", type=str, default=None, help="Optional run-id filter.")

    parser.add_argument("--pore-volume-min", type=float, default=None)
    parser.add_argument("--pore-volume-max", type=float, default=None)
    parser.add_argument("--pore-length-min", type=float, default=None)
    parser.add_argument("--pore-length-max", type=float, default=None)
    parser.add_argument("--pore-dmin-min", type=float, default=None)
    parser.add_argument("--pore-dmin-max", type=float, default=None)
    parser.add_argument("--pore-dmax-min", type=float, default=None)
    parser.add_argument("--pore-dmax-max", type=float, default=None)
    parser.add_argument("--num-pores-min", type=int, default=None)
    parser.add_argument("--num-pores-max", type=int, default=None)

    parser.add_argument("--pocket-volume-min", type=float, default=None)
    parser.add_argument("--pocket-volume-max", type=float, default=None)
    parser.add_argument("--pocket-length-min", type=float, default=None)
    parser.add_argument("--pocket-length-max", type=float, default=None)
    parser.add_argument("--pocket-dmin-min", type=float, default=None)
    parser.add_argument("--pocket-dmin-max", type=float, default=None)
    parser.add_argument("--pocket-dmax-min", type=float, default=None)
    parser.add_argument("--pocket-dmax-max", type=float, default=None)
    parser.add_argument("--num-pockets-min", type=int, default=None)
    parser.add_argument("--num-pockets-max", type=int, default=None)

    parser.add_argument("--cavity-volume-min", type=float, default=None)
    parser.add_argument("--cavity-volume-max", type=float, default=None)
    parser.add_argument("--cavity-length-min", type=float, default=None)
    parser.add_argument("--cavity-length-max", type=float, default=None)
    parser.add_argument("--cavity-dmin-min", type=float, default=None)
    parser.add_argument("--cavity-dmin-max", type=float, default=None)
    parser.add_argument("--cavity-dmax-min", type=float, default=None)
    parser.add_argument("--cavity-dmax-max", type=float, default=None)
    parser.add_argument("--num-cavities-min", type=int, default=None)
    parser.add_argument("--num-cavities-max", type=int, default=None)

    parser.add_argument("--hub-volume-min", type=float, default=None)
    parser.add_argument("--hub-volume-max", type=float, default=None)
    parser.add_argument("--hub-length-min", type=float, default=None)
    parser.add_argument("--hub-length-max", type=float, default=None)
    parser.add_argument("--hub-dmin-min", type=float, default=None)
    parser.add_argument("--hub-dmin-max", type=float, default=None)
    parser.add_argument("--hub-dmax-min", type=float, default=None)
    parser.add_argument("--hub-dmax-max", type=float, default=None)
    parser.add_argument("--num-hubs-min", type=int, default=None)
    parser.add_argument("--num-hubs-max", type=int, default=None)

    parser.add_argument("--chains-min", type=int, default=None)
    parser.add_argument("--chains-max", type=int, default=None)
    parser.add_argument("--residues-min", type=int, default=None)
    parser.add_argument("--residues-max", type=int, default=None)
    parser.add_argument("--seq-unique-chains-min", type=int, default=None)
    parser.add_argument("--seq-unique-chains-max", type=int, default=None)

    parser.add_argument("--frac-alpha-min", type=float, default=None)
    parser.add_argument("--frac-alpha-max", type=float, default=None)
    parser.add_argument("--frac-beta-min", type=float, default=None)
    parser.add_argument("--frac-beta-max", type=float, default=None)
    parser.add_argument("--frac-coil-min", type=float, default=None)
    parser.add_argument("--frac-coil-max", type=float, default=None)

    parser.add_argument("--pore-circularity-min", type=float, default=None)
    parser.add_argument("--pore-circularity-max", type=float, default=None)
    parser.add_argument("--pore-uniformity-min", type=float, default=None)
    parser.add_argument("--pore-uniformity-max", type=float, default=None)
    parser.add_argument("--pocket-circularity-min", type=float, default=None)
    parser.add_argument("--pocket-circularity-max", type=float, default=None)
    parser.add_argument("--pocket-uniformity-min", type=float, default=None)
    parser.add_argument("--pocket-uniformity-max", type=float, default=None)
    parser.add_argument("--cavity-circularity-min", type=float, default=None)
    parser.add_argument("--cavity-circularity-max", type=float, default=None)
    parser.add_argument("--cavity-uniformity-min", type=float, default=None)
    parser.add_argument("--cavity-uniformity-max", type=float, default=None)
    parser.add_argument("--hub-circularity-min", type=float, default=None)
    parser.add_argument("--hub-circularity-max", type=float, default=None)
    parser.add_argument("--hub-uniformity-min", type=float, default=None)
    parser.add_argument("--hub-uniformity-max", type=float, default=None)

    parser.add_argument("--limit", type=int, default=50)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument(
        "--sort-by",
        choices=sorted(_SORT_COLUMN_MAP.keys()),
        default=_DEFAULT_SORT_BY,
    )
    parser.add_argument(
        "--sort-dir",
        choices=("asc", "desc"),
        default=_DEFAULT_SORT_DIR,
    )
    parser.add_argument(
        "--compact",
        action="store_true",
        help="Emit compact JSON (single line) instead of pretty JSON.",
    )

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    result = query_gallery_index(
        db_path=Path(args.db),
        run_id=args.run_id,
        pore_volume_min=args.pore_volume_min,
        pore_volume_max=args.pore_volume_max,
        pore_length_min=args.pore_length_min,
        pore_length_max=args.pore_length_max,
        pore_dmin_min=args.pore_dmin_min,
        pore_dmin_max=args.pore_dmin_max,
        pore_dmax_min=args.pore_dmax_min,
        pore_dmax_max=args.pore_dmax_max,
        num_pores_min=args.num_pores_min,
        num_pores_max=args.num_pores_max,
        pocket_volume_min=args.pocket_volume_min,
        pocket_volume_max=args.pocket_volume_max,
        pocket_length_min=args.pocket_length_min,
        pocket_length_max=args.pocket_length_max,
        pocket_dmin_min=args.pocket_dmin_min,
        pocket_dmin_max=args.pocket_dmin_max,
        pocket_dmax_min=args.pocket_dmax_min,
        pocket_dmax_max=args.pocket_dmax_max,
        num_pockets_min=args.num_pockets_min,
        num_pockets_max=args.num_pockets_max,
        cavity_volume_min=args.cavity_volume_min,
        cavity_volume_max=args.cavity_volume_max,
        cavity_length_min=args.cavity_length_min,
        cavity_length_max=args.cavity_length_max,
        cavity_dmin_min=args.cavity_dmin_min,
        cavity_dmin_max=args.cavity_dmin_max,
        cavity_dmax_min=args.cavity_dmax_min,
        cavity_dmax_max=args.cavity_dmax_max,
        num_cavities_min=args.num_cavities_min,
        num_cavities_max=args.num_cavities_max,
        hub_volume_min=args.hub_volume_min,
        hub_volume_max=args.hub_volume_max,
        hub_length_min=args.hub_length_min,
        hub_length_max=args.hub_length_max,
        hub_dmin_min=args.hub_dmin_min,
        hub_dmin_max=args.hub_dmin_max,
        hub_dmax_min=args.hub_dmax_min,
        hub_dmax_max=args.hub_dmax_max,
        num_hubs_min=args.num_hubs_min,
        num_hubs_max=args.num_hubs_max,
        chains_min=args.chains_min,
        chains_max=args.chains_max,
        residues_min=args.residues_min,
        residues_max=args.residues_max,
        seq_unique_chains_min=args.seq_unique_chains_min,
        seq_unique_chains_max=args.seq_unique_chains_max,
        frac_alpha_min=args.frac_alpha_min,
        frac_alpha_max=args.frac_alpha_max,
        frac_beta_min=args.frac_beta_min,
        frac_beta_max=args.frac_beta_max,
        frac_coil_min=args.frac_coil_min,
        frac_coil_max=args.frac_coil_max,
        pore_circularity_min=args.pore_circularity_min,
        pore_circularity_max=args.pore_circularity_max,
        pore_uniformity_min=args.pore_uniformity_min,
        pore_uniformity_max=args.pore_uniformity_max,
        pocket_circularity_min=args.pocket_circularity_min,
        pocket_circularity_max=args.pocket_circularity_max,
        pocket_uniformity_min=args.pocket_uniformity_min,
        pocket_uniformity_max=args.pocket_uniformity_max,
        cavity_circularity_min=args.cavity_circularity_min,
        cavity_circularity_max=args.cavity_circularity_max,
        cavity_uniformity_min=args.cavity_uniformity_min,
        cavity_uniformity_max=args.cavity_uniformity_max,
        hub_circularity_min=args.hub_circularity_min,
        hub_circularity_max=args.hub_circularity_max,
        hub_uniformity_min=args.hub_uniformity_min,
        hub_uniformity_max=args.hub_uniformity_max,
        limit=args.limit,
        offset=args.offset,
        sort_by=args.sort_by,
        sort_dir=args.sort_dir,
    )

    if args.compact:
        print(json.dumps(result, separators=(",", ":")))
    else:
        print(json.dumps(result, indent=2))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
