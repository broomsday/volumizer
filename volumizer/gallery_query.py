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
    "num_pores": "a.num_pores",
    "largest_pore_volume": "a.largest_pore_volume_a3",
    "largest_pore_length": "a.largest_pore_length_a",
    "largest_pore_min_diameter": "a.largest_pore_min_diameter_a",
    "largest_pore_max_diameter": "a.largest_pore_max_diameter_a",
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
    chains_min: int | None = None,
    chains_max: int | None = None,
    residues_min: int | None = None,
    residues_max: int | None = None,
    seq_unique_chains_min: int | None = None,
    seq_unique_chains_max: int | None = None,
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

    _append_range_filter(where_clauses, params, "s.num_chains", chains_min, chains_max)
    _append_range_filter(where_clauses, params, "s.num_residues", residues_min, residues_max)
    _append_range_filter(
        where_clauses,
        params,
        "s.num_sequence_unique_chains",
        seq_unique_chains_min,
        seq_unique_chains_max,
    )

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
        "a.num_pores, a.largest_pore_volume_a3, a.largest_pore_length_a, "
        "a.largest_pore_min_diameter_a, a.largest_pore_max_diameter_a, "
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
            "chains_min": chains_min,
            "chains_max": chains_max,
            "residues_min": residues_min,
            "residues_max": residues_max,
            "seq_unique_chains_min": seq_unique_chains_min,
            "seq_unique_chains_max": seq_unique_chains_max,
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

    parser.add_argument("--chains-min", type=int, default=None)
    parser.add_argument("--chains-max", type=int, default=None)
    parser.add_argument("--residues-min", type=int, default=None)
    parser.add_argument("--residues-max", type=int, default=None)
    parser.add_argument("--seq-unique-chains-min", type=int, default=None)
    parser.add_argument("--seq-unique-chains-max", type=int, default=None)

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
        chains_min=args.chains_min,
        chains_max=args.chains_max,
        residues_min=args.residues_min,
        residues_max=args.residues_max,
        seq_unique_chains_min=args.seq_unique_chains_min,
        seq_unique_chains_max=args.seq_unique_chains_max,
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
