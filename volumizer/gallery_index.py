"""
Build a local SQLite index from volumizer analyze/cluster outputs for gallery filtering.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sqlite3
import sys
from typing import Any

import numpy as np

from volumizer import pdb, rcsb, utils


_VALID_VOLUME_KINDS = {"pore", "pocket", "cavity", "hub"}


def _safe_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _normalize_pdb_id(value: Any) -> str | None:
    if not isinstance(value, str):
        return None

    candidate = value.strip()
    if len(candidate) == 0:
        return None

    try:
        return rcsb.normalize_pdb_id(candidate)
    except ValueError:
        return None


def _infer_pdb_id(result: dict[str, Any], source_label: str, input_path: Path | None) -> str | None:
    direct = _normalize_pdb_id(result.get("pdb_id"))
    if direct is not None:
        return direct

    from_source = _normalize_pdb_id(source_label)
    if from_source is not None:
        return from_source

    if input_path is not None:
        from_input = _normalize_pdb_id(input_path.stem)
        if from_input is not None:
            return from_input

    return None


def _sort_row_key(value: str) -> tuple[int, str]:
    try:
        return (0, f"{int(value):09d}")
    except (TypeError, ValueError):
        return (1, str(value))


def _warn(message: str) -> None:
    print(f"[gallery-index] {message}", file=sys.stderr)


def _dedupe_paths(paths: list[Path]) -> list[Path]:
    deduped: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        deduped.append(path)
    return deduped


def _resolve_path(
    raw_path: str | None,
    summary_path: Path,
    output_dir_hint: Path,
) -> Path | None:
    if raw_path is None:
        return None

    path = Path(str(raw_path))
    if path.is_absolute():
        return path

    summary_dir = summary_path.parent
    candidates = [
        output_dir_hint / path,
        summary_dir / path,
        summary_dir.parent / path,
        Path.cwd() / path,
    ]

    if len(path.parts) > 1 and path.parts[0] == output_dir_hint.name:
        candidates.append(output_dir_hint / Path(*path.parts[1:]))
    if len(path.parts) > 1 and path.parts[0] == summary_dir.name:
        candidates.append(summary_dir / Path(*path.parts[1:]))

    candidates = _dedupe_paths(candidates)
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    return candidates[0].resolve()


def _records_from_dataframe_json(payload: dict[str, Any]) -> list[dict[str, Any]]:
    required_columns = ("id", "type", "volume", "x", "y", "z")
    if any(column not in payload for column in required_columns):
        return []

    if any(not isinstance(payload[column], dict) for column in required_columns):
        return []

    row_keys = sorted(
        set(payload["type"].keys()),
        key=_sort_row_key,
    )

    records: list[dict[str, Any]] = []
    for row_key in row_keys:
        records.append(
            {
                "id": payload["id"].get(row_key),
                "type": payload["type"].get(row_key),
                "volume": payload["volume"].get(row_key),
                "x": payload["x"].get(row_key),
                "y": payload["y"].get(row_key),
                "z": payload["z"].get(row_key),
            }
        )

    return records


def _parse_volume_records(annotation_payload: Any) -> list[dict[str, Any]]:
    if isinstance(annotation_payload, dict):
        volumes = annotation_payload.get("volumes")
        if isinstance(volumes, list):
            return [volume for volume in volumes if isinstance(volume, dict)]
        return _records_from_dataframe_json(annotation_payload)

    if isinstance(annotation_payload, list):
        return [volume for volume in annotation_payload if isinstance(volume, dict)]

    return []


def _normalize_volume_rows(
    volume_records: list[dict[str, Any]],
) -> dict[str, list[dict[str, float | int | str | None]]]:
    grouped: dict[str, list[dict[str, float | int | str | None]]] = {
        "pore": [],
        "pocket": [],
        "cavity": [],
        "hub": [],
    }

    for raw_row in volume_records:
        kind = str(raw_row.get("type", "")).strip().lower()
        if kind not in _VALID_VOLUME_KINDS:
            continue

        volume_value = _safe_float(raw_row.get("volume"))
        if volume_value is None:
            continue

        axial_lengths = [
            _safe_float(raw_row.get("x")) or 0.0,
            _safe_float(raw_row.get("y")) or 0.0,
            _safe_float(raw_row.get("z")) or 0.0,
        ]
        axial_lengths = sorted(axial_lengths, reverse=True)

        grouped[kind].append(
            {
                "kind": kind,
                "raw_id": raw_row.get("id"),
                "volume_a3": volume_value,
                "length_a": float(axial_lengths[0]),
                "max_diameter_a": float(axial_lengths[1]),
                "min_diameter_a": float(axial_lengths[2]),
                "centroid_x": None,
                "centroid_y": None,
                "centroid_z": None,
                "cross_section_circularity": _safe_float(raw_row.get("cross_section_circularity")),
                "cross_section_uniformity": _safe_float(raw_row.get("cross_section_uniformity")),
            }
        )

    for kind in grouped:
        grouped[kind].sort(key=lambda row: float(row["volume_a3"]), reverse=True)

    return grouped


def _extract_chain_sequences(structure) -> dict[str, tuple[str, ...]]:
    if len(structure) == 0:
        return {}

    chain_ids = np.asarray(structure.chain_id, dtype=object)
    residue_ids = np.asarray(structure.res_id, dtype=object)
    residue_names = np.asarray(structure.res_name, dtype=object)
    if hasattr(structure, "ins_code"):
        insertion_codes = np.asarray(structure.ins_code, dtype=object)
    else:
        insertion_codes = np.full(len(structure), "", dtype=object)

    chain_sequences: dict[str, list[str]] = {}
    previous_residue_key: dict[str, tuple[str, str]] = {}

    for chain_id, residue_id, insertion_code, residue_name in zip(
        chain_ids,
        residue_ids,
        insertion_codes,
        residue_names,
    ):
        chain_key = str(chain_id)
        residue_key = (str(residue_id), str(insertion_code))
        if previous_residue_key.get(chain_key) == residue_key:
            continue

        previous_residue_key[chain_key] = residue_key
        chain_sequences.setdefault(chain_key, []).append(str(residue_name))

    return {
        chain_key: tuple(sequence)
        for chain_key, sequence in chain_sequences.items()
    }


def _sequence_overlap(seq1: tuple[str, ...], seq2: tuple[str, ...]) -> float:
    """
    Compute fractional overlap between two residue-name sequences.

    Uses best-offset alignment to handle terminal truncation: slides the
    shorter sequence along the longer one and returns the best match ratio
    relative to the longer sequence length.
    """
    len1, len2 = len(seq1), len(seq2)
    max_len = max(len1, len2)
    if max_len == 0:
        return 1.0

    # Ensure seq1 is the longer one
    if len1 < len2:
        seq1, seq2 = seq2, seq1
        len1, len2 = len2, len1

    best_matches = 0
    # Slide seq2 along seq1 at each possible offset
    for offset in range(len1 - len2 + 1):
        matches = sum(1 for i in range(len2) if seq1[offset + i] == seq2[i])
        if matches > best_matches:
            best_matches = matches
            if best_matches == len2:
                break

    return best_matches / max_len


def _count_sequence_unique_chains(
    chain_sequences: dict[str, tuple[str, ...]],
    identity_threshold: float = 0.95,
) -> int:
    """
    Count sequence-unique chains using fuzzy matching.

    Chains with >= identity_threshold fractional overlap are considered
    the same, which handles minor terminal truncation differences between
    symmetry copies.
    """
    sequences = list(chain_sequences.values())
    if len(sequences) == 0:
        return 0

    # Greedy clustering: assign each sequence to the first matching cluster
    cluster_reps: list[tuple[str, ...]] = [sequences[0]]
    for seq in sequences[1:]:
        matched = False
        for rep in cluster_reps:
            if _sequence_overlap(seq, rep) >= identity_threshold:
                matched = True
                break
        if not matched:
            cluster_reps.append(seq)

    return len(cluster_reps)


def _compute_structure_metrics(
    input_path: Path,
    assembly_policy: str,
) -> tuple[int | None, int | None, int | None]:
    try:
        structure = pdb.load_structure(input_path, assembly_policy=assembly_policy)
        cleaned = pdb.clean_structure(structure)
    except Exception as error:  # pragma: no cover - defensive
        _warn(f"failed to load/clean {input_path}: {error}")
        return None, None, None

    chain_sequences = _extract_chain_sequences(cleaned)
    if len(chain_sequences) == 0:
        return 0, 0, 0

    num_chains = len(chain_sequences)
    num_residues = int(sum(len(sequence) for sequence in chain_sequences.values()))
    num_sequence_unique_chains = _count_sequence_unique_chains(chain_sequences)

    return num_chains, num_residues, num_sequence_unique_chains


_SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS runs (
    run_id TEXT PRIMARY KEY,
    created_at TEXT NOT NULL,
    source_summary_path TEXT NOT NULL,
    resolution REAL,
    assembly_policy TEXT
);

CREATE TABLE IF NOT EXISTS structures (
    structure_id INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id TEXT NOT NULL,
    source_label TEXT NOT NULL,
    pdb_id TEXT,
    input_path TEXT,
    annotated_cif_path TEXT,
    annotation_json_path TEXT,
    num_chains INTEGER,
    num_residues INTEGER,
    num_sequence_unique_chains INTEGER,
    frac_alpha REAL,
    frac_beta REAL,
    frac_coil REAL,
    UNIQUE(run_id, source_label),
    FOREIGN KEY(run_id) REFERENCES runs(run_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS volumes (
    volume_id INTEGER PRIMARY KEY AUTOINCREMENT,
    structure_id INTEGER NOT NULL,
    kind TEXT NOT NULL,
    rank_in_kind INTEGER NOT NULL,
    volume_a3 REAL NOT NULL,
    length_a REAL,
    min_diameter_a REAL,
    max_diameter_a REAL,
    centroid_x REAL,
    centroid_y REAL,
    centroid_z REAL,
    cross_section_circularity REAL,
    cross_section_uniformity REAL,
    FOREIGN KEY(structure_id) REFERENCES structures(structure_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS structure_aggregates (
    structure_id INTEGER PRIMARY KEY,
    num_pores INTEGER NOT NULL DEFAULT 0,
    largest_pore_volume_a3 REAL,
    largest_pore_length_a REAL,
    largest_pore_min_diameter_a REAL,
    largest_pore_max_diameter_a REAL,
    num_pockets INTEGER NOT NULL DEFAULT 0,
    largest_pocket_volume_a3 REAL,
    largest_pocket_length_a REAL,
    largest_pocket_min_diameter_a REAL,
    largest_pocket_max_diameter_a REAL,
    num_cavities INTEGER NOT NULL DEFAULT 0,
    largest_cavity_volume_a3 REAL,
    largest_cavity_length_a REAL,
    largest_cavity_min_diameter_a REAL,
    largest_cavity_max_diameter_a REAL,
    num_hubs INTEGER NOT NULL DEFAULT 0,
    largest_hub_volume_a3 REAL,
    largest_hub_length_a REAL,
    largest_hub_min_diameter_a REAL,
    largest_hub_max_diameter_a REAL,
    largest_pore_circularity REAL,
    largest_pore_uniformity REAL,
    largest_pocket_circularity REAL,
    largest_pocket_uniformity REAL,
    largest_cavity_circularity REAL,
    largest_cavity_uniformity REAL,
    largest_hub_circularity REAL,
    largest_hub_uniformity REAL,
    FOREIGN KEY(structure_id) REFERENCES structures(structure_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS renders (
    structure_id INTEGER PRIMARY KEY,
    x_png_path TEXT,
    y_png_path TEXT,
    z_png_path TEXT,
    render_style_hash TEXT,
    render_status TEXT NOT NULL,
    render_error TEXT,
    updated_at TEXT NOT NULL,
    FOREIGN KEY(structure_id) REFERENCES structures(structure_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_structures_filters
ON structures (num_chains, num_residues, num_sequence_unique_chains);

CREATE INDEX IF NOT EXISTS idx_structures_sse
ON structures (frac_alpha, frac_beta, frac_coil);

CREATE INDEX IF NOT EXISTS idx_structure_aggregates_pore
ON structure_aggregates (
    num_pores,
    largest_pore_volume_a3,
    largest_pore_length_a,
    largest_pore_min_diameter_a,
    largest_pore_max_diameter_a
);

CREATE INDEX IF NOT EXISTS idx_structure_aggregates_pocket
ON structure_aggregates (
    num_pockets,
    largest_pocket_volume_a3,
    largest_pocket_length_a,
    largest_pocket_min_diameter_a,
    largest_pocket_max_diameter_a
);

CREATE INDEX IF NOT EXISTS idx_structure_aggregates_cavity
ON structure_aggregates (
    num_cavities,
    largest_cavity_volume_a3,
    largest_cavity_length_a,
    largest_cavity_min_diameter_a,
    largest_cavity_max_diameter_a
);

CREATE INDEX IF NOT EXISTS idx_structure_aggregates_hub
ON structure_aggregates (
    num_hubs,
    largest_hub_volume_a3,
    largest_hub_length_a,
    largest_hub_min_diameter_a,
    largest_hub_max_diameter_a
);

CREATE INDEX IF NOT EXISTS idx_volumes_structure_kind
ON volumes (structure_id, kind, rank_in_kind);
"""


def _derive_default_run_id(summary_path: Path) -> str:
    candidate = summary_path.parent.name.strip()
    if len(candidate) == 0:
        candidate = summary_path.stem

    sanitized = "".join(
        char.lower() if (char.isalnum() or char in {"-", "_"}) else "_"
        for char in candidate
    ).strip("_")
    if len(sanitized) == 0:
        return "run"
    return sanitized


def build_gallery_index(
    summary_path: Path,
    db_path: Path,
    run_id: str | None = None,
    replace_run: bool = False,
    strict: bool = False,
) -> dict[str, int | str]:
    summary_path = summary_path.resolve()
    if not summary_path.is_file():
        raise FileNotFoundError(f"Summary file does not exist: {summary_path}")

    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("Summary payload must be a JSON object")

    config = payload.get("config")
    if not isinstance(config, dict):
        config = {}

    output_dir_hint = summary_path.parent
    config_output_dir = config.get("output_dir")
    if isinstance(config_output_dir, str):
        output_dir_hint = _resolve_path(config_output_dir, summary_path, summary_path.parent) or summary_path.parent

    assembly_policy = str(config.get("assembly_policy") or pdb.DEFAULT_ASSEMBLY_POLICY)
    resolution = _safe_float(config.get("resolution"))

    requested_keep_non_protein = bool(config.get("keep_non_protein", False))
    previous_keep_non_protein = bool(utils.KEEP_NON_PROTEIN)
    utils.set_non_protein(requested_keep_non_protein)

    summary_results = payload.get("results")
    if not isinstance(summary_results, list):
        raise ValueError("Summary payload is missing list field: results")

    db_path.parent.mkdir(parents=True, exist_ok=True)

    effective_run_id = run_id if run_id is not None else _derive_default_run_id(summary_path)
    if len(str(effective_run_id).strip()) == 0:
        raise ValueError("run_id must not be empty")

    now_utc = datetime.now(tz=timezone.utc).isoformat()
    indexed_structures = 0
    indexed_volumes = 0
    skipped_structures = 0

    input_metric_cache: dict[Path, tuple[int | None, int | None, int | None]] = {}

    try:
        with sqlite3.connect(db_path) as connection:
            connection.executescript(_SCHEMA_SQL)

            existing_run = connection.execute(
                "SELECT 1 FROM runs WHERE run_id = ?",
                (effective_run_id,),
            ).fetchone()
            if existing_run is not None and not replace_run:
                raise ValueError(
                    f"Run already exists in index: {effective_run_id}. "
                    "Use --replace-run to rebuild this run ID."
                )
            if existing_run is not None and replace_run:
                connection.execute(
                    "DELETE FROM runs WHERE run_id = ?",
                    (effective_run_id,),
                )

            connection.execute(
                """
                INSERT INTO runs (run_id, created_at, source_summary_path, resolution, assembly_policy)
                VALUES (?, ?, ?, ?, ?)
                """,
                (
                    effective_run_id,
                    now_utc,
                    str(summary_path),
                    resolution,
                    assembly_policy,
                ),
            )

            for result in summary_results:
                if not isinstance(result, dict):
                    skipped_structures += 1
                    continue

                source_label = str(result.get("source") or "structure")
                input_path = _resolve_path(
                    result.get("input_path"),
                    summary_path,
                    output_dir_hint,
                )
                annotation_path = _resolve_path(
                    result.get("annotation_output"),
                    summary_path,
                    output_dir_hint,
                )
                annotated_cif_path = _resolve_path(
                    result.get("structure_output"),
                    summary_path,
                    output_dir_hint,
                )

                if annotation_path is None or not annotation_path.is_file():
                    skipped_structures += 1
                    _warn(
                        f"skipping {source_label}: annotation output missing ({annotation_path})"
                    )
                    if strict:
                        raise FileNotFoundError(
                            f"Missing annotation output for {source_label}: {annotation_path}"
                        )
                    continue

                num_chains = None
                num_residues = None
                num_sequence_unique_chains = None
                if input_path is None or not input_path.is_file():
                    _warn(
                        f"structure metrics unavailable for {source_label}: input path missing ({input_path})"
                    )
                    if strict:
                        raise FileNotFoundError(
                            f"Missing input structure for {source_label}: {input_path}"
                        )
                else:
                    if input_path not in input_metric_cache:
                        input_metric_cache[input_path] = _compute_structure_metrics(
                            input_path,
                            assembly_policy=assembly_policy,
                        )
                    (
                        num_chains,
                        num_residues,
                        num_sequence_unique_chains,
                    ) = input_metric_cache[input_path]

                annotation_payload = json.loads(annotation_path.read_text(encoding="utf-8"))
                frac_alpha = _safe_float(
                    annotation_payload.get("frac_alpha") if isinstance(annotation_payload, dict) else None
                )
                frac_beta = _safe_float(
                    annotation_payload.get("frac_beta") if isinstance(annotation_payload, dict) else None
                )
                frac_coil = _safe_float(
                    annotation_payload.get("frac_coil") if isinstance(annotation_payload, dict) else None
                )
                volume_records = _parse_volume_records(annotation_payload)
                grouped_rows = _normalize_volume_rows(volume_records)
                pdb_id = _infer_pdb_id(
                    result,
                    source_label=source_label,
                    input_path=input_path,
                )

                cursor = connection.execute(
                    """
                    INSERT INTO structures (
                        run_id,
                        source_label,
                        pdb_id,
                        input_path,
                        annotated_cif_path,
                        annotation_json_path,
                        num_chains,
                        num_residues,
                        num_sequence_unique_chains,
                        frac_alpha,
                        frac_beta,
                        frac_coil
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        effective_run_id,
                        source_label,
                        pdb_id,
                        str(input_path) if input_path is not None else None,
                        str(annotated_cif_path) if annotated_cif_path is not None else None,
                        str(annotation_path),
                        num_chains,
                        num_residues,
                        num_sequence_unique_chains,
                        frac_alpha,
                        frac_beta,
                        frac_coil,
                    ),
                )
                structure_id = int(cursor.lastrowid)

                for kind, rows in grouped_rows.items():
                    for rank_in_kind, row in enumerate(rows, start=1):
                        connection.execute(
                            """
                            INSERT INTO volumes (
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
                            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                            """,
                            (
                                structure_id,
                                kind,
                                rank_in_kind,
                                float(row["volume_a3"]),
                                row["length_a"],
                                row["min_diameter_a"],
                                row["max_diameter_a"],
                                row["centroid_x"],
                                row["centroid_y"],
                                row["centroid_z"],
                                row["cross_section_circularity"],
                                row["cross_section_uniformity"],
                            ),
                        )
                        indexed_volumes += 1

                def _kind_aggregates(kind: str) -> tuple[int, float | None, float | None, float | None, float | None, float | None, float | None]:
                    rows = grouped_rows[kind]
                    if len(rows) > 0:
                        largest = rows[0]
                        return (
                            len(rows),
                            float(largest["volume_a3"]),
                            float(largest["length_a"]),
                            float(largest["min_diameter_a"]),
                            float(largest["max_diameter_a"]),
                            largest["cross_section_circularity"],
                            largest["cross_section_uniformity"],
                        )
                    return (0, None, None, None, None, None, None)

                pore_agg = _kind_aggregates("pore")
                pocket_agg = _kind_aggregates("pocket")
                cavity_agg = _kind_aggregates("cavity")
                hub_agg = _kind_aggregates("hub")

                connection.execute(
                    """
                    INSERT INTO structure_aggregates (
                        structure_id,
                        num_pores, largest_pore_volume_a3, largest_pore_length_a,
                        largest_pore_min_diameter_a, largest_pore_max_diameter_a,
                        largest_pore_circularity, largest_pore_uniformity,
                        num_pockets, largest_pocket_volume_a3, largest_pocket_length_a,
                        largest_pocket_min_diameter_a, largest_pocket_max_diameter_a,
                        largest_pocket_circularity, largest_pocket_uniformity,
                        num_cavities, largest_cavity_volume_a3, largest_cavity_length_a,
                        largest_cavity_min_diameter_a, largest_cavity_max_diameter_a,
                        largest_cavity_circularity, largest_cavity_uniformity,
                        num_hubs, largest_hub_volume_a3, largest_hub_length_a,
                        largest_hub_min_diameter_a, largest_hub_max_diameter_a,
                        largest_hub_circularity, largest_hub_uniformity
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                              ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        structure_id,
                        *pore_agg,
                        *pocket_agg,
                        *cavity_agg,
                        *hub_agg,
                    ),
                )

                connection.execute(
                    """
                    INSERT INTO renders (
                        structure_id,
                        x_png_path,
                        y_png_path,
                        z_png_path,
                        render_style_hash,
                        render_status,
                        render_error,
                        updated_at
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        structure_id,
                        None,
                        None,
                        None,
                        None,
                        "pending",
                        None,
                        now_utc,
                    ),
                )

                indexed_structures += 1

            connection.commit()
    finally:
        utils.set_non_protein(previous_keep_non_protein)

    return {
        "run_id": effective_run_id,
        "indexed_structures": indexed_structures,
        "indexed_volumes": indexed_volumes,
        "skipped_structures": skipped_structures,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build or refresh local gallery SQLite index from volumizer run.summary.json"
        )
    )
    parser.add_argument(
        "--summary",
        type=Path,
        required=True,
        help="Path to run.summary.json generated by volumizer CLI.",
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=Path("data") / "gallery.db",
        help="SQLite database path (default: data/gallery.db).",
    )
    parser.add_argument(
        "--run-id",
        type=str,
        default=None,
        help="Optional run identifier override (default: derived from summary parent dir).",
    )
    parser.add_argument(
        "--replace-run",
        action="store_true",
        help="Replace existing rows for this run-id if present.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail on missing input/annotation paths instead of skipping entries.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    result = build_gallery_index(
        summary_path=Path(args.summary),
        db_path=Path(args.db),
        run_id=args.run_id,
        replace_run=bool(args.replace_run),
        strict=bool(args.strict),
    )

    print(
        json.dumps(
            {
                "summary": str(Path(args.summary).resolve()),
                "db": str(Path(args.db).resolve()),
                **result,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
