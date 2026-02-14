"""
Helpers for comparing volumizer annotation outputs across backends.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = ("id", "type", "volume", "x", "y", "z")
NUMERIC_COLUMNS = ("volume", "x", "y", "z")


def load_annotation_json(path: Path) -> pd.DataFrame:
    """
    Load an annotation dataframe from JSON and normalize expected columns.
    """
    return normalize_annotation_dataframe(pd.read_json(path))


def normalize_annotation_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return a stable representation suitable for parity comparisons.
    """
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise AssertionError(f"Missing expected annotation columns: {missing}")

    normalized = df.loc[:, REQUIRED_COLUMNS].copy()
    normalized["id"] = pd.to_numeric(normalized["id"], errors="raise").astype(int)
    normalized["type"] = normalized["type"].astype(str)
    for column in NUMERIC_COLUMNS:
        normalized[column] = pd.to_numeric(normalized[column], errors="raise").astype(
            float
        )

    normalized = normalized.sort_values(
        by=["type", "volume", "x", "y", "z", "id"],
        ascending=[True, False, False, False, False, True],
    ).reset_index(drop=True)
    return normalized


def assert_annotation_dataframes_close(
    actual: pd.DataFrame,
    expected: pd.DataFrame,
    volume_tolerance: float = 1e-6,
    dimension_tolerance: float = 1e-2,
) -> None:
    """
    Assert two annotation dataframes are equal within explicit tolerances.
    """
    actual_normalized = normalize_annotation_dataframe(actual)
    expected_normalized = normalize_annotation_dataframe(expected)

    if len(actual_normalized) != len(expected_normalized):
        raise AssertionError(
            "Annotation row count mismatch: "
            f"actual={len(actual_normalized)} expected={len(expected_normalized)}"
        )

    mismatches: list[str] = []
    for idx in range(len(actual_normalized)):
        actual_row = actual_normalized.iloc[idx]
        expected_row = expected_normalized.iloc[idx]

        if actual_row["type"] != expected_row["type"]:
            mismatches.append(
                f"row {idx} type mismatch: actual={actual_row['type']} expected={expected_row['type']}"
            )

        volume_delta = abs(actual_row["volume"] - expected_row["volume"])
        if volume_delta > volume_tolerance:
            mismatches.append(
                f"row {idx} volume mismatch: actual={actual_row['volume']} "
                f"expected={expected_row['volume']} delta={volume_delta}"
            )

        for column in ("x", "y", "z"):
            delta = abs(actual_row[column] - expected_row[column])
            if delta > dimension_tolerance:
                mismatches.append(
                    f"row {idx} {column} mismatch: actual={actual_row[column]} "
                    f"expected={expected_row[column]} delta={delta}"
                )

    if mismatches:
        raise AssertionError("\n".join(mismatches))
