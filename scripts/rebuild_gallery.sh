#!/usr/bin/env bash
set -euo pipefail

DB_PATH="${VOLUMIZER_GALLERY_DB:-data/gallery.db}"
RAW_PATH="${1:?Usage: $0 <path/to/run-dir-or-summary.json> [run-id]}"
RUN_ID="${2:-}"

# Accept either a directory or a file
if [[ -d "$RAW_PATH" ]]; then
  SUMMARY_PATH="$RAW_PATH/run.summary.json"
  if [[ ! -f "$SUMMARY_PATH" ]]; then
    echo "Error: No run.summary.json found in $RAW_PATH" >&2
    exit 1
  fi
elif [[ -f "$RAW_PATH" ]]; then
  SUMMARY_PATH="$RAW_PATH"
else
  echo "Error: Path does not exist: $RAW_PATH" >&2
  exit 1
fi

RUN_ID_ARGS=()
if [[ -n "$RUN_ID" ]]; then
  RUN_ID_ARGS=(--run-id "$RUN_ID")
fi

echo "==> Removing existing gallery DB: $DB_PATH"
rm -f "$DB_PATH"

echo "==> Rebuilding gallery index from: $SUMMARY_PATH"
uv run python -m volumizer.gallery_index --summary "$SUMMARY_PATH" --db "$DB_PATH" "${RUN_ID_ARGS[@]}"

echo "==> Rendering thumbnails"
uv run python -m volumizer.gallery_render --db "$DB_PATH" "${RUN_ID_ARGS[@]}"

echo "==> Done"
