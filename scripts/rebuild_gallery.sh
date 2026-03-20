#!/usr/bin/env bash
set -euo pipefail

DB_PATH="${VOLUMIZER_GALLERY_DB:-data/gallery.db}"
SUMMARY_PATH="${1:?Usage: $0 <path/to/run.summary.json> [run-id]}"
RUN_ID="${2:-}"

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
