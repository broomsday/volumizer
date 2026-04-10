#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_VERSION="3.11"
WITH_GALLERY=true
WITH_NATIVE=false
COMPILE_C_HELPERS=true
PLAYWRIGHT_WITH_DEPS=false

log() {
  echo "[bootstrap-linux] $*" >&2
}

fail() {
  log "error: $*"
  exit 1
}

usage() {
  cat <<EOF
Usage: bash scripts/bootstrap_linux.sh [options]

Prepare the repository for local use on Linux after cloning.

Options:
  --python-version VERSION     Python version for uv (default: 3.11)
  --no-gallery                 Skip npm/Playwright gallery setup
  --with-native                Build the optional Rust native backend
  --skip-c-helpers             Skip optional local C helper compilation
  --playwright-with-deps       Run 'playwright install --with-deps chromium'
  -h, --help                   Show this help message
EOF
}

require_command() {
  local command_name="$1"
  if ! command -v "$command_name" >/dev/null 2>&1; then
    fail "missing required command: $command_name"
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --python-version)
      PYTHON_VERSION="${2:?--python-version requires a value}"
      shift 2
      ;;
    --no-gallery)
      WITH_GALLERY=false
      shift
      ;;
    --with-native)
      WITH_NATIVE=true
      shift
      ;;
    --skip-c-helpers)
      COMPILE_C_HELPERS=false
      shift
      ;;
    --playwright-with-deps)
      PLAYWRIGHT_WITH_DEPS=true
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      fail "unknown argument: $1"
      ;;
  esac
done

require_command uv

if [[ "$WITH_GALLERY" == true ]]; then
  require_command node
  require_command npm
fi

if [[ "$WITH_NATIVE" == true ]]; then
  require_command cargo
  require_command rustc
fi

cd "$REPO_ROOT"

log "installing Python ${PYTHON_VERSION} via uv"
uv python install "$PYTHON_VERSION"

UV_SYNC_ARGS=(
  --python "$PYTHON_VERSION"
  --group test
)
if [[ "$WITH_GALLERY" == true ]]; then
  UV_SYNC_ARGS+=(--group web)
fi
if [[ "$WITH_NATIVE" == true ]]; then
  UV_SYNC_ARGS+=(--group native)
fi

log "syncing Python dependencies"
uv sync "${UV_SYNC_ARGS[@]}"

if [[ "$COMPILE_C_HELPERS" == true ]]; then
  if command -v cc >/dev/null 2>&1; then
    log "building optional local C helpers"
    bash src/compile_c_libs.sh
  else
    log "warning: 'cc' not found; skipping optional C helpers"
  fi
fi

if [[ "$WITH_GALLERY" == true ]]; then
  if [[ -f package-lock.json ]]; then
    log "installing npm dependencies with npm ci"
    npm ci
  else
    log "installing npm dependencies with npm install"
    npm install
  fi

  if [[ "$PLAYWRIGHT_WITH_DEPS" == true ]]; then
    log "installing Playwright Chromium with system dependencies"
    npx playwright install --with-deps chromium
  else
    log "installing Playwright Chromium"
    npx playwright install chromium
  fi
fi

if [[ "$WITH_NATIVE" == true ]]; then
  log "building optional Rust native backend"
  uv run --python "$PYTHON_VERSION" maturin develop --manifest-path native/Cargo.toml
fi

log "bootstrap complete"
cat <<EOF

Next steps:
  Run the CLI:
    uv run --python ${PYTHON_VERSION} volumizer --version

  Run the gallery:
    ./gallery /path/to/run.summary.json

  Serve the gallery without rendering thumbnails first:
    ./gallery /path/to/run.summary.json --skip-thumbnails
EOF
