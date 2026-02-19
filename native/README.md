# Native Module

This directory contains the Rust extension module `volumizer_native` used by the optional native backend.

## Implemented Kernels

Current native kernels include:
- Fibonacci sphere point generation
- Batch atom-shell expansion
- Neighbor voxel detection
- BFS component expansion
- Exposed/buried solvent split
- Buried-component classification output

## Backend Selection

Runtime backend mode is controlled by `VOLUMIZER_BACKEND`:
- `python`: force Python/ctypes path
- `auto`: use native when importable, otherwise fallback to Python
- `native`: require native module and fail if unavailable

## Local Build

From repository root:

```bash
uv sync --python 3.11 --group test --group native
uv run --python 3.11 maturin develop --manifest-path native/Cargo.toml
```

Then verify import:

```bash
uv run --python 3.11 python -c "import volumizer_native; print(volumizer_native.backend_info())"
```

## Validation

Run the full test suite, or at minimum native-focused tests:

```bash
uv run --python 3.11 pytest -q tests/test_native_backend.py tests/test_voxel_native_backend.py
```
