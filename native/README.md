# Native Module Scaffold

This directory contains the Phase 1 Rust scaffold for `volumizer_native`.

## Purpose

- Establish a buildable `pyo3` extension module.
- Provide a first contract-aligned function (`fibonacci_sphere_points`).
- Keep Python fallback as default until native parity coverage is expanded.

## Local Build (when Rust toolchain is installed)

From repository root:

```bash
uv sync --python 3.11 --group test
uv run --python 3.11 pip install maturin
uv run --python 3.11 maturin develop --manifest-path native/Cargo.toml
```

Then validate import:

```bash
uv run --python 3.11 python -c "import volumizer_native; print(volumizer_native.backend_info())"
```
