"""
Utilities for selecting and loading the optional native Rust backend.
"""

from __future__ import annotations

from functools import lru_cache
import importlib
import importlib.util
import os
from pathlib import Path
import sys
from types import ModuleType


BACKEND_ENV = "VOLUMIZER_BACKEND"
BACKEND_AUTO = "auto"
BACKEND_NATIVE = "native"
BACKEND_PYTHON = "python"
VALID_BACKENDS = {BACKEND_AUTO, BACKEND_NATIVE, BACKEND_PYTHON}
DEFAULT_BACKEND = BACKEND_AUTO


def _get_local_native_artifact_candidates(root_dir: Path) -> list[Path]:
    """
    Return existing local build artifacts ordered by newest mtime first.
    """
    candidates = [
        root_dir / "native" / "target" / "debug" / "libvolumizer_native.so",
        root_dir / "native" / "target" / "release" / "libvolumizer_native.so",
        root_dir / "native" / "target" / "debug" / "volumizer_native.pyd",
        root_dir / "native" / "target" / "release" / "volumizer_native.pyd",
        root_dir / "native" / "target" / "debug" / "libvolumizer_native.dylib",
        root_dir / "native" / "target" / "release" / "libvolumizer_native.dylib",
    ]

    existing_candidates = [candidate for candidate in candidates if candidate.is_file()]
    return sorted(
        existing_candidates,
        key=lambda candidate: candidate.stat().st_mtime,
        reverse=True,
    )


def _get_native_source_mtime(root_dir: Path) -> float:
    """
    Return the newest mtime across native build inputs.
    """
    source_paths = [
        root_dir / "native" / "Cargo.toml",
        root_dir / "native" / "Cargo.lock",
    ]
    source_paths.extend((root_dir / "native" / "src").rglob("*.rs"))

    existing_source_paths = [path for path in source_paths if path.is_file()]
    if not existing_source_paths:
        return 0.0

    return max(path.stat().st_mtime for path in existing_source_paths)


def _find_stale_local_native_artifact() -> Path | None:
    """
    Return the newest stale local artifact when one exists.
    """
    root_dir = Path(__file__).resolve().parents[1]
    source_mtime = _get_native_source_mtime(root_dir)
    for candidate in _get_local_native_artifact_candidates(root_dir):
        if candidate.stat().st_mtime < source_mtime:
            return candidate

    return None


def _load_native_module_from_artifact_path(candidate: Path) -> ModuleType | None:
    """
    Load a native module from a specific local artifact path.
    """
    spec = importlib.util.spec_from_file_location("volumizer_native", candidate)
    if spec is None or spec.loader is None:
        return None

    module = importlib.util.module_from_spec(spec)
    sys.modules["volumizer_native"] = module
    spec.loader.exec_module(module)
    return module


def _load_native_module_from_local_artifact() -> ModuleType | None:
    """
    Load a locally-built native module artifact (debug/release target) when present.
    """
    root_dir = Path(__file__).resolve().parents[1]
    source_mtime = _get_native_source_mtime(root_dir)

    for candidate in _get_local_native_artifact_candidates(root_dir):
        if candidate.stat().st_mtime < source_mtime:
            continue

        module = _load_native_module_from_artifact_path(candidate)
        if module is not None:
            return module

    return None


def get_requested_backend() -> str:
    """
    Return requested backend mode from environment.
    """
    backend = os.getenv(BACKEND_ENV, DEFAULT_BACKEND).strip().lower()
    if backend not in VALID_BACKENDS:
        return BACKEND_AUTO
    return backend


@lru_cache(maxsize=3)
def _load_native_module(mode: str) -> ModuleType | None:
    """
    Attempt to import native backend based on mode.
    """
    if mode == BACKEND_PYTHON:
        return None

    try:
        return importlib.import_module("volumizer_native")
    except ImportError as error:
        local_module = _load_native_module_from_local_artifact()
        if local_module is not None:
            return local_module

        stale_artifact = _find_stale_local_native_artifact()
        if stale_artifact is not None and mode == BACKEND_NATIVE:
            raise RuntimeError(
                "Native backend requested but local native artifact is stale. "
                f"Rebuild `{stale_artifact.name}` or set VOLUMIZER_BACKEND=python."
            ) from error

        if mode == BACKEND_NATIVE:
            raise RuntimeError(
                "Native backend requested but `volumizer_native` is not importable. "
                "Install/build native module or set VOLUMIZER_BACKEND=python."
            ) from error
        return None


def clear_backend_cache() -> None:
    """
    Clear cached backend imports. Mainly useful for tests.
    """
    _load_native_module.cache_clear()


def get_native_module() -> ModuleType | None:
    """
    Return native module if available under selected mode.
    """
    return _load_native_module(get_requested_backend())


def get_native_module_for_mode(mode: str) -> ModuleType | None:
    """
    Return native module for an explicit backend mode.
    """
    normalized_mode = mode.strip().lower()
    if normalized_mode not in VALID_BACKENDS:
        normalized_mode = BACKEND_AUTO
    return _load_native_module(normalized_mode)


def using_native() -> bool:
    """
    Return True when native backend is selected and importable.
    """
    return get_native_module() is not None


def active_backend() -> str:
    """
    Return currently active backend.
    """
    if using_native():
        return BACKEND_NATIVE
    return BACKEND_PYTHON
