"""
Utilities for selecting and loading the optional native Rust backend.
"""

from __future__ import annotations

from functools import lru_cache
import importlib
import os
from types import ModuleType


BACKEND_ENV = "VOLUMIZER_BACKEND"
BACKEND_AUTO = "auto"
BACKEND_NATIVE = "native"
BACKEND_PYTHON = "python"
VALID_BACKENDS = {BACKEND_AUTO, BACKEND_NATIVE, BACKEND_PYTHON}


def get_requested_backend() -> str:
    """
    Return requested backend mode from environment.
    """
    backend = os.getenv(BACKEND_ENV, BACKEND_AUTO).strip().lower()
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
