"""
Helpers for tuning Mol* representations used by the gallery UI/renderers.
"""

from __future__ import annotations

import math


DEFAULT_VOLUME_REFERENCE_RESOLUTION = 3.0
DEFAULT_VOLUME_SURFACE_SMOOTHNESS = 1.5
VOLUME_SURFACE_STYLE_VERSION = "resolution-linked-v1"


def compute_volume_surface_radius_offset(
    voxel_resolution: float | None,
) -> float:
    """
    Increase gaussian-surface radius as voxel spacing gets coarser.

    The annotated volume pseudo-atoms sit at voxel centers. Around the default
    3.0 A spacing, Mol*'s physical radii already produce a reasonable surface,
    so no offset is needed. At coarser spacings the centers drift farther apart,
    so expand the Gaussian surface by half of the excess spacing.
    """
    if voxel_resolution is None:
        return 0.0

    resolution = float(voxel_resolution)
    if not math.isfinite(resolution):
        return 0.0

    if resolution <= DEFAULT_VOLUME_REFERENCE_RESOLUTION:
        return 0.0

    return round((resolution - DEFAULT_VOLUME_REFERENCE_RESOLUTION) / 2.0, 3)


def build_volume_surface_style(
    voxel_resolution: float | None,
) -> dict[str, float | str]:
    """
    Return Mol* gaussian-surface tuning values for volumizer pseudo-atoms.
    """
    return {
        "quality": "custom",
        "radius_offset": compute_volume_surface_radius_offset(voxel_resolution),
        "smoothness": DEFAULT_VOLUME_SURFACE_SMOOTHNESS,
    }
