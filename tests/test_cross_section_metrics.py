import numpy as np
import pytest

from volumizer.voxel import (
    _align_to_principal_axes,
    _compute_cross_section_metrics_from_indices,
)


class FakeVoxelGrid:
    """Minimal stand-in for pyntcloud VoxelGrid with voxel_centers and x_y_z."""

    def __init__(self, centers: np.ndarray, voxel_size: float = 1.0):
        self.voxel_centers = centers
        self.x_y_z = np.array([voxel_size, voxel_size, voxel_size])


def _cylinder_voxels(radius: float, length: float, voxel_size: float = 1.0) -> np.ndarray:
    """Generate voxel centers approximating a cylinder along the X axis."""
    centers = []
    for x in np.arange(-length / 2, length / 2 + voxel_size / 2, voxel_size):
        for y in np.arange(-radius, radius + voxel_size / 2, voxel_size):
            for z in np.arange(-radius, radius + voxel_size / 2, voxel_size):
                if y**2 + z**2 <= radius**2:
                    centers.append([x, y, z])
    return np.array(centers, dtype=np.float64)


def _elliptical_tube_voxels(
    radius_y: float, radius_z: float, length: float, voxel_size: float = 1.0
) -> np.ndarray:
    """Generate voxel centers approximating an elliptical tube along the X axis."""
    centers = []
    for x in np.arange(-length / 2, length / 2 + voxel_size / 2, voxel_size):
        for y in np.arange(-radius_y, radius_y + voxel_size / 2, voxel_size):
            for z in np.arange(-radius_z, radius_z + voxel_size / 2, voxel_size):
                if (y / radius_y) ** 2 + (z / radius_z) ** 2 <= 1.0:
                    centers.append([x, y, z])
    return np.array(centers, dtype=np.float64)


def test_align_to_principal_axes_returns_none_for_too_few_voxels():
    grid = FakeVoxelGrid(np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    indices = np.array([0, 1], dtype=np.int64)
    assert _align_to_principal_axes(indices, grid) is None


def test_align_to_principal_axes_basic():
    centers = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.5, 0.0]],
        dtype=np.float64,
    )
    grid = FakeVoxelGrid(centers)
    indices = np.array([0, 1, 2, 3], dtype=np.int64)
    aligned = _align_to_principal_axes(indices, grid)
    assert aligned is not None
    assert aligned.shape == (4, 3)


def test_cross_section_metrics_too_few_voxels():
    grid = FakeVoxelGrid(np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    indices = np.array([0, 1], dtype=np.int64)
    circ, unif = _compute_cross_section_metrics_from_indices(indices, grid)
    assert circ is None
    assert unif is None


def test_cross_section_metrics_degenerate_flat():
    centers = np.array(
        [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    grid = FakeVoxelGrid(centers)
    indices = np.array([0, 1, 2], dtype=np.int64)
    circ, unif = _compute_cross_section_metrics_from_indices(indices, grid)
    # All on same plane along primary axis => uniformity should be high (or 1.0)
    assert circ is None or 0.0 <= circ <= 1.0
    assert unif is None or 0.0 <= unif <= 1.0


def test_cylinder_has_high_circularity_and_uniformity():
    centers = _cylinder_voxels(radius=5.0, length=20.0, voxel_size=1.0)
    grid = FakeVoxelGrid(centers, voxel_size=1.0)
    indices = np.arange(len(centers), dtype=np.int64)

    circ, unif = _compute_cross_section_metrics_from_indices(indices, grid)

    assert circ is not None
    assert unif is not None
    # Cylinder should have high circularity (near 1.0)
    assert circ > 0.8, f"Expected circularity > 0.8 for cylinder, got {circ}"
    # Cylinder should have high uniformity (near 1.0; edge slices have fewer voxels)
    assert unif > 0.75, f"Expected uniformity > 0.75 for cylinder, got {unif}"


def test_elliptical_tube_has_lower_circularity():
    centers = _elliptical_tube_voxels(radius_y=8.0, radius_z=2.0, length=20.0, voxel_size=1.0)
    grid = FakeVoxelGrid(centers, voxel_size=1.0)
    indices = np.arange(len(centers), dtype=np.int64)

    circ, unif = _compute_cross_section_metrics_from_indices(indices, grid)

    assert circ is not None
    assert unif is not None
    # Elliptical cross-section should have lower circularity
    assert circ < 0.5, f"Expected circularity < 0.5 for 4:1 ellipse, got {circ}"
    # But uniform cross-section along length => high uniformity (edge effects lower it slightly)
    assert unif > 0.75, f"Expected uniformity > 0.75 for uniform tube, got {unif}"


def test_cone_has_low_uniformity():
    """A cone-like shape (radius varies linearly along axis) should have low uniformity."""
    centers = []
    for x in np.arange(0, 20, 1.0):
        radius = 1.0 + x * 0.5  # radius grows from 1 to 10.5
        for y in np.arange(-radius, radius + 0.5, 1.0):
            for z in np.arange(-radius, radius + 0.5, 1.0):
                if y**2 + z**2 <= radius**2:
                    centers.append([x, y, z])
    centers = np.array(centers, dtype=np.float64)
    grid = FakeVoxelGrid(centers, voxel_size=1.0)
    indices = np.arange(len(centers), dtype=np.int64)

    circ, unif = _compute_cross_section_metrics_from_indices(indices, grid)

    assert circ is not None
    assert unif is not None
    # Cone cross-sections are circular (small slices near tip lower the average)
    assert circ > 0.4, f"Expected circularity > 0.4 for cone, got {circ}"
    # But varying radius => lower uniformity
    assert unif < 0.7, f"Expected uniformity < 0.7 for cone, got {unif}"
