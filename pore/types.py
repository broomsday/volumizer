"""
Definition of custom types.
"""


from typing import NamedTuple, Optional

import numpy as np


class ComponentData(NamedTuple):
    component_id: str
    component_type: str


class VoxelGroup(NamedTuple):
    voxels: tuple[np.ndarray, np.ndarray, np.ndarray]
    indices: set[int]
    num_voxels: int
    voxel_type: Optional[str] = None
    volume: Optional[float] = None


class Annotation(NamedTuple):
    total_pore_volume: float
    total_cavity_volume: float
    largest_pore_volume: float
    largest_cavity_volume: float
    pore_volumes: dict[int, float]
    cavity_volumes: dict[int, float]
