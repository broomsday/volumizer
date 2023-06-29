"""
Definition of custom types.
"""


from typing import NamedTuple, Optional

import numpy as np


class ComponentData(NamedTuple):
    component_id: str
    component_type: str


class VoxelGroup(NamedTuple):
    voxels: tuple[np.ndarray, ...]
    indices: set[int]
    num_voxels: int
    surface_indices: set[int] = set()
    voxel_type: Optional[str] = None
    volume: float = 0.0
    center: Optional[np.ndarray] = None
    axial_lengths: list[float] = [0.0, 0.0, 0.0]


class Annotation(NamedTuple):
    total_hub_volume: float
    total_pore_volume: float
    total_cavity_volume: float
    total_pocket_volume: float
    largest_hub_volume: float
    largest_pore_volume: float
    largest_cavity_volume: float
    largest_pocket_volume: float
    num_hubs: int
    num_pores: int
    num_cavities: int
    num_pockets: int
    hub_volumes: dict[int, Optional[float]]
    pore_volumes: dict[int, Optional[float]]
    cavity_volumes: dict[int, Optional[float]]
    pocket_volumes: dict[int, Optional[float]]
    hub_dimensions: dict[int, list[float]]
    pore_dimensions: dict[int, list[float]]
    cavity_dimensions: dict[int, list[float]]
    pocket_dimensions: dict[int, list[float]]
