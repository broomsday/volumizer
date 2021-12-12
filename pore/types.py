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
    voxel_type: Optional[str] = None
    volume: Optional[float] = None


class Annotation(NamedTuple):
    total_pore_volume: float
    total_cavity_volume: float
    total_pocket_volume: float
    largest_pore_volume: float
    largest_cavity_volume: float
    largest_pocket_volume: float
    num_pores: int
    num_cavities: int
    num_pockets: int
    pore_volumes: dict[int, Optional[float]]
    cavity_volumes: dict[int, Optional[float]]
    pocket_volumes: dict[int, Optional[float]]
