"""
Definition of custom types.
"""


from typing import NamedTuple, Optional

import numpy as np


class ComponentData(NamedTuple):
    component_id: str
    component_type: str


class _VoxelGroupBase(NamedTuple):
    voxels: tuple[np.ndarray, ...]
    indices: set[int]
    num_voxels: int
    surface_indices: set[int] | np.ndarray | None = None
    voxel_type: Optional[str] = None
    volume: float = 0.0
    center: Optional[np.ndarray] = None
    axial_lengths: list[float] | None = None
    cross_section_circularity: Optional[float] = None
    cross_section_uniformity: Optional[float] = None


class VoxelGroup(_VoxelGroupBase):
    __slots__ = ()

    def __new__(
        cls,
        voxels: tuple[np.ndarray, ...],
        indices: set[int],
        num_voxels: int,
        surface_indices: set[int] | np.ndarray | None = None,
        voxel_type: Optional[str] = None,
        volume: float = 0.0,
        center: Optional[np.ndarray] = None,
        axial_lengths: list[float] | None = None,
        cross_section_circularity: Optional[float] = None,
        cross_section_uniformity: Optional[float] = None,
    ):
        if surface_indices is None:
            surface_indices = set()
        if axial_lengths is None:
            axial_lengths = [0.0, 0.0, 0.0]

        return super(VoxelGroup, cls).__new__(
            cls,
            voxels,
            indices,
            num_voxels,
            surface_indices,
            voxel_type,
            volume,
            center,
            axial_lengths,
            cross_section_circularity,
            cross_section_uniformity,
        )


class _AnnotationBase(NamedTuple):
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
    hub_cross_section_circularity: dict[int, Optional[float]] | None = None
    pore_cross_section_circularity: dict[int, Optional[float]] | None = None
    cavity_cross_section_circularity: dict[int, Optional[float]] | None = None
    pocket_cross_section_circularity: dict[int, Optional[float]] | None = None
    hub_cross_section_uniformity: dict[int, Optional[float]] | None = None
    pore_cross_section_uniformity: dict[int, Optional[float]] | None = None
    cavity_cross_section_uniformity: dict[int, Optional[float]] | None = None
    pocket_cross_section_uniformity: dict[int, Optional[float]] | None = None


class Annotation(_AnnotationBase):
    __slots__ = ()

    def __new__(
        cls,
        total_hub_volume: float,
        total_pore_volume: float,
        total_cavity_volume: float,
        total_pocket_volume: float,
        largest_hub_volume: float,
        largest_pore_volume: float,
        largest_cavity_volume: float,
        largest_pocket_volume: float,
        num_hubs: int,
        num_pores: int,
        num_cavities: int,
        num_pockets: int,
        hub_volumes: dict[int, Optional[float]],
        pore_volumes: dict[int, Optional[float]],
        cavity_volumes: dict[int, Optional[float]],
        pocket_volumes: dict[int, Optional[float]],
        hub_dimensions: dict[int, list[float]],
        pore_dimensions: dict[int, list[float]],
        cavity_dimensions: dict[int, list[float]],
        pocket_dimensions: dict[int, list[float]],
        hub_cross_section_circularity: dict[int, Optional[float]] | None = None,
        pore_cross_section_circularity: dict[int, Optional[float]] | None = None,
        cavity_cross_section_circularity: dict[int, Optional[float]] | None = None,
        pocket_cross_section_circularity: dict[int, Optional[float]] | None = None,
        hub_cross_section_uniformity: dict[int, Optional[float]] | None = None,
        pore_cross_section_uniformity: dict[int, Optional[float]] | None = None,
        cavity_cross_section_uniformity: dict[int, Optional[float]] | None = None,
        pocket_cross_section_uniformity: dict[int, Optional[float]] | None = None,
    ):
        return super(Annotation, cls).__new__(
            cls,
            total_hub_volume,
            total_pore_volume,
            total_cavity_volume,
            total_pocket_volume,
            largest_hub_volume,
            largest_pore_volume,
            largest_cavity_volume,
            largest_pocket_volume,
            num_hubs,
            num_pores,
            num_cavities,
            num_pockets,
            hub_volumes,
            pore_volumes,
            cavity_volumes,
            pocket_volumes,
            hub_dimensions,
            pore_dimensions,
            cavity_dimensions,
            pocket_dimensions,
            {} if hub_cross_section_circularity is None else hub_cross_section_circularity,
            {} if pore_cross_section_circularity is None else pore_cross_section_circularity,
            {} if cavity_cross_section_circularity is None else cavity_cross_section_circularity,
            {} if pocket_cross_section_circularity is None else pocket_cross_section_circularity,
            {} if hub_cross_section_uniformity is None else hub_cross_section_uniformity,
            {} if pore_cross_section_uniformity is None else pore_cross_section_uniformity,
            {} if cavity_cross_section_uniformity is None else cavity_cross_section_uniformity,
            {} if pocket_cross_section_uniformity is None else pocket_cross_section_uniformity,
        )
