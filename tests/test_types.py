import numpy as np

from volumizer.types import Annotation, VoxelGroup


def test_voxel_group_defaults_are_not_shared():
    empty_voxels = (
        np.array([], dtype=int),
        np.array([], dtype=int),
        np.array([], dtype=int),
    )
    first = VoxelGroup(voxels=empty_voxels, indices=set(), num_voxels=0)
    second = VoxelGroup(voxels=empty_voxels, indices=set(), num_voxels=0)

    first.surface_indices.add(1)
    first.axial_lengths.append(2.0)

    assert second.surface_indices == set()
    assert second.axial_lengths == [0.0, 0.0, 0.0]


def _empty_annotation() -> Annotation:
    return Annotation(
        total_hub_volume=0.0,
        total_pore_volume=0.0,
        total_cavity_volume=0.0,
        total_pocket_volume=0.0,
        largest_hub_volume=0.0,
        largest_pore_volume=0.0,
        largest_cavity_volume=0.0,
        largest_pocket_volume=0.0,
        num_hubs=0,
        num_pores=0,
        num_cavities=0,
        num_pockets=0,
        hub_volumes={},
        pore_volumes={},
        cavity_volumes={},
        pocket_volumes={},
        hub_dimensions={},
        pore_dimensions={},
        cavity_dimensions={},
        pocket_dimensions={},
    )


def test_annotation_metric_defaults_are_not_shared():
    first = _empty_annotation()
    second = _empty_annotation()

    first.hub_cross_section_circularity[1] = 0.5
    first.pocket_cross_section_uniformity[2] = 0.25

    assert second.hub_cross_section_circularity == {}
    assert second.pocket_cross_section_uniformity == {}
