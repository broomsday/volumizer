import numpy as np

from volumizer import pdb
from volumizer.types import VoxelGroup


class DummyVoxelGrid:
    def __init__(self, voxel_centers: np.ndarray):
        self.voxel_centers = voxel_centers


def _make_voxel_group(indices: set[int], surface_indices: set[int]) -> VoxelGroup:
    empty = np.array([], dtype=np.int64)
    return VoxelGroup(
        voxels=(empty, empty, empty),
        indices=indices,
        num_voxels=len(indices),
        surface_indices=surface_indices,
    )


def test_volume_to_structure_assigns_surface_and_metadata():
    voxel_centers = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ],
        dtype=np.float64,
    )

    structure = pdb.volume_to_structure(
        voxel_type="POR",
        voxel_group_index=7,
        voxel_indices={0, 2},
        surface_indices={2},
        voxel_grid_centers=voxel_centers,
    )

    assert len(structure) == 2

    atom_ids = np.asarray(structure.atom_id, dtype=np.int64)
    b_factors = np.asarray(structure.b_factor, dtype=np.float64)
    by_id_bfactor = {
        int(atom_id): float(beta)
        for atom_id, beta in zip(atom_ids, b_factors)
    }

    assert by_id_bfactor[0] == 0.0
    assert by_id_bfactor[2] == 50.0
    assert np.allclose(structure.coord, voxel_centers[atom_ids])

    assert set(np.asarray(structure.atom_name).tolist()) == {"O"}
    assert set(np.asarray(structure.res_name).tolist()) == {"POR"}
    assert set(np.asarray(structure.res_id, dtype=np.int64).tolist()) == {7}
    assert set(np.asarray(structure.chain_id).tolist()) == {"B"}
    assert set(np.asarray(structure.element).tolist()) == {"O"}


def test_volumes_to_structure_combines_voxel_groups():
    voxel_centers = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ],
        dtype=np.float64,
    )

    hubs = {1: _make_voxel_group(indices={0, 1}, surface_indices={1})}
    pores = {2: _make_voxel_group(indices={2}, surface_indices={2})}

    structure = pdb.volumes_to_structure(
        voxel_grid=DummyVoxelGrid(voxel_centers),
        hubs=hubs,
        pores=pores,
        pockets={},
        cavities={},
        occluded={},
    )

    assert len(structure) == 3

    expected = {
        0: ("HUB", "A", "F", "F", 1, 0.0),
        1: ("HUB", "A", "F", "F", 1, 50.0),
        2: ("POR", "B", "O", "O", 2, 50.0),
    }

    for index in range(len(structure)):
        atom_id = int(structure.atom_id[index])
        assert atom_id in expected

        res_name, chain_id, atom_name, element, res_id, b_factor = expected[atom_id]
        assert np.allclose(structure.coord[index], voxel_centers[atom_id])
        assert str(structure.res_name[index]) == res_name
        assert str(structure.chain_id[index]) == chain_id
        assert str(structure.atom_name[index]) == atom_name
        assert str(structure.element[index]) == element
        assert int(structure.res_id[index]) == res_id
        assert float(structure.b_factor[index]) == b_factor
