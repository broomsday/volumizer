from pathlib import Path

import numpy as np
import pytest
from biotite.structure.io import pdbx

from volumizer import pdb
from volumizer.paths import TEST_DIR
from volumizer.types import VoxelGroup


TEST_PDB_DIR = TEST_DIR / "pdbs"
TEST_IDENTITY_CIF = TEST_PDB_DIR / "4jpn.cif"
TEST_FALLBACK_CIF = TEST_PDB_DIR / "4jpp.cif"


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


def test_load_structure_rejects_invalid_assembly_policy():
    with pytest.raises(ValueError, match="Unsupported assembly policy"):
        pdb.load_structure(Path("dummy.cif"), assembly_policy="invalid")


def test_load_structure_asymmetric_uses_get_structure_for_cif(monkeypatch, tmp_path: Path):
    sentinel = object()
    dummy_file = object()
    called = {"structure": 0, "assembly": 0}

    monkeypatch.setattr(pdbx.PDBxFile, "read", lambda _: dummy_file)

    def _get_structure(file, model=1):
        called["structure"] += 1
        assert file is dummy_file
        assert model == 1
        return sentinel

    def _get_assembly(*args, **kwargs):
        called["assembly"] += 1
        raise AssertionError("get_assembly should not be used for asymmetric policy")

    monkeypatch.setattr(pdbx, "get_structure", _get_structure)
    monkeypatch.setattr(pdbx, "get_assembly", _get_assembly)

    stage_timings: dict[str, float] = {}
    result = pdb.load_structure(
        tmp_path / "dummy.cif",
        assembly_policy="asymmetric",
        stage_timings=stage_timings,
    )

    assert result is sentinel
    assert called == {"structure": 1, "assembly": 0}
    assert stage_timings["load_structure_parse_decode"] >= 0.0
    assert "load_structure_assembly_expand" not in stage_timings
    assert "load_structure_fallback" not in stage_timings


def test_can_use_identity_assembly_shortcut_detects_expected_examples():
    identity_file = pdbx.PDBxFile.read(TEST_IDENTITY_CIF)
    fallback_file = pdbx.PDBxFile.read(TEST_FALLBACK_CIF)

    assert pdb._can_use_identity_assembly_shortcut_cif(identity_file)
    assert not pdb._can_use_identity_assembly_shortcut_cif(fallback_file)


def test_load_structure_auto_shortcuts_identity_assembly(monkeypatch):
    sentinel = object()
    called = {"structure": 0, "assembly": 0}

    def _get_structure(file, model=1):
        called["structure"] += 1
        assert model == 1
        return sentinel

    def _get_assembly(*args, **kwargs):
        called["assembly"] += 1
        raise AssertionError("identity shortcut should skip get_assembly in auto mode")

    monkeypatch.setattr(pdbx, "get_structure", _get_structure)
    monkeypatch.setattr(pdbx, "get_assembly", _get_assembly)

    stage_timings: dict[str, float] = {}
    result = pdb.load_structure(
        TEST_IDENTITY_CIF,
        assembly_policy="auto",
        stage_timings=stage_timings,
    )

    assert result is sentinel
    assert called == {"structure": 1, "assembly": 0}
    assert stage_timings["load_structure_parse_decode"] >= 0.0
    assert "load_structure_assembly_expand" not in stage_timings
    assert "load_structure_fallback" not in stage_timings


def test_load_structure_records_fallback_stage_on_missing_assembly_metadata():
    stage_timings: dict[str, float] = {}

    structure = pdb.load_structure(
        TEST_FALLBACK_CIF,
        assembly_policy="biological",
        stage_timings=stage_timings,
    )

    assert len(structure) > 0
    assert stage_timings["load_structure_parse_decode"] > 0.0
    assert stage_timings["load_structure_assembly_expand"] > 0.0
    assert stage_timings["load_structure_fallback"] > 0.0


def test_load_structure_records_parse_and_assembly_stages_for_biological_cif():
    stage_timings: dict[str, float] = {}

    structure = pdb.load_structure(
        TEST_IDENTITY_CIF,
        assembly_policy="biological",
        stage_timings=stage_timings,
    )

    assert len(structure) > 0
    assert stage_timings["load_structure_parse_decode"] > 0.0
    assert stage_timings["load_structure_assembly_expand"] > 0.0
    assert "load_structure_fallback" not in stage_timings
