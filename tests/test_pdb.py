from pathlib import Path
import itertools

import biotite.structure as bts
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


def _make_test_protein_structure() -> bts.AtomArray:
    structure = bts.AtomArray(1)
    structure.coord = np.array([[10.0, 0.0, 0.0]], dtype=np.float64)
    structure.atom_name = np.array(["CA"], dtype=object)
    structure.res_name = np.array(["ALA"], dtype=object)
    structure.res_id = np.array([1], dtype=np.int64)
    structure.chain_id = np.array(["A"], dtype=object)
    structure.element = np.array(["C"], dtype=object)
    return structure


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
    assert {"atom_id", "b_factor"}.issubset(structure.get_annotation_categories())

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
    assert {"atom_id", "b_factor"}.issubset(structure.get_annotation_categories())

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


def test_volumes_to_structure_uses_display_type_overrides_for_labels():
    voxel_centers = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ],
        dtype=np.float64,
    )

    pockets = {5: _make_voxel_group(indices={0, 1}, surface_indices={1})}
    cavities = {9: _make_voxel_group(indices={2}, surface_indices={2})}

    structure = pdb.volumes_to_structure(
        voxel_grid=DummyVoxelGrid(voxel_centers),
        hubs={},
        pores={},
        pockets=pockets,
        cavities=cavities,
        occluded={},
        display_type_overrides={("POK", 5): "CAV"},
    )

    assert len(structure) == 3
    by_atom_id = {
        int(structure.atom_id[index]): (
            str(structure.res_name[index]),
            str(structure.chain_id[index]),
            str(structure.atom_name[index]),
            str(structure.element[index]),
            int(structure.res_id[index]),
            float(structure.b_factor[index]),
        )
        for index in range(len(structure))
    }

    assert by_atom_id == {
        0: ("CAV", "D", "S", "S", 5, 0.0),
        1: ("CAV", "D", "S", "S", 5, 50.0),
        2: ("CAV", "D", "S", "S", 9, 50.0),
    }
    assert "POK" not in set(np.asarray(structure.res_name).tolist())


def test_save_structure_cif_preserves_b_factor_after_concat(tmp_path: Path):
    protein_structure = _make_test_protein_structure()
    pdb.ensure_b_factor_annotation(protein_structure)

    voxel_centers = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ],
        dtype=np.float64,
    )
    volume_structure = pdb.volume_to_structure(
        voxel_type="POR",
        voxel_group_index=7,
        voxel_indices=np.array([0, 2], dtype=np.int64),
        surface_indices=np.array([2], dtype=np.int64),
        voxel_grid_centers=voxel_centers,
    )

    combined_structure = protein_structure + volume_structure
    out_path = tmp_path / "combined.annotated.cif"
    pdb.save_structure(combined_structure, out_path)

    pdbx_file = pdbx.PDBxFile.read(out_path)
    atom_site = pdbx_file.get_category("atom_site", expect_looped=True)

    assert "B_iso_or_equiv" in atom_site
    assert atom_site["B_iso_or_equiv"].tolist() == ["0.00", "0.00", "50.00"]


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


def test_compute_sse_fractions_returns_valid_fractions():
    structure = pdb.load_structure(TEST_PDB_DIR / "cavity.pdb")
    cleaned = pdb.clean_structure(structure)
    result = pdb.compute_sse_fractions(cleaned)

    assert set(result.keys()) == {"frac_alpha", "frac_beta", "frac_coil"}
    for key in ("frac_alpha", "frac_beta", "frac_coil"):
        value = result[key]
        assert value is None or (0.0 <= value <= 1.0), f"{key}={value}"

    non_none = [v for v in result.values() if v is not None]
    if non_none:
        assert abs(sum(non_none) - 1.0) < 0.01


class TestDeduplicateAssemblyChainIds:
    def _make_atom_array(self, chain_ids: list[str]) -> bts.AtomArray:
        """Build a minimal AtomArray with the given chain_id sequence."""
        n = len(chain_ids)
        atoms = bts.AtomArray(n)
        atoms.chain_id = np.array(chain_ids, dtype="U4")
        atoms.coord = np.zeros((n, 3), dtype=np.float32)
        atoms.atom_name = np.full(n, "CA", dtype="U4")
        atoms.res_name = np.full(n, "ALA", dtype="U4")
        atoms.res_id = np.arange(n, dtype=np.int32)
        atoms.element = np.full(n, "C", dtype="U2")
        return atoms

    def test_single_copy_is_noop(self):
        atoms = self._make_atom_array(["A", "A", "B", "B", "C"])
        result = pdb._deduplicate_assembly_chain_ids(atoms)
        assert list(result.chain_id) == ["A", "A", "B", "B", "C"]

    def test_two_copies_get_unique_ids(self):
        # Simulate 2 copies of A,B
        atoms = self._make_atom_array(["A", "A", "B", "B", "A", "A", "B", "B"])
        result = pdb._deduplicate_assembly_chain_ids(atoms)
        chain_ids = list(result.chain_id)
        # Copy 0 keeps A, B; copy 1 gets new unique IDs
        assert chain_ids[:4] == ["A", "A", "B", "B"]
        assert chain_ids[4] != "A" and chain_ids[4] != "B"
        assert chain_ids[6] != "A" and chain_ids[6] != "B"
        assert chain_ids[4] == chain_ids[5]  # same new ID within block
        assert chain_ids[6] == chain_ids[7]
        assert chain_ids[4] != chain_ids[6]  # different for A-copy vs B-copy
        # Total unique IDs = 4
        assert len(set(chain_ids)) == 4

    def test_three_copies_trimer(self):
        # Simulate 3 copies of A,B,C,D (like 2ZBT pattern)
        pattern = ["A"] * 3 + ["B"] * 3 + ["C"] * 3 + ["D"] * 3
        atoms = self._make_atom_array(pattern * 3)
        result = pdb._deduplicate_assembly_chain_ids(atoms)
        chain_ids = list(result.chain_id)
        # Should have 12 unique chain IDs
        assert len(set(chain_ids)) == 12
        # Copy 0 retains original IDs
        assert chain_ids[0] == "A"
        assert chain_ids[3] == "B"
        assert chain_ids[6] == "C"
        assert chain_ids[9] == "D"

    def test_empty_structure(self):
        atoms = self._make_atom_array([])
        result = pdb._deduplicate_assembly_chain_ids(atoms)
        assert len(result) == 0

    def test_original_not_mutated(self):
        atoms = self._make_atom_array(["A", "B", "A", "B"])
        original_ids = list(atoms.chain_id)
        pdb._deduplicate_assembly_chain_ids(atoms)
        assert list(atoms.chain_id) == original_ids

    def _make_single_chain_copies(self, atoms_per_copy: int, num_copies: int) -> bts.AtomArray:
        """Build an AtomArray simulating multiple copies of a single chain."""
        total = atoms_per_copy * num_copies
        atoms = bts.AtomArray(total)
        atoms.chain_id = np.full(total, "A", dtype="U4")
        atoms.coord = np.zeros((total, 3), dtype=np.float32)
        atoms.atom_name = np.full(total, "CA", dtype="U4")
        atoms.res_name = np.full(total, "ALA", dtype="U4")
        atoms.element = np.full(total, "C", dtype="U2")
        # Each copy restarts res_id from 1
        res_ids = np.tile(np.arange(1, atoms_per_copy + 1), num_copies)
        atoms.res_id = res_ids.astype(np.int32)
        return atoms

    def test_single_chain_copies_get_unique_ids(self):
        # Simulate 6 copies of a single chain (like 3R1K pattern)
        atoms = self._make_single_chain_copies(atoms_per_copy=100, num_copies=6)
        result = pdb._deduplicate_assembly_chain_ids(atoms)
        chain_ids = list(result.chain_id)
        assert len(set(chain_ids)) == 6
        # Copy 0 keeps original chain ID
        assert chain_ids[0] == "A"
        # Other copies get new IDs
        assert len(set(chain_ids[100:])) == 5

    def test_single_chain_copies_8x(self):
        # Simulate 8 copies of a single chain (like 6L7D pattern)
        atoms = self._make_single_chain_copies(atoms_per_copy=50, num_copies=8)
        result = pdb._deduplicate_assembly_chain_ids(atoms)
        assert len(set(result.chain_id)) == 8

    def test_single_chain_no_copies_is_noop(self):
        atoms = self._make_single_chain_copies(atoms_per_copy=100, num_copies=1)
        result = pdb._deduplicate_assembly_chain_ids(atoms)
        assert list(result.chain_id) == ["A"] * 100

    def test_many_copies_do_not_exhaust_chain_id_pool(self):
        pattern = []
        for chain_id in [
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "J",
            "K",
            "L",
            "M",
            "N",
            "O",
            "P",
        ]:
            pattern.extend([chain_id] * 2)
        atoms = self._make_atom_array(pattern * 120)
        result = pdb._deduplicate_assembly_chain_ids(atoms)
        assert len(set(result.chain_id)) == 16 * 120

    def test_chain_id_pool_order_is_stable(self):
        observed = list(itertools.islice(pdb._iter_chain_id_pool(), 66))
        assert observed[:5] == ["A", "B", "C", "D", "E"]
        assert observed[60:66] == ["8", "9", "AA", "AB", "AC", "AD"]


RCSB_2ZBT = Path("data/runs/rcsb70/downloads/2ZBT.cif")
RCSB_6L7D = Path("data/runs/rcsb70/downloads/6L7D.cif")
RCSB_3R1K = Path("data/runs/rcsb70/downloads/3R1K.cif")


@pytest.mark.skipif(not RCSB_2ZBT.exists(), reason="2ZBT test data not available")
class TestBiologicalAssembly2ZBT:
    def test_biological_assembly_has_12_unique_chains(self):
        structure = pdb.load_structure(RCSB_2ZBT, assembly_policy="biological")
        cleaned = pdb.clean_structure(structure)
        chain_ids = set(cleaned.chain_id)
        assert len(chain_ids) == 12

    def test_biological_assembly_has_3x_asymmetric_atoms(self):
        bio = pdb.load_structure(RCSB_2ZBT, assembly_policy="biological")
        asym = pdb.load_structure(RCSB_2ZBT, assembly_policy="asymmetric")
        # Assembly should be exactly 3x the asymmetric unit
        assert len(bio) == 3 * len(asym)

    def test_asymmetric_unit_has_4_chains(self):
        structure = pdb.load_structure(RCSB_2ZBT, assembly_policy="asymmetric")
        cleaned = pdb.clean_structure(structure)
        chain_ids = set(cleaned.chain_id)
        assert len(chain_ids) == 4


@pytest.mark.skipif(not RCSB_6L7D.exists(), reason="6L7D test data not available")
class TestBiologicalAssembly6L7D:
    """6L7D has 1 protein chain in the asymmetric unit and 8 symmetry copies."""

    def test_biological_assembly_has_8_unique_chains(self):
        structure = pdb.load_structure(RCSB_6L7D, assembly_policy="biological")
        cleaned = pdb.clean_structure(structure)
        chain_ids = set(cleaned.chain_id)
        assert len(chain_ids) == 8

    def test_biological_assembly_has_8x_asymmetric_atoms(self):
        bio = pdb.load_structure(RCSB_6L7D, assembly_policy="biological")
        asym = pdb.load_structure(RCSB_6L7D, assembly_policy="asymmetric")
        bio_protein = bio[bts.filter_amino_acids(bio)]
        asym_protein = asym[bts.filter_amino_acids(asym)]
        assert len(bio_protein) == 8 * len(asym_protein)


@pytest.mark.skipif(not RCSB_3R1K.exists(), reason="3R1K test data not available")
class TestBiologicalAssembly3R1K:
    """3R1K has 1 protein chain in the asymmetric unit and 6 symmetry copies."""

    def test_biological_assembly_has_6_unique_chains(self):
        structure = pdb.load_structure(RCSB_3R1K, assembly_policy="biological")
        cleaned = pdb.clean_structure(structure)
        chain_ids = set(cleaned.chain_id)
        assert len(chain_ids) == 6

    def test_biological_assembly_has_6x_asymmetric_atoms(self):
        bio = pdb.load_structure(RCSB_3R1K, assembly_policy="biological")
        asym = pdb.load_structure(RCSB_3R1K, assembly_policy="asymmetric")
        bio_protein = bio[bts.filter_amino_acids(bio)]
        asym_protein = asym[bts.filter_amino_acids(asym)]
        assert len(bio_protein) == 6 * len(asym_protein)
