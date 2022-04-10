"""
Functions for parsing, cleaning, and modifying PDBs.
"""

from pathlib import Path
import tempfile
import itertools
from typing import Optional
import warnings
import subprocess

from Bio.PDB import PDBParser, Select, Structure, DSSP
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.PDBIO import PDBIO
from Bio import pairwise2
from pyntcloud.structures.voxelgrid import VoxelGrid
import pandas as pd
import numpy as np

from pore.constants import (
    VOXEL_ATOM_NAMES,
    VOXEL_TYPE_CHAIN_MAP,
    VOXEL_TYPE_ATOM_MAP,
    RESIDUE_LETTER_CONVERSION,
    SEQUENCE_IDENTITY_CUTOFF,
    STRIDE_CODES,
)
from pore.paths import PREPARED_PDB_DIR, ANNOTATED_PDB_DIR
from pore.types import VoxelGroup, Annotation
from pore import utils


PDB_IN = PDBParser()
PDB_OUT = PDBIO()


class ProteinSelect(Select):
    def __init__(self, components):
        Select.__init__(self)
        self.components = components

    def accept_model(self, model):
        if (model.serial_num == 0) or (model.serial_num == 1):
            return True
        return utils.KEEP_MODELS

    def accept_residue(self, residue):
        if residue.resname in self.components:
            return True
        return utils.KEEP_NON_PROTEIN

    def accept_atom(self, atom):
        if (not atom.is_disordered()) or (atom.get_altloc() == "A"):
            atom.set_altloc(" ")

        if atom.element != "H":
            return True
        return utils.KEEP_HYDROGENS


def save_pdb(structure: Structure.Structure, pdb_file: Path, remarks: Optional[str] = None) -> None:
    """
    Save a biopython PDB structure to a PDB file.
    """
    PDB_OUT.set_structure(structure)
    PDB_OUT.save(str(pdb_file))

    # if we are adding remarks, read in the just written file, add the remarks and overwrite
    if remarks is not None:
        with open(pdb_file, mode="r", encoding="utf-8") as in_file:
            lines = in_file.readlines()

        with open(pdb_file, mode="w", encoding="utf-8") as out_file:
            out_file.write(remarks + "".join(lines))


def load_pdb(pdb_file: Path, show_warnings: bool = False) -> Structure.Structure:
    """
    Load a PDB file as a biopython PDB structure.
    """
    if not show_warnings:
        warnings.filterwarnings("ignore", category=PDBConstructionWarning)
    structure = PDB_IN.get_structure(pdb_file.stem, pdb_file)

    # return the warnings back to their default state for future use
    warnings.filterwarnings("default", category=PDBConstructionWarning)

    return structure


def clean_structure(structure: Structure.Structure, protein_components: set[str]) -> Structure.Structure:
    """
    Clean a PDB by saving it to a temporary file with the select class and then reloading.
    """
    PDB_OUT.set_structure(structure)
    with tempfile.SpooledTemporaryFile(mode="w+") as pdb_file:
        PDB_OUT.save(pdb_file, select=ProteinSelect(protein_components))
        pdb_file.seek(0)
        structure = PDB_IN.get_structure(structure.id, pdb_file)

    return structure


def get_structure_coords(structure: Structure.Structure, canonical_protein_atoms_only: bool = False) -> pd.DataFrame:
    """
    Load the PDB file using Biopython.

    Return the coordinates of backbone heavy atoms and CB atoms as a dataframe.
    """
    if canonical_protein_atoms_only:
        coordinates = [atom.coord for atom in structure.get_atoms() if atom.name in VOXEL_ATOM_NAMES]
        elements = [atom.element for atom in structure.get_atoms() if atom.name in VOXEL_ATOM_NAMES]
    else:
        coordinates = [atom.coord for atom in structure.get_atoms()]
        elements = [atom.element for atom in structure.get_atoms()]

    coordinates = pd.DataFrame(coordinates, columns=["x", "y", "z"])
    coordinates["element"] = elements

    return coordinates


def make_atom_line(
    voxel_type: str,
    resnum: int,
    voxel_index: int,
    point: np.ndarray,
    beta: float,
) -> str:
    """
    Make a single atom line for a single voxel
    """
    return f"ATOM  {voxel_index:>5d}  {VOXEL_TYPE_ATOM_MAP[voxel_type]}   {voxel_type} {VOXEL_TYPE_CHAIN_MAP[voxel_type]}{resnum:>4d}    {point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00{beta:>6.2f}           {VOXEL_TYPE_ATOM_MAP[voxel_type]}"


def make_atom_lines(
    voxel_type: str,
    resnum: int,
    voxel_indices: set[int],
    surface_indices: set[int],
    voxel_grid_centers: np.ndarray,
) -> list[str]:
    """
    Make all atom lines for a given set of voxel indices
    """
    return [
        make_atom_line(voxel_type, resnum, voxel_index, voxel_grid_centers[voxel_index], 50.0)
        if voxel_index in surface_indices
        else make_atom_line(voxel_type, resnum, voxel_index, voxel_grid_centers[voxel_index], 0.0)
        for voxel_index in voxel_indices
    ]


def points_to_pdb(
    voxel_grid: VoxelGrid,
    hubs: dict[int, VoxelGroup],
    pores: dict[int, VoxelGroup],
    pockets: dict[int, VoxelGroup],
    cavities: dict[int, VoxelGroup],
    occluded: dict[int, VoxelGroup],
) -> list[str]:
    """
    Write out voxels of our volumes as though they were atoms in a PDB file.
    """
    return list(
        itertools.chain(
            *[
                make_atom_lines("HUB", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in hubs.items()
            ],
            *[
                make_atom_lines("POR", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in pores.items()
            ],
            *[
                make_atom_lines("POK", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in pockets.items()
            ],
            *[
                make_atom_lines("CAV", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in cavities.items()
            ],
            *[
                make_atom_lines("OCC", i, voxel_group.indices, voxel_group.surface_indices, voxel_grid.voxel_centers)
                for i, voxel_group in occluded.items()
            ],
        )
    )


def save_annotated_pdb(pdb_name: str, annotated_lines: list[str]) -> None:
    """
    Save a PDB formatted coordinate file of the voxels.
    Individual atoms/voxels are labelled according to type
    """
    # get the original PDB lines so that our annotation can be appended
    with open(PREPARED_PDB_DIR / f"{pdb_name}.pdb", mode="r", encoding="utf-8") as input_pdb_file:
        pdb_lines = input_pdb_file.readlines()
    pdb_lines = [line.rstrip("\n") for line in pdb_lines]

    # remove the terminal END line (NOTE: also removes terminal ENDMDL line)
    for line in pdb_lines[::-1]:
        if line.strip() == "":
            pdb_lines.pop()
        if ("ENDMDL" not in line) and ("END" in line):
            pdb_lines.pop()
            break

    # add the annotated lines and save
    with open(
        ANNOTATED_PDB_DIR / f"{pdb_name}.{str(utils.VOXEL_SIZE)}.pdb", mode="w", encoding="utf-8"
    ) as annotated_pdb_file:
        annotated_pdb_file.write("\n".join([*pdb_lines, "END", *annotated_lines]))


def generate_rotation_translation_remarks(rotation: np.ndarray, translation: np.ndarray) -> str:
    """
    Given a rotation matrix and translation vector, generate a string allowing them to be
    written as REMARKS to a PDB file.
    """
    return (
        f"REMARK TRANSLATION {str(translation)}\n"
        f"REMARK ROTATION [{str(rotation[0])}, {str(rotation[1])}, {str(rotation[2])}]\n"
    )


def generate_volume_remarks(annotation: Annotation) -> list[str]:
    """
    Return REMARKS for a PDB file to hold all volume information from an annotation object.
    """
    hub_remarks = [
        f"REMARK HUB VOLUME {annotation.hub_volumes[i]} DIMENSIONS {annotation.hub_dimensions[i][0], annotation.hub_dimensions[i][1], annotation.hub_dimensions[i][2]}"
        for i in annotation.hub_volumes
    ]

    pore_remarks = [
        f"REMARK PORE VOLUME {annotation.pore_volumes[i]} DIMENSIONS {annotation.pore_dimensions[i][0], annotation.pore_dimensions[i][1], annotation.pore_dimensions[i][2]}"
        for i in annotation.pore_volumes
    ]

    cavity_remarks = [
        f"REMARK CAVITY VOLUME {annotation.cavity_volumes[i]} DIMENSIONS {annotation.cavity_dimensions[i][0], annotation.cavity_dimensions[i][1], annotation.cavity_dimensions[i][2]}"
        for i in annotation.cavity_volumes
    ]

    pocket_remarks = [
        f"REMARK POCKET VOLUME {annotation.pocket_volumes[i]} DIMENSIONS {annotation.pocket_dimensions[i][0], annotation.pocket_dimensions[i][1], annotation.pocket_dimensions[i][2]}"
        for i in annotation.pocket_volumes
    ]

    return hub_remarks + pore_remarks + cavity_remarks + pocket_remarks


def generate_resolution_remarks() -> list[str]:
    """
    Add a REMARK line with the resolution used fo rthis annotation
    """
    return [f"REMARK RESOLUTION {str(utils.VOXEL_SIZE)}"]


def get_stoichiometry(pdb: Path, match_cutoff: float = SEQUENCE_IDENTITY_CUTOFF) -> dict[int, int]:
    """
    Example return = {1: 10, 2: 5} for a heteromultimer with 10 of one chain and 5 of the other
    """
    # load cleaned structure using biopython
    structure = load_pdb(pdb)

    # get the sequence of each using biopython
    sequences = []
    for chain in structure.get_chains():
        sequence = "".join([RESIDUE_LETTER_CONVERSION.get(residue.resname, "X") for residue in chain.get_residues()])
        sequences.append(sequence)

    # compute pairwise sequence identity using biopython to assign sequence clusters
    cluster_count = 0
    clusters = {cluster_count: [sequences.pop()]}
    while len(sequences) > 0:
        query_sequence = sequences.pop()
        joined_cluster = False
        for cluster_id, cluster_sequences in clusters.items():
            # only check against the first sequence in the cluster
            alignment = pairwise2.align.globalxx(query_sequence, cluster_sequences[0], one_alignment_only=True)[0]
            string_alignment = pairwise2.format_alignment(*alignment).split("\n")
            identity = string_alignment[1].count("|") / min(len(string_alignment[0]), len(string_alignment[2]))
            # if this is a match, add it to this cluster
            if identity >= match_cutoff:
                clusters[cluster_id].append(query_sequence)
                joined_cluster = True
                break

        # if we don't find any matches, start a new cluster
        if not joined_cluster:
            cluster_count += 1
            clusters[cluster_count] = [query_sequence]

    # compute the stoichiometry
    return {cluster_id: len(cluster_sequences) for cluster_id, cluster_sequences in clusters.items()}


def clean_stride(stride_line: str) -> str:
    """
    Cleans extra info from a STRIDE summary line
    """
    return stride_line[10:60].strip()


def get_secondary_structure(pdb: Path) -> dict[str, float]:
    """
    Compute the fraction of basic secondary structures, helix, strand, loop
    using the stride program.
    NOTE: requires local installation of STRIDE
    See http://webclu.bio.wzw.tum.de/stride/
    """
    process = subprocess.run(["stride", "-o", f"{pdb}"], capture_output=True)
    output = process.stdout.decode("UTF-8")

    primary_structure = "".join([clean_stride(line) for line in output.split("\n") if line[:3] == "SEQ"])
    secondary_structure = "".join([clean_stride(line) for line in output.split("\n") if line[:3] == "STR"])

    helix_count = sum([secondary_structure.count(code) for code in STRIDE_CODES["helix"]])
    strand_count = sum([secondary_structure.count(code) for code in STRIDE_CODES["strand"]])
    total_count = len(primary_structure)

    if total_count == 0:
        return {
            "helix": 0.0,
            "strand": 0.0,
            "coil": 1.0,
        }

    return {
        "helix": helix_count / total_count,
        "strand": strand_count / total_count,
        "coil": (total_count - helix_count - strand_count) / total_count,
    }
