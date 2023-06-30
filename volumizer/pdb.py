"""
Functions for parsing, cleaning, and modifying PDBs.
"""

from pathlib import Path
import itertools
from typing import Optional
import subprocess

from Bio import pairwise2
import biotite.structure as bts
from biotite.structure.io import load_structure, save_structure
from pyntcloud.structures.voxelgrid import VoxelGrid
import pandas as pd
import numpy as np

from volumizer.constants import (
    VOXEL_ATOM_NAMES,
    VOXEL_TYPE_CHAIN_MAP,
    VOXEL_TYPE_ATOM_MAP,
    RESIDUE_LETTER_CONVERSION,
    SEQUENCE_IDENTITY_CUTOFF,
    STRIDE_CODES,
)
from volumizer.paths import PREPARED_PDB_DIR, ANNOTATED_PDB_DIR
from volumizer.types import VoxelGroup, Annotation
from volumizer import utils


def save_pdb_lines(pdb_lines: list[str], output: Path) -> None:
    """
    Save individual pdb lines to a file.
    """
    with open(output, mode="w", encoding="utf-8") as out_file:
        out_file.writelines("%s\n" % line for line in pdb_lines)


def save_pdb(structure: bts.AtomArray, pdb_file: Path, remarks: Optional[str] = None) -> None:
    """
    Save a biotite structure to a PDB file.
    """
    # TODO: simplify or even remove this function, should just use `save_pdb_lines` or biotite's `save_structure`
    save_structure(pdb_file, structure)

    # if we are adding remarks, read in the just written file, add the remarks and overwrite
    if remarks is not None:
        with open(pdb_file, mode="r", encoding="utf-8") as in_file:
            lines = in_file.readlines()

        with open(pdb_file, mode="w", encoding="utf-8") as out_file:
            out_file.write(remarks + "".join(lines))


def load_pdb(pdb_file: Path) -> bts.AtomArray:  # TODO: can probably just delete this and import `load_structure` directly from `pdb` in other functions
    """
    Load a PDB file as a biotite AtomArray.
    """
    return load_structure(pdb_file)


def clean_structure(structure: bts.AtomArray, protein_components: set[str]) -> bts.AtomArray:
    """
    Clean the AtomArray of a PDB based on selected preferences.
    """
    if isinstance(structure, bts.AtomArrayStack):
        if not utils.KEEP_MODELS:
            structure = structure[0]
        elif structure.stack_depth() > 1:
            new_structure = structure[0]
            for model in structure[1:]:
                new_structure += model
            structure = new_structure

    if not utils.KEEP_NON_PROTEIN:
        structure = structure[np.isin(structure.res_name, list(protein_components))]

    if not utils.KEEP_HYDROGENS:
        structure = structure[~np.isin(structure.element, ["H"])]

    return structure


def get_structure_coords(structure: bts.AtomArray, canonical_protein_atoms_only: bool = False) -> pd.DataFrame:
    """
    Return the coordinates of an atom array as a dataframe.

    # TODO: shouldn't we just keep the coords as an np array?  why are we converting to a less performant format?
    # TODO: in fact, we could just keep the whole structure and thereby keep the elements along with it
    """
    if canonical_protein_atoms_only:
        structure = structure[np.isin(structure.atom_name, VOXEL_ATOM_NAMES)]

    coordinates = structure.coord
    elements = structure.element

    coordinates = pd.DataFrame(coordinates, columns=["x", "y", "z"])
    coordinates["element"] = elements

    return coordinates


def structure_to_lines(structure: bts.AtomArray) -> list[str]:
    """
    Convert an AtomArray into individual PDB lines
    """
    try:
        return [make_atom_line(atom.atom_name, atom.res_name, atom.res_id, atom.chain_id, idx, atom.coord, atom.element, beta=atom.b_factor) for idx, atom in enumerate(structure)]
    except AttributeError:
        return [make_atom_line(atom.atom_name, atom.res_name, atom.res_id, atom.chain_id, idx, atom.coord, atom.element) for idx, atom in enumerate(structure)]


def format_atom_name(atom_name: str) -> str:
    """
    Format an atom-name to occupy 4 characters the weird way PDBs do it.
    """
    if len(atom_name) == 1:
        return f" {atom_name}  "
    elif len(atom_name) == 2:
        return f" {atom_name} "
    elif len(atom_name) == 3:
        return f" {atom_name}"
    return atom_name[:4]


def make_atom_line(
    atom_name: str,
    res_name: str,
    res_num: int,
    chain_id: str,
    atom_index: int,
    coord: np.ndarray,
    element: str,
    beta: float = 0.0,
) -> str:
    """
    Make a single atom line for a single voxel
    """
    return f"ATOM  {atom_index:>5d}  {format_atom_name(atom_name)}{res_name} {chain_id}{res_num:>4d}    {coord[0]:>8.3f}{coord[1]:>8.3f}{coord[2]:>8.3f}  1.00{beta:>6.2f}           {element}"
 

def make_atom_lines(    # TODO: rename this to indicate it's specifically for voxels
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
        make_atom_line(VOXEL_TYPE_ATOM_MAP[voxel_type], voxel_type, resnum, VOXEL_TYPE_CHAIN_MAP[voxel_type], voxel_index, voxel_grid_centers[voxel_index], VOXEL_TYPE_ATOM_MAP[voxel_type], 50.0)
        if voxel_index in surface_indices
        else make_atom_line(VOXEL_TYPE_ATOM_MAP[voxel_type], voxel_type, resnum, VOXEL_TYPE_CHAIN_MAP[voxel_type], voxel_index, voxel_grid_centers[voxel_index], VOXEL_TYPE_ATOM_MAP[voxel_type], 0.0)
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


def save_annotated_pdb(prepared_pdb_lines: list[str], annotation_lines: list[str], output_file: Path) -> None:
    """
    Save a PDB formatted coordinate file of the voxels plus the prepared input pdb lines.
    Individual atoms/voxels are labelled according to type
    """
    # TODO: possibly get rid of this in favour of several smaller functions (merge_pdb_lines and save_pdb_lines)
    prepared_pdb_lines = [line.rstrip("\n") for line in prepared_pdb_lines]

    # remove the terminal END line from the pdb.
    for line in prepared_pdb_lines[::-1]:
        if line.strip() == "":
            prepared_pdb_lines.pop()
        if ("ENDMDL" not in line) and ("END" in line):
            prepared_pdb_lines.pop()
            break

    full_lines = "\n".join([*prepared_pdb_lines, "END", *annotation_lines])
    with open(output_file, mode="w", encoding="utf-8") as annotated_pdb_file:
        annotated_pdb_file.write(full_lines)

    # add the annotated lines and save
    #with open(
    #    ANNOTATED_PDB_DIR / f"{pdb_name}.{str(utils.VOXEL_SIZE)}.pdb", mode="w", encoding="utf-8"
    #) as annotated_pdb_file:
    #    annotated_pdb_file.write("\n".join([*pdb_lines, "END", *annotation_lines]))


def generate_rotation_translation_remarks(rotation: np.ndarray, translation: np.ndarray) -> list[str]:
    """
    Given a rotation matrix and translation vector, generate a string allowing them to be
    written as REMARKS to a PDB file.
    """
    return [
        f"REMARK TRANSLATION {str(translation)}",
        f"REMARK ROTATION [{str(rotation[0])}, {str(rotation[1])}, {str(rotation[2])}]",
    ]


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


def merge_pdb_lines(preparation_remarks: list[str], prepared_pdb: list[str], annotation_lines: list[str]) -> list[str]:
    """
    Merge together the various PDB-line elements into the final output file.
    """
    final_lines = preparation_remarks
    final_lines.extend(prepared_pdb)
    final_lines.extend(annotation_lines)

    return final_lines


def get_stoichiometry(pdb: Path, match_cutoff: float = SEQUENCE_IDENTITY_CUTOFF) -> dict[int, int]:
    """
    Example return = {1: 10, 2: 5} for a heteromultimer with 10 of one chain and 5 of the other
    """
    # TODO: this function and some below don't need to be in main volumizer package, are only used in helper scripts
    structure = load_pdb(pdb)

    sequences = []
    for chain in bts.chain_iter(structure):
        sequence = "".join([RESIDUE_LETTER_CONVERSION.get(residue.res_name[0], "X") for residue in bts.residue_iter(chain)])
        sequences.append(sequence)

    # compute pairwise sequence identity using biopython to assign sequence clusters
    cluster_count = 0
    clusters = {cluster_count: [sequences.pop()]}
    while len(sequences) > 0:
        query_sequence = sequences.pop()
        joined_cluster = False
        for cluster_id, cluster_sequences in clusters.items():
            # only check against the first sequence in the cluster
            # TODO: replace these functions with biotite?
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
    # TODO: replace this with biotite?
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
