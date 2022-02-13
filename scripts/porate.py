"""
Command-line entry-point to find pores and cavities in PDBs.
"""

from pathlib import Path

import typer

from pore import pore, cli, utils, pdb
from pore.constants import VOXEL_SIZE


def porate_pdb_id(pdb_id: str) -> None:
    """
    Download the given PDB ID and then porate it.
    """
    if utils.have_annotation(pdb_id):
        print(f"Already completed: {pdb_id}")
        annotation_df = utils.load_annotation_df(pdb_id)
    else:
        pdb_path = pore.download_pdb_file(pdb_id)
        if pdb_path is None:
            return None

        print(f"Working on: {pdb_id}")
        annotation, annotated_pdb_lines = pore.process_pdb_file(pdb_path)
        annotation_df = utils.make_annotation_dataframe(annotation)

        utils.save_annotation_dataframe(pdb_id, annotation_df)
        pdb.save_annotated_pdb(pdb_id, annotated_pdb_lines)

    print(annotation_df)


def porate_pdb_file(pdb_file: Path) -> None:
    """
    Operate directly on the given PDB file.
    """
    if utils.have_annotation(pdb_file.stem):
        print(f"Already completed: {pdb_file.stem}")
        annotation_df = utils.load_annotation_df(pdb_file.stem)
    else:
        print(f"Working on: {pdb_file.stem}")
        annotation, annotated_pdb_lines = pore.process_pdb_file(pdb_file)
        annotation_df = utils.make_annotation_dataframe(annotation)

        utils.save_annotation_dataframe(Path(pdb_file).stem, annotation_df)
        pdb.save_annotated_pdb(Path(pdb_file).stem, annotated_pdb_lines)

    print(annotation_df)


def main(
    porate_input: str = typer.Argument(
        ..., help="PDB ID, PDB file, file with one PDB ID per line, or folder containing PDB files"
    ),
    resolution: float = typer.Option(VOXEL_SIZE, help="Edge-length of voxels used to discretize the structure."),
    non_protein: bool = typer.Option(False, help="Include non-protein residues during the calculations."),
):
    """
    Find pores and cavities in the supplied PDB files.
    """

    utils.setup_dirs()
    utils.ensure_protein_components_file()

    utils.set_resolution(resolution)
    utils.set_non_protein(non_protein)
    input_type = cli.guess_input_type(porate_input)

    if input_type == "pdb_id":
        porate_pdb_id(porate_input)
    elif input_type == "pdb_file":
        pdb_file = Path(porate_input)
        porate_pdb_file(pdb_file)
    elif input_type == "id_file":
        with open(porate_input, mode="r", encoding="utf-8") as id_file:
            pdb_ids = [line.strip() for line in id_file.readlines()]
        for pdb_id in pdb_ids:
            porate_pdb_id(pdb_id)
    elif input_type == "pdb_dir":
        pdb_files = Path(porate_input).glob("*.pdb")
        for pdb_file in pdb_files:
            porate_pdb_file(pdb_file)
    else:
        raise RuntimeError("File mode not implemented")


if "__main__" in __name__:
    typer.run(main)
