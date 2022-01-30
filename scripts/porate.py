"""
Command-line entry-point to find pores and cavities in PDBs.
"""

from pathlib import Path
from typing import Optional

import typer
import pandas as pd

from pore import pore, cli, utils, pdb
from pore.types import Annotation


def porate_pdb_id(pdb_id: str) -> tuple[Optional[Annotation], Optional[pd.DataFrame]]:
    """
    Download the given PDB ID and then porate it.
    """
    pdb_path = pore.download_pdb_file(pdb_id)
    if pdb_path is None:
        return None, None

    annotation, annotated_pdb_lines = pore.process_pdb_file(pdb_path)
    annotation_df = utils.make_annotation_dataframe(annotation)

    utils.save_annotation_dataframe(pdb_id, annotation_df)
    pdb.save_annotated_pdb(pdb_id, annotated_pdb_lines)

    return annotation, annotation_df


def porate_pdb_file(pdb_file: Path) -> tuple[Annotation, pd.DataFrame]:
    """
    Operate directly on the given PDB file.
    """
    annotation, annotated_pdb_lines = pore.process_pdb_file(pdb_file)
    annotation_df = utils.make_annotation_dataframe(annotation)

    utils.save_annotation_dataframe(Path(pdb_file).stem, annotation_df)
    pdb.save_annotated_pdb(Path(pdb_file).stem, annotated_pdb_lines)

    return annotation, annotation_df


def main(
    porate_input: str = typer.Argument(
        ..., help="PDB ID, PDB file, file with one PDB ID per line, or folder containing PDB files"
    ),
    resolution: float = typer.Option(2.0, help="Edge-length of voxels used to discretize the structure."),
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
        annotation, annotation_df = porate_pdb_id(porate_input)
        utils.print_annotation(annotation, annotation_df)
    elif input_type == "pdb_file":
        pdb_file = Path(porate_input)
        annotation, annotation_df = porate_pdb_file(pdb_file)
        utils.print_annotation(annotation, annotation_df)
    elif input_type == "id_file":
        with open(porate_input, mode="r", encoding="utf-8") as id_file:
            ids = [line.strip() for line in id_file.readlines()]
        for id in ids:
            print(id)
            annotation, annotation_df = porate_pdb_id(id)
            utils.print_annotation(annotation, annotation_df)
    elif input_type == "pdb_dir":
        pdb_files = Path(porate_input).glob("*.pdb")
        for pdb_file in pdb_files:
            print(pdb_file.stem)
            annotation, annotation_df = porate_pdb_file(pdb_file)
            utils.print_annotation(annotation, annotation_df)
    else:
        raise RuntimeError("File mode not implemented")


if "__main__" in __name__:
    typer.run(main)
