"""
Command-line entry-point to find pores and cavities in PDBs.
"""

from pathlib import Path

import typer

from pore import pore, cli, utils
from pore.paths import ANNOTATED_PDB_DIR


def main(
    input: str = typer.Argument(..., help="PDB ID, PDB file, file with one PDB ID per line, or folder containing PDB files"),
    resolution: float = typer.Option(2.0, help="Edge-length of voxels used to discretize the structure.")
):
    """
    Find pores and cavities in the supplied PDB files.
    """

    utils.setup_dirs()
    utils.ensure_protein_components_file()

    utils.set_resolution(resolution)

    input_type = cli.guess_input_type(input)

    if input_type == "pdb_id":
        pdb_path = pore.download_pdb_file(input)
        annotation, pdb_lines = pore.process_pdb_file(pdb_path)
        utils.print_annotation(annotation)
        utils.save_annotated_pdb(input, pdb_lines)
    elif input_type == "pdb_file":
        pdb_path = Path(input)
        annotation, pdb_lines = pore.process_pdb_file(pdb_path)
        utils.print_annotation(annotation)
        utils.save_annotated_pdb(Path(input).stem, pdb_lines)
    else:
        raise RuntimeError("File mode not implemented")




if "__main__" in __name__:
    typer.run(main)