"""
Command-line entry-point to find pores and cavities in PDBs.
"""

from pathlib import Path

import typer

from pore import pore, cli, utils


def main(
    input: str = typer.Argument(..., help="PDB ID, PDB file, file a PDB ID per line, or folder containing PDB files")
):
    """
    Find pores and cavities in the supplied PDB files.
    """

    utils.setup_dirs()
    utils.ensure_protein_components_file()

    input_type = cli.guess_input_type(input)

    if input_type == "pdb_id":
        pdb_path = pore.download_pdb_file(input)
        annotation = pore.process_pdb_file(pdb_path)
        print(annotation)
    elif input_type == "pdb_file":
        pdb_path = Path(input)
        annotation = pore.process_pdb_file(pdb_path)
        print(annotation)
    else:
        raise RuntimeError("File mode not implemented")


if "__main__" in __name__:
    typer.run(main)