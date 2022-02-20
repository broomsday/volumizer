"""
Helper script to take a list of PDB IDs and return only those whose size or related metrics
fall within given ranges.
"""


from pathlib import Path
from time import sleep
import warnings

import typer
from tqdm import tqdm

from pore import analysis, rcsb, utils


def main(
    input_list: Path = typer.Argument(..., help=""),
    output_list: Path = typer.Argument(..., help=""),
    min_atoms: int = typer.Option(1, help=""),
    max_atoms: int = typer.Option(None, help=""),
    min_residues: int = typer.Option(1, help=""),
    max_residues: int = typer.Option(None, help=""),
    min_chains: int = typer.Option(1, help=""),
    max_chains: int = typer.Option(None, help=""),
):
    """
    Subset a PDB list based on some size metrics
    """
    # we'll be saving some data so make sure directories are available
    utils.setup_dirs()

    metrics = {
        "min_atoms": min_atoms,
        "max_atoms": max_atoms,
        "min_residues": min_residues,
        "max_residues": max_residues,
        "min_chains": min_chains,
        "max_chains": max_chains,
    }

    # get the list of PDBs we need to check
    with open(input_list, mode="r", encoding="utf-8") as in_file:
        # NOTE: split around '.' to ignore any resolution suffixes
        pdb_ids = [line.rstrip().split(".")[0] for line in in_file.readlines()]

    # check the PDBs
    satisfied_pdb_ids = []
    for pdb_id in tqdm(pdb_ids):
        if utils.have_pdb_size_metrics_on_file(pdb_id):
            pdb_size_metrics = utils.load_pdb_size_metrics(pdb_id)
        else:
            pdb_size_metrics = rcsb.get_pdb_size_metrics(pdb_id)
            utils.save_pdb_size_metrics(pdb_id, pdb_size_metrics)
            sleep(0.25)

        if pdb_size_metrics is None:
            warnings.warn(f"No PDB size metrics for: {pdb_id}")
        elif analysis.pdb_satisfies_metrics(pdb_size_metrics, metrics):
            satisfied_pdb_ids.append(pdb_id)

    with open(output_list, mode="w", encoding="utf-8") as out_file:
        out_file.writelines([f"{pdb}\n" for pdb in satisfied_pdb_ids])

    print(f"Original number of PDBs: {len(pdb_ids)}")
    print(f"Final number of PDBs: {len(satisfied_pdb_ids)}")


if "__main__" in __name__:
    typer.run(main)