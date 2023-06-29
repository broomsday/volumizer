"""
For each PDB file in a list, read in the file, determine the sequence of all chains
and then guess the stoichiometry based on sequence alignment.
"""


from pathlib import Path
import warnings

import typer
from tqdm import tqdm

from volumizer import utils, pdb, paths, analysis


def main(
    input_list: Path = typer.Argument(..., help=""),
    output_list: Path = typer.Argument(..., help=""),
    min_chain_repeats: int = typer.Option(1, help=""),
    max_chain_repeats: int = typer.Option(None, help=""),
    min_unique_chains: int = typer.Option(1, help=""),
    max_unique_chains: int = typer.Option(None, help=""),
    stoichiometry_factorable: bool = typer.Option(
        False,
        help="If True, only structures where stoichiometry for all chains are factors of one another can pass, e.g. 8-4-2",
    ),
):
    """
    For each PDB file in a list, read in the file, determine the sequence of all chains
    and then guess the stoichiometry based on sequence alignment.
    """
    # we'll be saving some data so make sure directories are available
    utils.setup_dirs()

    metrics = {
        "min_chain_repeats": min_chain_repeats,
        "max_chain_repeats": max_chain_repeats,
        "min_unique_chains": min_unique_chains,
        "max_unique_chains": max_unique_chains,
        "stoichiometry_factorable": stoichiometry_factorable,
    }

    # get the list of PDBs we need to check
    with open(input_list, mode="r", encoding="utf-8") as in_file:
        # NOTE: split around '.' to ignore any resolution suffixes
        pdb_ids = [line.rstrip().split(".")[0] for line in in_file.readlines()]

    satisfied_pdb_ids = []
    for pdb_id in tqdm(pdb_ids):
        if utils.have_stoichiometry_on_file(pdb_id):
            stoichiometry = utils.load_stoichiometry(pdb_id)
        else:
            prepared_pdb = paths.PREPARED_PDB_DIR / f"{pdb_id}.pdb"
            stoichiometry = pdb.get_stoichiometry(prepared_pdb)
            utils.save_stoichiometry(pdb_id, stoichiometry)

        if stoichiometry is None:
            warnings.warn(f"No stoichiometry: {pdb_id}")
        elif analysis.pdb_satisfies_stoichiometry(stoichiometry, metrics):
            satisfied_pdb_ids.append(pdb_id)

    with open(output_list, mode="w", encoding="utf-8") as out_file:
        out_file.writelines([f"{pdb}\n" for pdb in satisfied_pdb_ids])

    print(f"Original number of PDBs: {len(pdb_ids)}")
    print(f"Final number of PDBs: {len(satisfied_pdb_ids)}")


if "__main__" in __name__:
    typer.run(main)