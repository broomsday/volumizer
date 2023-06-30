"""
For each PDB file in a list, read in the file, determine the secondary structure
and return a list of only those matching given secondary structure ranges.
"""


from pathlib import Path
import warnings

import typer
from tqdm import tqdm

from volumizer import utils, pdb, paths, analysis


def main(
    input_list: Path = typer.Argument(..., help=""),
    output_list: Path = typer.Argument(..., help=""),
    min_helix: float = typer.Option(0.0, help=""),
    max_helix: float = typer.Option(1.0, help=""),
    min_strand: float = typer.Option(0.0, help=""),
    max_strand: float = typer.Option(1.0, help=""),
    min_coil: float = typer.Option(0.0, help=""),
    max_coil: float = typer.Option(1.0, help=""),
):
    """
    For each PDB file in a list, read in the file, determine the secondary structure
    and return a list of only those matching given secondary structure ranges.
    """
    # we'll be saving some data so make sure directories are available
    utils.setup_dirs()

    metrics = {
        "min_helix": min_helix,
        "max_helix": max_helix,
        "min_strand": min_strand,
        "max_strand": max_strand,
        "min_coil": min_coil,
        "max_coil": max_coil,
    }

    # get the list of PDBs we need to check
    with open(input_list, mode="r", encoding="utf-8") as in_file:
        # NOTE: split around '.' to ignore any resolution suffixes
        pdb_ids = [line.rstrip().split(".")[0] for line in in_file.readlines()]

    satisfied_pdb_ids = []
    for pdb_id in tqdm(pdb_ids):
        if utils.have_secondary_structure_on_file(pdb_id):
            secondary_structure = utils.load_secondary_structure(pdb_id)
        else:
            prepared_pdb = paths.PREPARED_PDB_DIR / f"{pdb_id}.pdb"
            prepared_structure = pdb.load_structure(prepared_pdb)
            secondary_structure = pdb.get_secondary_structure(prepared_structure)
            utils.save_secondary_structure(pdb_id, secondary_structure)

        if secondary_structure is None:
            warnings.warn(f"No secondary structure: {pdb_id}")
        elif analysis.pdb_satisfies_secondary_structure(secondary_structure, metrics):
            satisfied_pdb_ids.append(pdb_id)

    with open(output_list, mode="w", encoding="utf-8") as out_file:
        out_file.writelines([f"{pdb}\n" for pdb in satisfied_pdb_ids])

    print(f"Original number of PDBs: {len(pdb_ids)}")
    print(f"Final number of PDBs: {len(satisfied_pdb_ids)}")


if "__main__" in __name__:
    typer.run(main)