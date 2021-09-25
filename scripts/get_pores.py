"""
Download biological assemblies from the RCSB and find those that are pores capable of supporting de novo enzyme design.
"""


from pathlib import Path

import typer

from pore import rcsb, pdb, mongo, paths


def main(
    cluster_file: Path = typer.Option(paths.RCSB_CLUSTER_FILE, help="Custom file containing one PDB ID per line")
) -> None:
    """
    Download biological assemblies, clean and process them, then find pores.
    """
    if cluster_file == paths.RCSB_CLUSTER_FILE:
        rcsb.get_rcsb_cluster_file()

    pdbs = rcsb.build_pdb_set(cluster_file)
    print(pdbs)


if "__main__" in __name__:
    typer.run(main)