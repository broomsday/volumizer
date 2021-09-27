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
    db = mongo.fetch_database(pdbs)
    rcsb.download_biological_assemblies(db)

    print("Downloaded", db.pdbs.find({"downloaded": True}).count())
    print("Not downloaded:", db.pdbs.find({"downloaded": False}).count())

    pdb.process_all_pdbs(db)


if "__main__" in __name__:
    typer.run(main)