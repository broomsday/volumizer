"""
Download biological assemblies from the RCSB and find those that are pores capable of supporting de novo enzyme design.
"""


from pathlib import Path

import typer

from pore import rcsb, pdb, mongo, paths


def main(
    cluster_file: Path = typer.Option(paths.RCSB_CLUSTER_FILE, help="Custom file containing one PDB ID per line"),
    skip_download: bool = typer.Option(False, help="Do not attempt to download any PDBs"),
    skip_process: bool = typer.Option(False, help="Do not process raw downloaded PDBs"),
) -> None:
    """
    Download biological assemblies, clean and process them, then find pores.
    """
    if cluster_file == paths.RCSB_CLUSTER_FILE:
        rcsb.get_rcsb_cluster_file()

    pdbs = rcsb.build_pdb_set(cluster_file)
    db = mongo.fetch_database(pdbs)

    if not skip_download:
        rcsb.download_biological_assemblies(db)

    print("Downloaded", db.pdbs.count_documents({"downloaded": True}))
    print("Not downloaded:", db.pdbs.count_documents({"downloaded": False}))

    if not skip_process:
        pdb.process_all_pdbs(db)


if "__main__" in __name__:
    typer.run(main)