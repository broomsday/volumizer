"""
Download biological assemblies from the RCSB and find those that are pores capable of supporting de novo enzyme design.
"""


from pathlib import Path

import typer

from pore import rcsb, pdb, mongo, paths


def main(
    cluster_file: Path = typer.Option(paths.RCSB_CLUSTER_FILE, help="Custom file containing one PDB ID per line"),
    component_file: Path = typer.Option(paths.RCSB_CCD_FILE, help="Custom file containing one component name per line.  Only these components will be included in the cleaned PDBs."),
    skip_download: bool = typer.Option(False, help="Do not attempt to download any PDBs"),
    skip_clean: bool = typer.Option(False, help="Do not clean raw downloaded PDBs"),
    quick_connect: bool = typer.Option(False),
) -> None:
    """
    Download biological assemblies, clean and process them, then find pores.
    """
    if cluster_file == paths.RCSB_CLUSTER_FILE:
        print("Downloading cluster file...")
        rcsb.get_cluster_file()
    if component_file == paths.RCSB_CCD_FILE:
        print("Downloading component file...")
        rcsb.get_component_file()

    if quick_connect:
        db = mongo.quick_connect_database()
        components = []
    else:
        print("Building PDB set...")
        pdbs = rcsb.build_pdb_set(cluster_file)
        print("Building component set...")
        components = rcsb.build_component_set(component_file)
        print("Getting the database...")
        db = mongo.fetch_database(pdbs, set([component.component_id for component in components]))

    if not skip_download:
        rcsb.download_biological_assemblies(db)

    print(f"Downloaded: {db.pdbs.count_documents({'downloaded': True})} PDBs")
    print(f"Not downloaded: {db.pdbs.count_documents({'downloaded': False})} PDBs")

    if not skip_clean:
        rcsb.process_all_components(db, components)
        pdb.clean_all_pdbs(db)

    pdb.process_all_pdbs(db)



if "__main__" in __name__:
    typer.run(main)