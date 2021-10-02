"""
Give some summary information on the database
"""


import typer

from pore import mongo


def main():
    """
    Report the number of downloaded, cleaned, and processed PDBs.
    Report the nubmer of PDBs that are pores.
    Report the number of downloaded, and processed components.
    Report the number of components that are protein.
    """
    client = mongo.start_mongo_client()

    print(f"Downloaded: {client.status.pdbs.count_documents({'downloaded': True})} PDBs")
    print(f"Cleaned: {client.status.pdbs.count_documents({'cleaned': True})} PDBs")

    print(f"\nProcessed: {client.status.components.count_documents({'processed': True})} components")
    print(f"{client.status.components.count_documents({'protein': True})} components are protein")

    print(f"\nProcessed: {client.status.pdbs.count_documents({'processed': True})} PDBs")
    print(f"{client.status.pdbs.count_documents({'pore': True})} PDBs are pores")


if "__main__" in __name__:
    typer.run(main)