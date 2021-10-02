"""
Drop a collection from the status database.

This is intended to be used as a development tool when the collection has become dirty.
"""


import typer

from pore.constants import MONGO_CONNECTION_STRING
from pore.mongo import start_mongo_client


def main(
    collection: str = typer.Argument(..., help="The collection you with to reset [pdbs, components]"),
) -> None:
    """
    Completey drop the given collection from the status database.
    """
    client = start_mongo_client()
    client.status.drop_collection(collection)


if "__main__" in __name__:
    typer.run(main)