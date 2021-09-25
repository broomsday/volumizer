"""
Module for interacting with the mongo database
"""


from pymongo import MongoClient, database

from pore.constants import MONGO_CONNECTION_STRING
from pore.types import PdbStatus


def start_mongo_client() -> MongoClient:
    """
    Connect to the mongodb.
    """
    return MongoClient(MONGO_CONNECTION_STRING)


def database_exists(client: MongoClient) -> bool:
    """
    Check to see if we have already created the database.
    """
    return client.status.command("dbstats")["objects"] > 0


def get_init_entry(pdb: str) -> PdbStatus:
    """
    Get an initial PdbStatus for this pdb.
    """
    return PdbStatus(
        pdb_id=pdb,
        downloaded=False,
        cleaned=False,
        processed=False,
        pore=None,
    )


def initialize_database(client: MongoClient, pdbs: set[str]) -> None:
    """
    Initialize the database with PDB IDs and blank status fields.
    """
    [client.status.pdbs.insert_one(get_init_entry(pdb)._asdict()) for pdb in pdbs]


def update_database_pdbs(db: database.Database, pdbs: set[str]) -> None:
    """
    If a database exists, ensure it has the PDBs we want to work with.
    Delete PDBs that are not part of our set, and add those that are.
    """
    db_pdbs = set(db.pdbs.distinct("pdb_id"))

    pdbs_to_remove = db_pdbs - pdbs
    [db.pdbs.remove({"pdb_id": pdb}) for pdb in pdbs_to_remove]

    pdbs_to_add = pdbs - db_pdbs
    [db.pdbs.insert_one(get_init_entry(pdb)._asdict()) for pdb in pdbs_to_add]


def fetch_database(pdbs: set[str]) -> database.Database:
    """
    If a database already exists, ensure it is up-to-date in terms of PDB IDs, and return it.

    If it does not exist, initialize it.
    """
    client = start_mongo_client()
    if not database_exists(client):
        initialize_database(client, pdbs)

    update_database_pdbs(client.status, pdbs)

    return client.status


def update_downloaded(db: database.Database, pdb: dict, downloaded: bool) -> None:
    """
    Update the downloaded field 
    """
    update_result = db.pdbs.update_one({"pdb_id": pdb["pdb_id"]}, {"$set": {"downloaded": downloaded}}, upsert=False)
    assert update_result.matched_count == 1
