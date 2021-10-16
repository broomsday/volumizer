"""
Module for interacting with the mongo database
"""


from pymongo import MongoClient, database

from pore.constants import MONGO_CONNECTION_STRING
from pore.types import PdbStatus, ComponentStatus, ComponentData


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


def get_init_pdb_entry(pdb: str) -> PdbStatus:
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


def get_init_component_entry(component: str) -> ComponentStatus:
    """
    Get an initial ComponentStatus for this pdb.
    """
    return ComponentStatus(
        component_id=component,
        processed=False,
        protein=None,
    )


def initialize_pdb_database(client: MongoClient, pdbs: set[str]) -> None:
    """
    Initialize the database with PDB IDs and blank status fields.
    """
    [client.status.pdbs.insert_one(get_init_pdb_entry(pdb)._asdict()) for pdb in pdbs]


def initialize_component_database(client: MongoClient, components: set[str]) -> None:
    """
    Initialize the database with PDB IDs and blank status fields.
    """
    [client.status.components.insert_one(get_init_component_entry(component)._asdict()) for component in components]


def update_database_pdbs(db: database.Database, pdbs: set[str]) -> None:
    """
    If a database exists, ensure it has the PDBs we want to work with.
    Delete PDBs that are not part of our set, and add those that are.
    """
    db_pdbs = set(db.pdbs.distinct("pdb_id"))

    pdbs_to_remove = db_pdbs - pdbs
    [db.pdbs.remove({"pdb_id": pdb}) for pdb in pdbs_to_remove]

    pdbs_to_add = pdbs - db_pdbs
    [db.pdbs.insert_one(get_init_pdb_entry(pdb)._asdict()) for pdb in pdbs_to_add]


def update_database_components(db: database.Database, components: set[str]) -> None:
    """
    If a database exists, ensure it has the components we want to work with.
    Delete components that are not part of our set, and add those that are.
    """
    db_components = set(db.components.distinct("component_id"))

    components_to_remove = db_components - components
    [db.components.remove({"component_id": component}) for component in components_to_remove]

    components_to_add = components - db_components
    [db.components.insert_one(get_init_component_entry(component)._asdict()) for component in components_to_add]


def quick_connect_database() -> database.Database:
    """
    Assume the database already exists and return it without updating.
    """
    return start_mongo_client().status


def fetch_database(pdbs: set[str], components: set[str]) -> database.Database:
    """
    If a database already exists, ensure it is up-to-date in terms of PDB IDs, and return it.

    If it does not exist, initialize it.
    """
    client = start_mongo_client()
    if not database_exists(client):
        initialize_pdb_database(client, pdbs)
        initialize_component_database(client, components)

    update_database_pdbs(client.status, pdbs)
    update_database_components(client.status, components)

    return client.status


def update_pdb_downloaded(db: database.Database, pdb: dict, downloaded: bool) -> None:
    """
    Update the downloaded field of the pdb collection
    """
    update_result = db.pdbs.update_one({"pdb_id": pdb["pdb_id"]}, {"$set": {"downloaded": downloaded}}, upsert=False)
    assert update_result.matched_count == 1


def update_pdb_cleaned(db: database.Database, pdb: dict, cleaned: bool) -> None:
    """
    Update the cleaned field of the pdb collection
    """
    update_result = db.pdbs.update_one({"pdb_id": pdb["pdb_id"]}, {"$set": {"cleaned": cleaned}}, upsert=False)
    assert update_result.matched_count == 1


def update_component_is_protein(db: database.Database, component: ComponentData, protein: bool) -> None:
    """
    Update the protein field of the components collection
    """
    update_result = db.components.update_one({"component_id": component.component_id}, {"$set": {"protein": protein}}, upsert=False)
    assert update_result.matched_count == 1

    update_result = db.components.update_one({"component_id": component.component_id}, {"$set": {"cleaned": True}}, upsert=False)
    assert update_result.matched_count == 1


def get_protein_components(db: database.Database) -> set[str]:
    """
    Return a set of all component IDs that are protein.
    """
    return {component["component_id"] for component in db.components.find({"protein": True})}
