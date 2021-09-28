"""
Functions for parsing, cleaning, and modifying PDBs.
"""


from pymongo import database
from tqdm import tqdm

from pore import utils, mongo, rcsb


def assemble_biological_unit():
    """
    Alter the chain IDs in order to make the biological unit follow expected PDB rules.
    """


def remove_alternate_conformers():
    """
    Remove alternate conformers and remove the alt-loc identifier from the remaining conformer.
    """


def remove_hydrogens():
    """
    Remove hydrogen atoms.
    """


def remove_non_protein():
    """
    Remove anything that isn't an L- or D- protein component
    TODO: check RCSB components to get these definitions
    """


def renumber_residues():
    """
    Renumber residues, collapsing fake gaps but maintaining true ones.
    """


def process_one_pdb(pdb_id: str) -> bool:
    """
    Decompress and cleanup a single PDB.

    Return False if the final processed PDB could not be created.
    """
    # TODO: assemble the biological unit in a temporary space
    # TODO: remove waters, salts and truly non-protein residues (TODO: need list of protein residues from RCSB-components)
    # TODO: remove alternate conformers
    # TODO: remove hydrogens
    # TODO: renumber the residues

    return False


def process_all_pdbs(db: database.Database) -> None:
    """
    Decompress and cleanup all downloaded PDBs, and align to their major axis.
    
    Enter their state in the database.
    """
    # TODO: get an initial set of protein components
    #protein_residue_names = rcsb.
    print(db.command("dbstats"))
    quit()
    for pdb in tqdm(list(db.pdbs.find()), "Processing PDBs"):
        if utils.is_pdb_processed(pdb["pdb_id"]):
            mongo.update_processed(db, pdb, True)
        else:
            mongo.update_downloaded(db, pdb, process_one_pdb(pdb["pdb_id"]))