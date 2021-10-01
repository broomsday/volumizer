"""
Functions for parsing, cleaning, and modifying PDBs.
"""

import warnings

from pymongo import database
from tqdm import tqdm

from Bio.PDB import PDBParser, Select
from Bio.PDB.PDBIO import PDBIO

from pore import utils, mongo
from pore.paths import CLEANED_PDB_DIR


PDB_IN = PDBParser()
PDB_OUT = PDBIO()


class ProteinSelect(Select):
    def __init__(self, components):
        Select.__init__(self)
        self.components = components

    def accept_model(self, model):
        if model.serial_num == 1:
            return True
        return False
    def accept_residue(self, residue):
        if residue.resname in self.components:
            return True
        return False
    def accept_atom(self, atom):
        if atom.element != "H":
            return True
        if (not atom.is_disordered()) or (atom.get_altloc() == "A"):
            atom.set_altloc(" ")
            return True
        return False


def clean_one_pdb(pdb_id: str, protein_components: set[str]) -> bool:
    """
    Decompress and cleanup a single PDB.

    Return False if the final processed PDB could not be created.
    """
    utils.decompress_pdb(pdb_id)
    structure = PDB_IN.get_structure(pdb_id, CLEANED_PDB_DIR / f"{pdb_id}.pdb")

    PDB_OUT.set_structure(structure)
    PDB_OUT.save(str(CLEANED_PDB_DIR / f"{pdb_id}.pdb"), select=ProteinSelect(protein_components))

    return utils.is_pdb_cleaned(pdb_id)


def clean_all_pdbs(db: database.Database) -> None:
    """
    Decompress and cleanup all downloaded PDBs, and align to their major axis.
    
    Enter their state in the database.
    """
    CLEANED_PDB_DIR.mkdir(exist_ok=True, parents=True)

    protein_components = mongo.get_protein_components(db)

    # biopython tends to have many warnings about PDB construction which can be ignored
    warnings.filterwarnings("ignore")
    for pdb in tqdm(list(db.pdbs.find({"downloaded": True})), "Cleaning PDBs"):
        if utils.is_pdb_cleaned(pdb["pdb_id"]):
            mongo.update_pdb_cleaned(db, pdb, True)
        else:
            mongo.update_pdb_cleaned(db, pdb, clean_one_pdb(pdb["pdb_id"], protein_components))
    warnings.filterwarnings("default")