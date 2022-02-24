"""
Functions to interpret CLI commands.
"""


from pathlib import Path
from typing import Optional


def guess_input_type(input: str) -> Optional[str]:
    """
    Based on the input string, guess if this input is:

    1. a single PDB ID
    2. a single PDB file
    3. a text file containing multiple PDB IDs
    4. a directory containing multiple PDB files.
    """

    if Path(input).is_file():
        if "pdb" in Path(input).suffix:
            return "pdb_file"
        else:
            return "id_file"
    elif Path(input).is_dir():
        return "pdb_dir"
    elif len(input) == 4:
        return "pdb_id"
    else:
        return None


def guess_analysis_input_type(input: str) -> Optional[str]:
    """
    Based on the input string, guess if this input is:

    1. a text containing PDB IDs
    2. a directory containing Dataframes
    """

    if Path(input).is_file():
        return "file"
    elif Path(input).is_dir():
        return "dir"
    else:
        return None