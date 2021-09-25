"""
Definition of custom types.
"""


from typing import NamedTuple, Optional


class PdbStatus(NamedTuple):
    pdb_id: str
    downloaded: bool
    cleaned: bool
    processed: bool
    pore: Optional[bool]
