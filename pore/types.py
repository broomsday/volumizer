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


class ComponentData(NamedTuple):
    component_id: str
    component_type: str


class ComponentStatus(NamedTuple):
    component_id: str
    processed: bool
    protein: Optional[bool]
