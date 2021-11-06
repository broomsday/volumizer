"""
Definition of custom types.
"""


from typing import NamedTuple


class ComponentData(NamedTuple):
    component_id: str
    component_type: str


class Annotation(NamedTuple):
    total_pore_volume: float
    total_cavity_volume: float
    largest_pore_volume: float
    largest_cavity_volume: float
    pore_volumes: dict[int, float]
    cavity_volumes: dict[int, float]
