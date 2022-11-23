#!/usr/bin/env python3
from typing import List, Dict

FIELD_TYPES: Dict[str, int] = {'u': 0, 'v': 1, 'w': 2, 'p': 3, 'T': 4, 'C': 5}
PATCHES: Dict[str, int] = {'front': 0, 'back': 1, 'bottom': 2, 'top': 3, 'left': 4, 'right': 5}
NEUMANN = 'neumann'
DIRICHLET = 'dirichlet'


class BoundaryData:
    def __init__(self, field_type: str):
        self.boundary_conditions: List[str] = [""] * len(PATCHES)
        self.values: List[float] = [0] * len(PATCHES)
        self.empty = True
        self.type = field_type

    def add_boundary_data(self, patches: List[str], boundary_condition: str, value: float):
        self.empty = False
        for p in patches:
            self.boundary_conditions[PATCHES[p]] = boundary_condition
            self.values[PATCHES[p]] = value


class Obstacle:
    def __init__(self, name: str):
        self.name = name
        self.geometry: List[float]
        self.boundary: List[BoundaryData] = [BoundaryData(f) for f in FIELD_TYPES.keys()]

    def add_boundary(self, fields: List[str], patches: List[str], boundary_condition: str, value: float):
        for f in fields:
            self.boundary[FIELD_TYPES[f]].add_boundary_data(patches, boundary_condition, value)

    def add_geometry(self, geometry: List[float]):
        self.geometry = geometry
