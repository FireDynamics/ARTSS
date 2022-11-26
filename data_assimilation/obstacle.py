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
    def __init__(self, name: str, state: str = 'unmodified'):
        self.name = name
        self.geometry: Dict[str, float] = {}
        self.boundary: List[BoundaryData] = [BoundaryData(f) for f in FIELD_TYPES.keys()]
        self.state = state
        self.stride: Dict[str, float] = {}

    def add_boundary(self, fields: List[str], patches: List[str], boundary_condition: str, value: float):
        for f in fields:
            self.boundary[FIELD_TYPES[f]].add_boundary_data(patches, boundary_condition, value)

    def add_boundary_line(self, boundary: Dict[str, str]):
        fields = boundary['field'].split(",")
        patches = boundary['patch'].split(",")
        self.add_boundary(fields=fields, patches=patches, boundary_condition=boundary['type'], value=float(boundary['value']))

    def add_geometry(self, geometry: List[float]):
        keys = ['ox1', 'ox2', 'oy1', 'oy2', 'oz1', 'oz2']
        self.geometry = dict(zip(keys, geometry))

    def add_geometry_line(self, geometry: Dict[str, str]):
        self.geometry = dict(zip(geometry.keys(), [float(a) for a in geometry.values()]))

    def calculate_details(self):
        self.stride['stride x'] = self.geometry['ox2'] - self.geometry['ox1'] + 1
        self.stride['stride y'] = self.geometry['oy2'] - self.geometry['oy1'] + 1
        self.stride['stride z'] = self.geometry['oz2'] - self.geometry['oz1'] + 1
