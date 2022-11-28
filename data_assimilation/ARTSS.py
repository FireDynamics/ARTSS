#!/usr/bin/env python3
import os
import xml.etree.ElementTree as ET
from typing import List, Dict, Union, Set

from obstacle import Obstacle, PATCHES, STATE


class XML:
    def __init__(self, filename: str, path: str = '.'):
        self.filename: str = filename
        self.xml_tree: ET = ET.parse(os.path.join(path, self.filename))
        self.has_obstacles: bool = False
        self.obstacles: List[Obstacle] = []
        self.computational_domain: bool = False
        self.domain: Dict[str, float] = {}
        self.temperature_source: Dict[str, Dict[str, str]] = {}

    def read_xml(self):
        root = self.xml_tree.getroot()
        domain_param = root.find('domain_parameters')
        self.computational_domain = domain_param.attrib == "Yes"
        for child in domain_param:
            self.domain[child.tag] = float(child.text)

        obstacles = root.find('obstacles')
        if obstacles.attrib['enabled'] == 'Yes':
            self.has_obstacles = True
            for child in obstacles:
                state = 'XML'
                if 'state' in child.attrib:
                    state = STATE[child.attrib['state']]
                obstacle = Obstacle(name=child.attrib['name'], state=state)
                for content in child:
                    if content.tag == 'boundary':
                        obstacle.add_boundary_line(content.attrib)

                geometry = child.find('geometry')
                obstacle.add_geometry_line(geometry.attrib)
                self.obstacles.append(obstacle)

    def get_output_time_step(self) -> float:
        root = self.xml_tree.getroot()
        da_param = root.find('data_assimilation')
        dict_da = {}
        for key in da_param.attrib:
            dict_da[key] = (da_param.attrib[key])
        if 'write_output' in dict_da.keys():
            return float(dict_da['write_output'])
        return 1

    def get_temperature_source(self) -> dict:
        if self.temperature_source is None:
            # check if temperature source is even there ? or crash
            root = self.xml_tree.getroot()
            source_tree = root.find('solver').find('temperature').find('source')
            source = {}
            for key in source_tree.attrib:
                source[key] = source_tree.attrib[key].strip()

            source_params = {}
            for elem in source_tree:
                source_params[elem.tag] = elem.text.strip()
            random_params = {}
            if source['random'] == 'Yes':
                del source_params['random']  # remove empty random element
                random = source_tree.find('random')
                for key in random.attrib:
                    random_params[key] = random.attrib[key].strip()
                for elem in random:
                    random_params[elem.tag] = elem.text.strip()
            self.temperature_source = {'source_type': source, 'temperature_source': source_params,
                                       'random': random_params}
        return self.temperature_source

    def get_dt(self) -> float:
        pp = self.xml_tree.getroot().find('physical_parameters')
        for elem in pp:
            if elem.tag == 'dt':
                return float(elem.text)


class Domain:
    def __init__(self, domain_param: Dict[str, float], enable_computational_domain: bool, obstacles: List[Obstacle]):
        self.domain_param: Dict[str, Union[int, float]] = {}
        self.obstacles: List[Obstacle] = obstacles
        self.domain_boundary_cells: Dict[str, Set[int]] = {}
        self.domain_inner_cells: Dict[str, Set[int]] = {}
        self.obstacle_cells: Dict[str, Dict[str, Set[int]]] = {}
        self.computational_domain(domain_param, enable_computational_domain)
        self.calculate_domain()
        self.calculate_indices()

    def computational_domain(self, domain_param: Dict, enable_computational_domain: bool):
        self.domain_param = dict(domain_param)
        if enable_computational_domain:
            self.domain_param['x1'] = domain_param['X1']
            self.domain_param['y1'] = domain_param['Y1']
            self.domain_param['z1'] = domain_param['Z1']
            self.domain_param['x2'] = domain_param['X2']
            self.domain_param['y2'] = domain_param['Y2']
            self.domain_param['z2'] = domain_param['Z2']

    def calculate_domain(self):
        # from Domain.cpp/h
        self.domain_param['Lx'] = abs(self.domain_param['X2'] - self.domain_param['X1'])
        self.domain_param['Ly'] = abs(self.domain_param['Y2'] - self.domain_param['Y1'])
        self.domain_param['Lz'] = abs(self.domain_param['Z2'] - self.domain_param['Z1'])

        self.domain_param['lx'] = abs(self.domain_param['x2'] - self.domain_param['x1'])
        self.domain_param['ly'] = abs(self.domain_param['y2'] - self.domain_param['y1'])
        self.domain_param['lz'] = abs(self.domain_param['z2'] - self.domain_param['z1'])

        self.domain_param['dx'] = self.domain_param['lx'] / float(self.domain_param['nx'])
        self.domain_param['dy'] = self.domain_param['ly'] / float(self.domain_param['ny'])
        self.domain_param['dz'] = self.domain_param['lz'] / float(self.domain_param['nz'])

        self.domain_param['nx'] = int(self.domain_param['nx'])
        self.domain_param['ny'] = int(self.domain_param['ny'])
        self.domain_param['nz'] = int(self.domain_param['nz'])
        self.domain_param['Nx'] = int(round(self.domain_param['Lx'] / self.domain_param['dx'] + 2))
        self.domain_param['Ny'] = int(round(self.domain_param['Ly'] / self.domain_param['dy'] + 2))
        self.domain_param['Nz'] = int(round(self.domain_param['Lz'] / self.domain_param['dz'] + 2))

        self.domain_param['start x index PD'] = 1
        self.domain_param['start y index PD'] = 1
        self.domain_param['start z index PD'] = 1
        self.domain_param['end x index PD'] = self.domain_param['Nx'] - 2
        self.domain_param['end y index PD'] = self.domain_param['Ny'] - 2
        self.domain_param['end z index PD'] = self.domain_param['Nz'] - 2

        self.domain_param['start x index CD'] = int(
            round((self.domain_param['x1'] - self.domain_param['X1']) / self.domain_param['dx'])) + 1
        self.domain_param['start y index CD'] = int(
            round((self.domain_param['y1'] - self.domain_param['Y1']) / self.domain_param['dy'])) + 1
        self.domain_param['start z index CD'] = int(
            round((self.domain_param['z1'] - self.domain_param['Z1']) / self.domain_param['dz'])) + 1
        self.domain_param['end x index CD'] = self.domain_param['start x index CD'] + self.domain_param['nx'] - 1
        self.domain_param['end y index CD'] = self.domain_param['start y index CD'] + self.domain_param['ny'] - 1
        self.domain_param['end z index CD'] = self.domain_param['start z index CD'] + self.domain_param['nz'] - 1

    def match_grid(self, obstacle_coordinate: float, direction: str) -> int:
        # from function get_matching_index in Obstacle.cpp for obstacles only
        return int(round((-self.domain_param[direction.upper() + '1'] + obstacle_coordinate) / self.domain_param[
            'd' + direction.lower()]))

    def get_index(self, i: int, j: int, k: int) -> int:
        return int(i + self.domain_param['Nx'] * j + self.domain_param['Nx'] * self.domain_param['Ny'] * k)

    def calculate_indices(self):
        for o in self.obstacles:
            self.calculate_obstacle_cells(o)
        self.calculate_domain_boundaries()
        self.calculate_domain_inner()

    def calculate_obstacle_cells(self, o: Obstacle):
        # indices of inner cells
        inner_cells: List[int] = []
        index_x1: int = self.match_grid(o.geometry['ox1'], 'x')
        index_x2: int = self.match_grid(o.geometry['ox2'], 'x')
        index_y1: int = self.match_grid(o.geometry['oy1'], 'y')
        index_y2: int = self.match_grid(o.geometry['oy2'], 'y')
        index_z1: int = self.match_grid(o.geometry['oz1'], 'z')
        index_z2: int = self.match_grid(o.geometry['oz2'], 'z')
        for i in range(index_x1 + 1, index_x2):
            for j in range(index_y1 + 1, index_y2):
                for k in range(index_z1 + 1, index_z2):
                    inner_cells.append(self.calculate_index(i, j, k))
        self.obstacle_cells[o.name]['inner cells'] = set(inner_cells)

        # indices of obstacle boundary cells
        front = []
        back = []
        for i in range(index_x1, index_x2 + 1):
            for j in range(index_y1, index_y2 + 1):
                front.append(self.calculate_index(i, j, index_z1))
                back.append(self.calculate_index(i, j, index_z2))
        self.obstacle_cells[o.name]['front'] = set(front)
        self.obstacle_cells[o.name]['back'] = set(back)

        bottom = []
        top = []
        for i in range(index_x1, index_x2 + 1):
            for k in range(index_z1, index_z2 + 1):
                bottom.append(self.calculate_index(i, index_y1, k))
                top.append(self.calculate_index(i, index_y2, k))
        self.obstacle_cells[o.name]['bottom'] = set(bottom)
        self.obstacle_cells[o.name]['top'] = set(top)

        left = []
        right = []
        for j in range(index_y1, index_y2 + 1):
            for k in range(index_z1, index_z2 + 1):
                left.append(self.calculate_index(index_x1, j, k))
                right.append(self.calculate_index(index_x2, j, k))
        self.obstacle_cells[o.name]['left'] = set(left)
        self.obstacle_cells[o.name]['right'] = set(right)

    def calculate_domain_boundaries(self):
        front = []
        back = []
        for i in range(self.domain_param['start x index CD'] - 1, self.domain_param['end x index CD'] + 2):
            for j in range(self.domain_param['start y index CD'] - 1, self.domain_param['end y index CD'] + 2):
                front.append(self.calculate_index(i, j, self.domain_param['start z index CD'] - 1))
                back.append(self.calculate_index(i, j, self.domain_param['end z index CD'] + 1))
        self.domain_boundary_cells['front'] = set(front)
        self.domain_boundary_cells['back'] = set(back)

        bottom = []
        top = []
        for i in range(self.domain_param['start x index CD'] - 1, self.domain_param['end x index CD'] + 2):
            for k in range(self.domain_param['start z index CD'] - 1, self.domain_param['end z index CD'] + 2):
                bottom.append(self.calculate_index(i, self.domain_param['start y index CD'] - 1, k))
                top.append(self.calculate_index(i, self.domain_param['end y index CD'] + 1, k))
        self.domain_boundary_cells['bottom'] = set(bottom)
        self.domain_boundary_cells['top'] = set(top)

        left = []
        right = []
        for j in range(self.domain_param['start y index CD'] - 1, self.domain_param['end y index CD'] + 2):
            for k in range(self.domain_param['start z index CD'] - 1, self.domain_param['end z index CD'] + 2):
                left.append(self.calculate_index(self.domain_param['start x index CD'] - 1, j, k))
                right.append(self.calculate_index(self.domain_param['end x index CD'] + 1, j, k))
        self.domain_boundary_cells['left'] = set(left)
        self.domain_boundary_cells['right'] = set(right)

    def calculate_domain_inner(self):
        inner = []
        for i in range(self.domain_param['start x index CD'], self.domain_param['end x index CD']):
            for j in range(self.domain_param['start y index CD'], self.domain_param['end y index CD']):
                for k in range(self.domain_param['start z index CD'], self.domain_param['end z index CD']):
                    inner.append(self.calculate_index(i, j, k))
        inner = set(inner)
        for oc in self.obstacle_cells:
            for type_of_cell in self.obstacle_cells[oc]:
                inner = inner - self.obstacle_cells[oc][type_of_cell]
        self.domain_inner_cells = inner

    def calculate_index(self, i, j, k):
        return i + j * self.domain_param['Nx'] + k * self.domain_param['Nx'] * self.domain_param['Ny']

    def print_info(self):
        # alike to info output in logger
        print(
            f"Domain size inner cells: {self.domain_param['nx']} {self.domain_param['ny']} {self.domain_param['nz']}\n"
            f"step size (x|y|z): ({self.domain_param['dx']}|{self.domain_param['dy']}|{self.domain_param['dz']})")
        for o in self.obstacles:
            print(f"-- Obstacle {o.name}\n"
                  f"   strides (x y z): {o.stride['stride x']} {o.stride['stride y']} {o.stride['stride z']}\n"
                  f"   size of slices (Front|Back Bottom|Top Left|Right): "
                  f"{len(self.obstacle_cells[o.name]['front'])}|"
                  f"{len(self.obstacle_cells[o.name]['back'])} "
                  f"{len(self.obstacle_cells[o.name]['bottom'])}|"
                  f"{len(self.obstacle_cells[o.name]['top'])} "
                  f"{len(self.obstacle_cells[o.name]['left'])}|"
                  f"{len(self.obstacle_cells[o.name]['right'])}\n"
                  f"   size of Obstacle: {o.stride['stride x'] * o.stride['stride y'] * o.stride['stride z']}\n"
                  f"   coords (x y z): "
                  f"({o.geometry['ox1']}|{o.geometry['ox2']}) "
                  f"({o.geometry['oy1']}|{o.geometry['oy2']}) "
                  f"({o.geometry['oz1']}|{o.geometry['oz2']})")

    def print_debug(self):
        # alike to debug output in logger
        print(f"Nx: {self.domain_param['Nx']}, Ny: {self.domain_param['Ny']}, Nz: {self.domain_param['Nz']}")
        print(f"nx: {self.domain_param['nx']}, ny: {self.domain_param['ny']}, nz: {self.domain_param['nz']}")
        print(f"X: ({self.domain_param['X1']}|{self.domain_param['X2']}) "
              f"Y: ({self.domain_param['Y1']}|{self.domain_param['Y2']}) "
              f"Z: ({self.domain_param['Z1']}|{self.domain_param['Z2']})")
        print(f"Lx: {self.domain_param['Lx']}, Ly: {self.domain_param['Ly']}, Lz: {self.domain_param['Lz']}")
        print(f"lx: {self.domain_param['lx']}, ly: {self.domain_param['ly']}, lz: {self.domain_param['lz']}")
        print(f"dx: {self.domain_param['dx']}, dy: {self.domain_param['dy']}, dz: {self.domain_param['dz']}")
        print(f"index X CD: ({self.domain_param['start x index CD']}|{self.domain_param['end x index CD']}) "
              f"index Y CD: ({self.domain_param['start y index CD']}|{self.domain_param['end y index CD']}) "
              f"index Z CD: ({self.domain_param['start z index CD']}|{self.domain_param['end z index CD']})")
        print(f"index X PD: ({self.domain_param['start x index PD']}|{self.domain_param['end x index PD']}) "
              f"index Y PD: ({self.domain_param['start y index PD']}|{self.domain_param['end y index PD']}) "
              f"index Z PD: ({self.domain_param['start z index PD']}|{self.domain_param['end z index PD']})")

    def get_ijk_from_xyz(self, coord_x, coord_y, coord_z):
        return self.get_i_from_x(coord_x), self.get_j_from_y(coord_y), self.get_k_from_z(coord_z)

    def get_i_from_x(self, coord_x):
        return int((coord_x - self.domain_param['X1']) / self.domain_param['dx'])

    def get_j_from_y(self, coord_y):
        return int((coord_y - self.domain_param['Y1']) / self.domain_param['dy'])

    def get_k_from_z(self, coord_z):
        return int((coord_z - self.domain_param['Z1']) / self.domain_param['dz'])

    def get_i(self, index, k, j):
        return int(float(index) - k * self.domain_param['Nx'] * self.domain_param['Ny'] - j * self.domain_param['Nx'])

    def get_j(self, index, k):
        return int((float(index) - k * self.domain_param['Nx'] * self.domain_param['Ny']) / self.domain_param['Nx'])

    def get_k(self, index):
        return int(float(index) / (self.domain_param['Nx'] * self.domain_param['Ny']))

    def get_type(self, index):
        k = self.get_k(index)
        j = self.get_j(index, k)
        i = self.get_i(index, k, j)
        matches = []
        for p in PATCHES.keys():
            if index in self.domain_boundary_cells[p]:
                matches.append(f'domain boundary cell ({p})')
        for oc in self.obstacle_cells:
            for p in PATCHES.keys():
                if index in self.obstacle_cells[oc][p]:
                    matches.append(f'obstacle boundary cell ({oc}/{p})')
            if index in self.obstacle_cells['inner cells']:
                matches.append(f'obstacle inner cell ({oc})')
        if len(matches) == 0:
            return "cell type unknown"
        return f'cell {index} ({i}|{j}|{k}) was found in {matches}'
