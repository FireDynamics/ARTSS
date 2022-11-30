#!/usr/bin/env python3
import os
import re
import struct
from datetime import datetime
import locale
from typing import List, Dict

import h5py
import numpy as np
from retry import retry
import xml.etree.ElementTree as ET
import xml.dom.minidom as md

from obstacle import Obstacle, FIELD_TYPES, PATCHES


def create_message(t_cur: float, config_file_name: str) -> bin:
    print(f'send message with time step "{t_cur}" and xml "{config_file_name}"')
    package = DAPackage(t_cur, config_file_name)
    return package.pack()


def is_float(x: str) -> bool:
    try:
        float(x)
        return True
    except ValueError:
        return False


def get_date_now() -> str:
    return datetime.now().strftime('%a %b %d %H:%M:%S %Y')


def change_heat_source(source_type: dict, temperature_source: dict, random: dict, changes: dict) -> [dict, dict, dict]:
    source_type_changes = changes['source_type']
    temperature_source_changes = changes['temperature_source']
    random_changes = changes['random']

    new_source_type = source_type.copy()
    new_temperature_source = temperature_source.copy()
    new_random = random.copy()

    for key in source_type_changes:
        new_source_type[key] = source_type_changes[key]

    for key in temperature_source_changes:
        new_temperature_source[key] = temperature_source_changes[key]

    if not source_type['random']:
        new_random = {}
    else:
        for key in random_changes:
            new_random[key] = random_changes[key]

    return new_source_type, new_temperature_source, new_random


class DAPackage:
    def __init__(self, time: float, file_name: str):
        self._time = time
        self._file_name = file_name

    def pack(self) -> bin:
        return struct.pack('di', self._time, len(self._file_name)) + self._file_name.encode('UTF-8')


class DAFile:
    def __init__(self):
        self.xml_root = ET.Element('ARTSS')
        self.data_assimilation: ET.SubElement = ET.SubElement(self.xml_root, 'data_assimilation')

    def write_obstacle_changes(self, list_obstacles: List[Obstacle], obstacle_enabled: bool):
        # create obstacle part
        # <data_assimilation>
        #   <class name="ObstacleChanger" tag="obstacle_changer" />
        # </data_assimilation>
        # <obstacle_changer>
        #   <obstacle name="name1" state="unmodified"/>
        #   <obstacle name="name2" state="modified" >
        #       <geometry ox1 = "0.0273" ox2 = "0.964" oy1 = "0.0078" oy2 = "0.992" oz1 = "-0.492" oz2 = "0.4785" />
        #       <boundary field = "u,v,w" patch = "front,back,left,right,bottom,top" type = "dirichlet" value = "0.0" />
        #       <boundary field = "p" patch = "front,back,left,right,bottom,top" type = "neumann" value = "0.0" />
        #   </obstacle>
        #   <obstacle name="name3" state="new" >
        #       <geometry ox1 = "0.0273" ox2 = "0.964" oy1 = "0.0078" oy2 = "0.992" oz1 = "-0.492" oz2 = "0.4785" />
        #       <boundary field = "u,v,w" patch = "front,back,left,right,bottom,top" type = "dirichlet" value = "0.0" />
        #       <boundary field = "p" patch = "front,back,left,right,bottom,top" type = "neumann" value = "0.0" />
        #   </obstacle>
        # </obstacle_changer>
        tag = 'obstacle_changer'
        ET.SubElement(self.data_assimilation, 'class', {'name': 'ObstacleChanger', 'tag': tag})
        obstacle_root = ET.SubElement(self.xml_root, tag)
        for obstacle in list_obstacles:
            single_obstacle = ET.SubElement(obstacle_root, 'obstacle', name=obstacle.name, state=obstacle.state)
            if obstacle.state != 'unmodified':
                # convert from Dict[str, float] to Dict[str, str]
                geometry = dict(zip(obstacle.geometry.keys(), [str(a) for a in obstacle.geometry.values()]))
                ET.SubElement(single_obstacle, 'geometry', geometry)
                for count, boundary in enumerate(obstacle.boundary):
                    field_type = list(FIELD_TYPES.keys())[count]
                    if not boundary.empty:
                        for p in PATCHES.keys():
                            ET.SubElement(single_obstacle, 'boundary', field=field_type, patch=p,
                                          type=boundary.boundary_conditions[PATCHES[p]],
                                          value=str(boundary.values[PATCHES[p]]))

    def create_obstacle_changes(self, list_obstacles: List[Obstacle], obstacle_enabled: bool):
        # create obstacle part
        # <data_assimilation>
        #   <class name="ObstacleChanger" tag="obstacle_changer" />
        # </data_assimilation>
        # <obstacle_changer>
        #   <obstacle name = "name" state="changed">
        #       <geometry ox1 = "0.0273" ox2 = "0.964" oy1 = "0.0078" oy2 = "0.992" oz1 = "-0.492" oz2 = "0.4785" />
        #       <boundary field = "u,v,w" patch = "front,back,left,right,bottom,top" type = "dirichlet" value = "0.0" />
        #       <boundary field = "p" patch = "front,back,left,right,bottom,top" type = "neumann" value = "0.0" />
        #   </obstacle>
        # </obstacle_changer>
        tag = 'obstacle_changer'
        ET.SubElement(self.data_assimilation, 'class', {'name': 'ObstacleChanger', 'tag': tag})
        obstacle_root = ET.SubElement(self.xml_root, tag)
        for obstacle in list_obstacles:
            single_obstacle = ET.SubElement(obstacle_root, 'obstacle', name=obstacle.name)
            # convert from Dict[str, float] to Dict[str, str]
            geometry = dict(zip(obstacle.geometry.keys(), [str(a) for a in obstacle.geometry.values()]))
            ET.SubElement(single_obstacle, 'geometry', geometry)
            for count, boundary in enumerate(obstacle.boundary):
                field_type = list(FIELD_TYPES.keys())[count]
                if not boundary.empty:
                    for p in PATCHES.keys():
                        ET.SubElement(single_obstacle, 'boundary', field=field_type, patch=p,
                                      type=boundary.boundary_conditions[PATCHES[p]],
                                      value=str(boundary.values[PATCHES[p]]))

    def create_temperature_source_changes(self, source_type: dict,
                                          temperature_source: dict,
                                          random: dict):
        # create temperature source part
        # <data_assimilation>
        #   <class name="TemperatureSourceChanger" tag="temperature_source" />
        # </data_assimilation>
        # <temperature_source>
        #   <source type="ExplicitEuler" dir="y" temp_fct="Gauss" dissipation="No" random="Yes">
        #     <HRR> 25000. </HRR> <!-- Total heat release rate( in kW) -->
        #     <cp> 1.023415823 </cp> <!-- specific heat capacity( in kJ / kgK)-->
        #     <x0> 40. </x0>
        #     <y0> -3. </y0>
        #     <z0> 0. </z0>
        #     <sigma_x> 1.0 </sigma_x>
        #     <sigma_y> 1.5 </sigma_y>
        #     <sigma_z> 1.0 </sigma_z>
        #     <tau> 5. </tau>
        #     <random absolute="Yes" custom_seed="Yes" custom_steps="Yes">
        #       <seed> 0 </seed>
        #       <step_size> 0.1 </step_size>
        #       <range> 1 </range>
        #     </random>
        #   </source>
        # </temperature_source>
        tag = 'temperature_source'
        ET.SubElement(self.data_assimilation, 'class', {'name': 'TemperatureSourceChanger', 'tag': tag})
        source_root = ET.SubElement(self.xml_root, tag)
        source = ET.SubElement(source_root, 'source', type=source_type['type'], dir=source_type['dir'],
                               temp_fct=source_type['temp_fct'],
                               dissipation='Yes' if source_type['dissipation'] == 'Yes' else 'No',
                               random='Yes' if source_type['random'] == 'Yes' else 'No')
        for key in temperature_source:
            ET.SubElement(source, key).text = str(temperature_source[key])

        if source_type['random'] == 'Yes':
            random_tree = ET.SubElement(source, 'random', absolute='Yes' if random['absolute'] == 'Yes' else 'No',
                                        custom_seed='Yes' if random['custom_seed'] == 'Yes' else 'No',
                                        custom_steps='Yes' if random['custom_steps'] == 'Yes' else 'No')
            ET.SubElement(random_tree, 'range').text = str(random['range'])
            if random['custom_steps'] == 'Yes':
                ET.SubElement(random_tree, 'step_size').text = str(random['step_size'])
            if random['custom_seed'] == 'Yes':
                ET.SubElement(random_tree, 'seed').text = str(random['seed'])

    def create_field_changes(self, fields: Dict[str, bool], field_file_name=''):
        # create field changes part
        # <data_assimilation>
        #   <class name="FieldChanger" tag="field_changes" />
        # </data_assimilation>
        # <field_changes>
        #   <fields_changed u="No" v="No" w="No" p="No" T="Yes" concentration="No" filename="field_file_name"/>
        # </field_changes>
        tag: str = 'field_changes'
        keys = ['u', 'v', 'w', 'p', 'T', 'C']
        ET.SubElement(self.data_assimilation, 'class', {'name': 'FieldChanger', 'tag': tag})
        field_values = {}
        for key in keys:
            if key in fields.keys():
                field_values[key] = 'Yes' if fields[key] else 'No'
            else:
                field_values[key] = 'No'
        subsection = ET.SubElement(self.xml_root, tag)
        if 'Yes' in field_values.values():
            ET.SubElement(subsection, 'fields_changed', u=field_values['u'], v=field_values['v'],
                          w=field_values['w'], p=field_values['p'],
                          T=field_values['T'], concentration=field_values['C'], filename=field_file_name)
        else:
            ET.SubElement(subsection, 'fields_changed', u=field_values['u'], v=field_values['v'],
                          w=field_values['w'], p=field_values['p'],
                          T=field_values['T'], concentration=field_values['C'])

    # write xml to file
    def write_xml(self, config_file_name: str, pretty_print=False):
        tree = ET.ElementTree(self.xml_root)
        tree.write(config_file_name)
        print(config_file_name)
        if pretty_print:
            dom = md.parse(config_file_name)
            pretty_format = dom.toprettyxml()
            f = open(config_file_name, 'w')
            f.write(pretty_format)
            f.close()


class FieldReader:
    def __init__(self, time_step, path: str = '.'):
        self.t: float = time_step
        self.dt: float = 0
        self.xml_file_name: str = ''
        self.grid_resolution: Dict[str, int] = {}
        self.fields: List[str] = []
        self.date: datetime
        self.header: str
        self.path: str = path
        self.read_header()

    @retry(delay=1, tries=6)
    def read_header(self):
        file = os.path.join(self.path, '.vis', f'{self.t:.5e}')
        with h5py.File(file, 'r') as inp:
            metadata = inp['metadata']
            self.dt = metadata['dt'][()]
            domain = list(metadata['domain'][:])
            self.grid_resolution = {'Nx': domain[0], 'Ny': domain[1], 'Nz': domain[2]}
            self.fields = list(metadata['fields'].asstr()[:])
            loc = locale.getlocale()
            locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
            self.date = datetime.strptime(metadata['date'].asstr()[()][0].strip(), '%a %b %d %H:%M:%S %Y')
            locale.setlocale(locale.LC_ALL, loc)
            self.xml_file_name = metadata['xml'][()][0]

    def print_header(self):
        print(f'dt: {self.dt}')
        print(f'grid resolution: {self.grid_resolution}')
        print(f'fields: {self.fields}')
        print(f'date: {self.date}')
        print(f'xml file name: {self.xml_file_name}')

    @staticmethod
    @retry(delay=1, tries=6)
    def get_t_current(path: str = '.') -> float:
        with open(os.path.join(path, '.vis/meta'), 'r') as inp:
            t = float([x for x in inp.readlines() if x.startswith('t:')][0][2:])
        return t

    @staticmethod
    @retry(delay=1, tries=6)
    def get_all_time_steps(path: str = '.') -> List[float]:
        pattern = re.compile('[0-9]+\.[0-9]{5}e[+|-][0-9]+')
        f = os.listdir(os.path.join(path, '.vis'))
        files = []
        for p in f:
            if pattern.match(p):
                files.append(p)
        files = [float(x) for x in files]
        files.sort()
        return files

    @staticmethod
    def get_xml_file_name(path: str = '.') -> str:
        fpath = os.path.join(path, '.vis/meta')
        with open(fpath, 'r') as inp:
            xml_file_name = [x for x in inp.readlines() if x.startswith('xml_name:')][0][len('xml_name:'):]
        return xml_file_name.strip()

    def get_grid_resolution(self) -> Dict[str, int]:
        return self.grid_resolution

    def get_fields(self) -> list:
        return self.fields

    def read_field_data(self) -> Dict[str, np.ndarray]:
        t_cur = self.get_t_current(path=self.path)
        if self.t > t_cur:
            print(f'cannot read time step {self.t} as the current time step is {t_cur}')
            return {}
        else:
            fields = {}
            print(f'read time step {self.t}')
            with h5py.File(os.path.join(self.path, '.vis', f'{self.t:.5e}'), 'r') as inp:
                inp = inp[f'/{self.t:.5e}']
                for i in inp:
                    fields[i] = np.array(inp[i][0])
                    print(f'read field: {i} {fields[i].shape}, {sum(fields[i])}')

            return fields

    @staticmethod
    def write_field_data_keys(file_name: str, data: Dict[str, np.ndarray], field_keys: List[str]):
        with h5py.File(file_name, 'w') as out:
            for key in field_keys:
                if key not in data.keys():
                    continue

                field = data[key]
                dset = out.create_dataset(key, (len(field),), dtype='d')
                dset[:] = field

    def write_field_data(self, file_name: str, data: Dict[str, np.ndarray]):
        FieldReader.write_field_data_keys(file_name=file_name, data=data, field_keys=self.fields)
