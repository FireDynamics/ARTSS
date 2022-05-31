#!/usr/bin/env python3
import os
from datetime import datetime

import h5py
import numpy as np


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


class FieldReader:
    def __init__(self, time_step, path: str = '.'):
        self.t: float = time_step
        self.dt: float
        self.xml_file_name: str
        self.grid_resolution: dict
        self.fields: list
        self.date: datetime
        self.header: str
        self.path = path
        self.read_header()

    def read_header(self):
        with h5py.File(os.path.join(self.path, '.vis', f'{self.t:.5e}'), 'r') as inp:
            metadata = inp['metadata']
            self.dt = metadata['dt'][()]
            domain = list(metadata['domain'][:])
            self.grid_resolution = {'Nx': domain[0], 'Ny': domain[1], 'Nz': domain[2]}
            self.fields = list(metadata['fields'].asstr()[:])
            self.date = datetime.strptime(metadata['date'].asstr()[()][0].strip(), '%a %b %d %H:%M:%S %Y')
            self.xml_file_name = metadata['xml'][()][0]

    def print_header(self):
        print(f'dt: {self.dt}')
        print(f'grid resolution: {self.grid_resolution}')
        print(f'fields: {self.fields}')
        print(f'date: {self.date}')
        print(f'xml file name: {self.xml_file_name}')

    @staticmethod
    def get_t_current(path: str = '.') -> float:
        with open(os.path.join(path, '.vis/meta'), 'r') as inp:
            t = float([x for x in inp.readlines() if x.startswith('t:')][0][2:])
        return t

    @staticmethod
    def get_xml_file_name(path: str = '.') -> str:
        fpath = os.path.join(path, '.vis/meta')
        with open(fpath, 'r') as inp:
            xml_file_name = [x for x in inp.readlines() if x.startswith('xml_name:')][0][len('xml_name:'):]
        return xml_file_name.strip()

    def get_grid_resolution(self) -> dict:
        return self.grid_resolution

    def get_fields(self) -> list:
        return self.fields

    def read_field_data(self, time_step: float) -> dict:
        t_cur = self.get_t_current(path=self.path)
        if time_step > t_cur:
            print(f'cannot read time step {time_step} as the current time step is {t_cur}')
            return {}
        else:
            fields = {}
            print(f'read time step {time_step}')
            with h5py.File(os.path.join(self.path, '.vis', f'{time_step:.5e}'), 'r') as inp:
                inp = inp[f'/{time_step:.5e}']
                for i in inp:
                    fields[i] = np.array(inp[i][0])
                    print(f'read field: {i} {fields[i].shape}, {sum(fields[i])}')

            return fields

    @staticmethod
    def write_field_data_keys(file_name: str, data: dict, field_keys: list):
        with h5py.File(file_name, 'w') as out:
            for key in field_keys:
                if key not in data.keys():
                    continue

                field = data[key]
                out.create_dataset(key, (len(field),), dtype='d')

    def write_field_data(self, file_name: str, data: dict):
        FieldReader.write_field_data_keys(file_name=file_name, data=data, field_keys=self.fields)
