#!/usr/bin/env python3

import os
import struct
import tempfile
from datetime import datetime
from typing import Type

import h5py
import wsgiref.validate
import numpy as np
import pandas as pd

from ARTSS import Domain


def is_float(x: str) -> bool:
    try:
        float(x)
        return True
    except ValueError:
        return False


def is_float(x: str) -> bool:
    try:
        float(x)
        return True
    except ValueError:
        return False


def get_date_now() -> str:
    return datetime.now().strftime('%a %b %d %H:%M:%S %Y')



class FieldReader:
    def __init__(self, time_step):
        self.t: float = time_step
        self.dt: float
        self.xml_file_name: str
        self.grid_resolution: dict
        self.fields: list
        self.date: datetime
        self.header: str
        self.read_header()

    def read_header(self):
        with h5py.File(f'./.vis/{self.t:.5e}', 'r') as inp:
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
    def get_t_current() -> float:
        with open('./.vis/meta', 'r') as inp:
            t = float([x for x in inp.readlines() if x.startswith('t:')][0][2:])
        return t

    @staticmethod
    def get_xml_file_name() -> str:
        with open('./.vis/meta', 'r') as inp:
            xml_file_name = [x for x in inp.readlines() if x.startswith('xml_name:')][0][len('xml_name:'):]
        return xml_file_name.strip()

    def get_grid_resolution(self) -> dict:
        return self.grid_resolution

    def get_fields(self) -> list:
        return self.fields

    def read_field_data(self, time_step: float) -> dict:
        t_cur = self.get_t_current()
        if time_step > t_cur:
            print(f'cannot read time step {time_step} as the current time step is {t_cur}')
            return {}
        else:
            fields = {}
            print(f'read time step {time_step}')
            with h5py.File(f'./.vis/{time_step:.5e}', 'r') as inp:
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
                out.create_dataset(key, (len(field), ), dtype='d')

    def write_field_data(self, file_name: str, data: dict):
        FieldReader.write_field_data_keys(file_name=file_name, data=data, field_keys=self.fields)



def gradient_based_optimisation(sensor_data: pd.DataFrame, domain: Type[Domain], field_reader: Type[FieldReader],
                                t_cur: float):
    if not comparison_sensor_simulation_data(sensor_data, field_reader, t_cur):
        # initiate gradient based optimisation, simulation data differs to much from simulation data
        pass
    return False


def comparison_sensor_simulation_data(sensor_data: pd.DataFrame, field_reader: Type[FieldReader], t_cur: float):
    accurate = True
    nabla: dict[str, float] = {}

    def compare(val_sen: float, val_sim: float, type_sen: str):
        if type_sen == 'T':
            diff = 10
        elif type_sen == 'C':
            diff = 11
        else:
            diff = 3
        return abs(val_sen - val_sim) < diff

    fields_sim = field_reader.read_field_data(t_cur)
    for key in sensor_data:
        # sensor data
        type_sensor = sensor_data[key]['type']
        index_sensor = sensor_data[key]['index']
        value_sensor = sensor_data[key]['value']

        value_sim = fields_sim[type][index_sensor]

        if not compare(value_sensor, value_sim, type_sensor):
            accurate = False

    return accurate, nabla
