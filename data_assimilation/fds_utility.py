import os
from typing import Dict, List

import fdsreader
import pandas
import pandas as pd

from ARTSS import Domain


def read_fds_data(input_path: str, input_file_name: str, artss: Domain) -> [Dict[str, Dict[str, any]], pd.DataFrame]:
    full_name = os.path.join(input_path, input_file_name)
    str_devc = read_devc_from_input_file(full_name + '.fds')
    dict_devc = parse_devc(str_devc)
    chid_file_name = get_chid_from_input_file(full_name + '.fds')

    fds_data = read_devc_from_csv_file(os.path.join(input_path, chid_file_name + '_devc.csv'))

    devices: Dict[str, Dict[str, any]] = {}
    devc_temperature: Dict[str, any] = {}
    devc_velocity: Dict[str, any] = {}
    devc_pressure: Dict[str, any] = {}
    devc_visibility: Dict[str, any] = {}
    for d in dict_devc:
        dict_devc[d]['type']: str = 'T'

        ijk = artss.get_ijk_from_xyz(dict_devc[d]['XYZ'][0], dict_devc[d]['XYZ'][2], dict_devc[d]['XYZ'][1])
        dict_devc[d]['index']: int = artss.get_index(*ijk)
        d_lower = d.lower()
        if d_lower.startswith('temperatur'):
            devc_temperature[d] = dict_devc[d]
        elif d_lower.startswith('velo'):
            devc_velocity[d] = dict_devc[d]
        elif d_lower.startswith('vis'):
            devc_visibility[d] = dict_devc[d]
        elif d_lower.startswith('pressure'):
            devc_pressure[d] = dict_devc[d]

    devices['temperature'] = devc_temperature
    devices['velocity'] = devc_velocity
    devices['pressure'] = devc_pressure
    devices['visibility'] = devc_visibility
    return devices, fds_data


def read_devc_from_csv_file(file_name: str) -> pd.DataFrame:
    df = pd.read_csv(file_name, skiprows=1)
    df.index = df.iloc[:, 0]
    df.drop('Time', axis=1, inplace=True)
    return df


def get_chid_from_input_file(input_file: str) -> str:
    with open(input_file, 'r') as inp:
        for line in inp:
            if line.startswith('&HEAD'):
                line = line[len('&HEAD'):].strip()
                tmp_arr = line.split('=')
                index = tmp_arr.index('CHID')
                chid_str = tmp_arr[index + 1].strip()
                index_end = chid_str[1:].index("'")
                chid = chid_str[1:index_end + 1]
                break
    return chid


def parse_devc(str_devc: list, start='&DEVC') -> dict:
    devc = {}
    for line in str_devc:
        line = line[len(start):-1]
        line.strip()
        tmp_arr = line.split(',')
        parameters = []
        for i in tmp_arr:
            if type(i) == str and '=' in i:
                parameters.append(i)
            else:
                parameters[-1] += f',{i}'

        param_dict = {}
        for p in parameters:
            tmp_arr = p.split('=')
            param_dict[tmp_arr[0].strip()] = tmp_arr[1].strip().strip("'").strip('"')
        id = param_dict.pop('ID')
        if 'XYZ' in param_dict.keys():
            param_dict['XYZ'] = [float(i) for i in param_dict['XYZ'].split(',')]
        devc[id] = param_dict
    return devc


def read_devc_from_input_file(fds_file_name: str, start="&DEVC", end="/") -> list:
    """
    Reads FDS input out of a list of string and collects namelist definitions
    of a given type, e.g. "&DEVC ".
    Reads all definitions line by line, without processing.

    :param fds_file_name: List of string, of a single namelist definition.
    :param start: String, indicate the namelist, default "&DEVC ".
    :param end: String, indicate the end of the namelist, default "/".
    :return: list with lines
    inspired by Tristan Hehnen
    """

    namelist_lines = []

    collect_line = False
    fds_input_content = open(fds_file_name, 'r')
    for line in fds_input_content:
        line = line.strip()
        if collect_line is True:
            namelist_lines[-1] += f' {line}'
        elif line.startswith(start):
            namelist_lines.append(line)
            collect_line = True
        collect_line = collect_line and (end not in line)
    fds_input_content.close()
    return namelist_lines


def read_fds_file(data_path: str, artss: Domain) -> pd.DataFrame:
    data = fdsreader.Simulation(data_path)
    devc = data.devices
    return_val: pd.DataFrame = pd.DataFrame([], index=devc['Time'].data)
    return_val.index.names = ['Time']
    return_val.columns.names = ['Position of Thermocouple']
    for key in devc:
        if key.startswith('Thermocouple'):
            if devc[key] is not None:
                print(key, devc[key].xyz, *devc[key].xyz, devc[key].xyz[1])
                # be careful, in FDS the third parameter indicates the height, but in ARTSS it is the second
                i, j, k = artss.get_ijk_from_xyz(devc[key].xyz[0], devc[key].xyz[2], devc[key].xyz[1])
                index = artss.get_index(i, j, k)
                print(f'index: {index} of {i}|{j}|{k}')
                if index in return_val:
                    print(f'index  {index} already exists, take average')
                    tmp = (devc[key].data + return_val[index]) / 2.
                    return_val[index] = tmp
                else:
                    return_val[index] = devc[key].data
            else:
                print(f'{key} is None, skipped')
    return return_val


def get_starting_time(fds_data : pd.DataFrame, threshold: float, keys: List[str]):
    columns = []
    for i in fds_data:
        for key in keys:
            if key in i.lower():
                columns.append(i)
    data_extract = fds_data[columns]
    for index, row in data_extract.iterrows():
        if any(row > threshold):
            return index
