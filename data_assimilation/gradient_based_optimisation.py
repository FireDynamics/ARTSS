import os
import time
import fdsreader
import pandas as pd

from ARTSS import Domain, XML
from data_assimilation import FieldReader


def start(fds_data_path: str, fds_input_file_name: str, artss_data_path: str):
    xml = XML(FieldReader.get_xml_file_name(artss_data_path), path=artss_data_path)
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()

    devc_info_temperature, devc_info_thermocouple, fds_data = read_fds_data(fds_data_path, fds_input_file_name, domain)
    print(devc_info_temperature)

    sensor_times = fds_data.index

    nabla: dict = {}
    for t_sensor in sensor_times[1:]:
        t_cur = FieldReader.get_t_current(path=artss_data_path)
        print(t_cur)
        while t_cur < t_sensor:
            time.sleep(10)
        t_artss = get_time_step_artss(t_sensor, os.path.join(artss_data_path, '.vis'))
        field_reader = FieldReader(t_artss, path=artss_data_path)
        diff = comparison_sensor_simulation_data(devc_info_thermocouple, fds_data, field_reader, t_artss, t_sensor)
        if not nabla:
            # first change
            pass
        else:
            # create nabla
            calc_nabla()

def calc_nabla(old_diff, new_diff, old_params, new_params) -> list:
    # calc nabla
    pass

def get_time_step_artss(t_sensor: float, artss_data_path: str) -> float:
    files = os.listdir(artss_data_path)
    files.remove('meta')
    files = [float(x) for x in files]
    files.sort()

    estimated_index = int(t_sensor * len(files) / files[-1])

    t_artss: float = -1
    if abs(files[estimated_index] - t_sensor) < 1e-10:
        t_artss = files[estimated_index]
    elif files[estimated_index] < t_sensor:
        # ascending
        while files[estimated_index] < len(files):
            if files[estimated_index] >= t_sensor:
                t_artss = files[estimated_index]
                break
            else:
                estimated_index += 1
    else:
        # descending
        while files[estimated_index] >= 0:
            if files[estimated_index] <= t_sensor:
                t_artss = files[estimated_index + 1]
                break
            else:
                estimated_index -= 1
    return t_artss


def read_fds_data(input_path: str, input_file_name: str, artss: Domain) -> [dict, dict, pd.DataFrame]:
    full_name = os.path.join(input_path, input_file_name)
    str_devc = read_devc_from_input_file(full_name + '.fds')
    dict_devc = parse_devc(str_devc)
    chid_file_name = get_chid_from_input_file(full_name + '.fds')

    fds_data = read_devc_from_csv_file(os.path.join(input_path, chid_file_name + '_devc.csv'))
    # fds_data = read_fds_file(input_path, artss)

    devc_temperature = {}
    devc_thermocouple = {}
    for d in dict_devc:
        dict_devc[d]['type']: str = 'T'
        dict_devc[d]['index']: int = artss.get_index(*dict_devc[d]['XYZ'])
        if d.startswith('Temperatur'):
            devc_temperature[d] = dict_devc[d]
        elif d.startswith('Thermocouple'):
            devc_thermocouple[d] = dict_devc[d]

    return devc_temperature, devc_thermocouple, fds_data


def read_devc_from_csv_file(file_name: str) -> pd.DataFrame:
    df = pd.read_csv(file_name, skiprows=1)
    df.index = df.iloc[:, 0]
    df.drop('Time', axis=1, inplace=True)
    return df


def get_chid_from_input_file(input: str) -> str:
    with open(input, 'r') as inp:
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
                i, j, k = artss.get_ijk_from_xyz(devc[key].xyz[0], devc[key].xyz[1], devc[key].xyz[2])
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


def gradient_based_optimisation(sensor_data: pd.DataFrame, domain: Domain, field_reader: FieldReader, t_cur: float):
    if not comparison_sensor_simulation_data(sensor_data, field_reader, t_cur):
        # initiate gradient based optimisation, simulation data differs to much from simulation data
        pass
    return False


def comparison_sensor_simulation_data(devc_info: dict, sensor_data: pd.DataFrame, field_reader: FieldReader, t_artss: float, t_sensor: float) -> dict:
    diff: dict = {'T': [], 'C': []}
    fields_sim = field_reader.read_field_data(t_artss)
    for key in devc_info:
        # sensor data
        type_sensor: str = devc_info[key]['type']
        index_sensor: int = devc_info[key]['index']
        value_sensor: float = sensor_data[key][t_sensor]

        value_sim = fields_sim[type_sensor][index_sensor]
        diff[type_sensor].append(value_sim - value_sensor)

    return diff
