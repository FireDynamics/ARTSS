import os
import time
import fdsreader
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from copy import copy
from pprint import pprint

import TCP_client
import data_assimilation
from ARTSS import Domain, XML, DAFile
from data_assimilation import FieldReader
from main import create_message


def start(fds_data_path: str, fds_input_file_name: str, artss_data_path: str):
    client = TCP_client.TCPClient()
    client.set_server_address('172.17.0.2', 7777)
    client.connect()
    cwd = os.getcwd()

    xml = XML(FieldReader.get_xml_file_name(artss_data_path), path=artss_data_path)
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    devc_info_temperature, devc_info_thermocouple, fds_data = read_fds_data(fds_data_path, fds_input_file_name, domain)
    print(devc_info_temperature)

    sensor_times = fds_data.index

    source_type, temperature_source, random = xml.get_temperature_source()

    keys = ['HRR', 'x0'] #, 'y0', 'z0']
    delta = {
        'HRR': float(temperature_source['HRR']) * 0.1,
        'x0': domain.domain_param['dx'] * 5,
        'y0': domain.domain_param['dy'] * 1,
        'z0': domain.domain_param['dz'] * 1
    }

    cur = {
        'HRR': float(temperature_source['HRR']),
        'x0': float(temperature_source['x0']),
        'y0': float(temperature_source['y0']),
        'z0': float(temperature_source['z0'])
    }

    def f(t_end, param):
        t_artss, t_revert = get_time_step_artss(t_end, artss_data_path)
        config_file_name = change_artss(
            param,
            [source_type, temperature_source, random],
            f'{t_artss}_f',
            path=cwd
        )
        client.send_message(create_message(t_revert, config_file_name))
        wait_artss(t_end, artss_data_path)

        field_reader = FieldReader(t_artss, path=artss_data_path)
        ret = comparison_sensor_simulation_data(
            devc_info_thermocouple,
            fds_data,
            field_reader,
            t_artss,
            t_end
        )

        return ret['T']

    # print('org: {f(t_artss, t_revert, dict())}')
    # t = sensor_times[2]
    # pprint(cur)
    # wait_artss(t, artss_data_path)

    # t_artss, t_revert = get_time_step_artss(t, artss_data_path)
    # pprint((t, t_artss, t_revert))
    # field_reader = FieldReader(t_artss, path=artss_data_path)
    # diff_orig = comparison_sensor_simulation_data(devc_info_temperature, fds_data, field_reader, t_artss, t)
    # print(f'org: {diff_orig["T"]}')
    # print(t_revert)

    # with open('map_2.csv', 'w') as out:
    #     for x in range(20, 85, 10):
    #         for y in [delta['y0'] * t for t in range(8)]:
    #             ret = f(t, {'x0': x, 'y0': y})
    #             print(f'map: {x},{y},{ret}')
    #             out.write(f'{x},{y},{ret}\n')

    # return

    for t_sensor in sensor_times[2:]:
        pprint(cur)
        wait_artss(t_sensor, artss_data_path)

        t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path)
        pprint((t_sensor, t_artss, t_revert))
        field_reader = FieldReader(t_artss, path=artss_data_path)
        diff_orig = comparison_sensor_simulation_data(devc_info_thermocouple, fds_data, field_reader, t_artss, t_sensor)
        print(f'org: {diff_orig["T"]}')
        print(t_revert)

        nabla = {}
        for p in keys:
            config_file_name = change_artss(
                {p: cur[p] + delta[p]},
                [source_type, temperature_source, random],
                f'{t_artss}_{p}',
                path=cwd
            )
            client.send_message(create_message(t_revert, config_file_name))
            wait_artss(t_sensor, artss_data_path)

            field_reader = FieldReader(t_artss, path=artss_data_path)
            diff_cur = comparison_sensor_simulation_data(
                devc_info_thermocouple,
                fds_data,
                field_reader,
                t_artss,
                t_sensor
            )

            # calc new nabla
            nabla[p] = (diff_cur['T'] - diff_orig['T']) / delta[p]
            print(f'cur: {diff_cur["T"]}')
            print(nabla[p])

        pprint(nabla)
        hrr_max = 1.0 if nabla['HRR'] < 0 else cur['HRR'] / nabla['HRR']

        if nabla['x0'] < 0:
            x_max = (domain.domain_param['X2'] - cur['x0']) / -nabla['x0']
        elif nabla['x0'] > 0:
            x_max = (cur['x0'] - domain.domain_param['X1']) / nabla['x0']
        else:
            x_max = 0

        # if nabla['y0'] < 0:
        #     y_max = (domain.domain_param['Y2'] - cur['y0']) / -nabla['y0']
        # else:
        #     y_max = (cur['y0'] - domain.domain_param['Y1']) / nabla['y0']

        # if nabla['z0'] < 0:
        #     z_max = (domain.domain_param['Z2'] - cur['z0']) / -nabla['z0']
        # else:
        #     z_max = (cur['z0'] - domain.domain_param['Z1']) / nabla['z0']

        # search direction
        nabla = np.asarray([nabla[x] for x in keys])
        D = np.eye(len(keys))  # -> hesse matrix
        d = np.dot(-D, nabla)

        # calc alpha
        n = 5
        sigma = 0.01
        alpha = 0.9 * min([1.0, hrr_max, x_max]) #, y_max, z_max])
        print(('mins', [1.0, hrr_max, x_max])) #, y_max, z_max]))
        print(f'alpha = {alpha}')
        while n > 0:
            x = np.asarray([cur[x] for x in keys])
            new_para = x + (alpha * d)
            new_para = {
                p: new_para[i]
                for i, p in enumerate(keys)
            }
            pprint(new_para)

            config_file_name = change_artss(
                new_para,
                [source_type, temperature_source, random],
                f'{t_artss}_{n}',
                path=cwd
            )
            print(f'make adjustment with {new_para}')
            client.send_message(create_message(t_revert, config_file_name))
            wait_artss(t_sensor, artss_data_path)

            field_reader = FieldReader(t_artss, path=artss_data_path)
            diff_cur = comparison_sensor_simulation_data(
                devc_info_thermocouple,
                fds_data,
                field_reader,
                t_artss,
                t_sensor
            )

            if diff_cur['T'] < diff_orig['T'] + np.dot(alpha * sigma * nabla, d):
                print(f'found alpha: {alpha}')
                print(f'found x_k+1: {new_para}')
                print(f"al1: {np.dot(alpha * sigma * nabla, d)}")
                print(f"als: {diff_cur['T']} < {diff_orig['T'] + np.dot(alpha * sigma * nabla, d)}")
                cur = copy(new_para)
                break

            alpha = 0.7 * alpha
            n = n - 1
            print(f'ls: {diff_cur["T"]}')


def wait_artss(t, artss_data_path):
    time.sleep(30)  # needed for artss to rewrite meta file
    t_cur = FieldReader.get_t_current(path=artss_data_path)
    while t_cur <= t:
        time.sleep(1)
        t_new = FieldReader.get_t_current(path=artss_data_path)
        if t_new != t_cur:
            t_cur = t_new
            print(f't_cur: {t_cur} wait for {t}')


def change_artss(change: dict, source: list, file_name: str, path='.') -> str:
    # initiate rollback
    source_type, temperature_source, random = data_assimilation.change_heat_source(*source,
                                                                                   changes={'source_type': {},
                                                                                            'temperature_source': change,
                                                                                            'random': {}})
    da = DAFile()
    da.create_config({'u': False, 'v': False, 'w': False, 'p': False, 'T': False, 'C': False})
    da.create_temperature_source_changes(
        source_type=source_type,
        temperature_source=temperature_source,
        random=random)

    config_file_name = os.path.join(path, f'config_{file_name}.xml')
    da.write_xml(config_file_name)
    return config_file_name


def get_time_step_artss(t_sensor: float, artss_data_path: str) -> [float, float]:
    files = FieldReader.get_all_time_steps(path=artss_data_path)
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
    t_revert = files[max(0, estimated_index - 10)]
    return t_artss, t_revert


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

        ijk = artss.get_ijk_from_xyz(dict_devc[d]['XYZ'][0], dict_devc[d]['XYZ'][2], dict_devc[d]['XYZ'][1])
        dict_devc[d]['index']: int = artss.get_index(*ijk)
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


def comparison_sensor_simulation_data(devc_info: dict, sensor_data: pd.DataFrame, field_reader: FieldReader, t_artss: float, t_sensor: float) -> dict:
    diff: dict = {'T': [], 'C': []}
    fields_sim = field_reader.read_field_data()
    for key in devc_info:
        # sensor data
        type_sensor: str = devc_info[key]['type']
        index_sensor: int = devc_info[key]['index']
        value_sensor: float = sensor_data[key][t_sensor]

        value_sim = fields_sim[type_sensor][index_sensor]
        diff[type_sensor].append(kelvin_to_celsius(value_sim) - value_sensor)
    print(f'diff: {diff["T"]}')
    for key in diff:
        if len(diff[key]) == 0:
            continue
        diff[key] = np.sqrt(sum(np.array(diff[key])**2))
    return diff


def calc_single_differences(devc_info: dict, sensor_data: pd.DataFrame, field_reader: FieldReader, t_artss: float, t_sensor: float) -> dict:
    diff: dict = {}
    fields_sim = field_reader.read_field_data()
    tmp_x = []
    tmp_y = []
    tmp_y2 = []
    index = []
    for key in devc_info:
        type_sensor: str = devc_info[key]['type']
        index_sensor: int = devc_info[key]['index']
        value_sensor: float = sensor_data[key][t_sensor]

        value_sim = fields_sim[type_sensor][index_sensor]
        diff[key] = (kelvin_to_celsius(value_sim) - value_sensor)
        tmp_x.append(key)
        index.append(index_sensor)
        tmp_y.append(value_sim)
        tmp_y2.append(value_sensor)

    # print(f' calc single diff {t_artss}, {t_sensor}')
    # print(index)
    # print(tmp_x)
    # print(tmp_y)
    # print(tmp_y2)
    return diff


def kelvin_to_celsius(kelvin):
    return kelvin - 273.5


def plot_differences(fds_data_path: str, fds_input_file_name: str, artss_data_path: str):
    cwd = os.getcwd()
    a_path = os.path.join(artss_data_path, '30')
    xml = XML(FieldReader.get_xml_file_name(path=a_path), path=a_path)
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    devc_info_temperature, devc_info_thermocouple, fds_data = read_fds_data(fds_data_path, fds_input_file_name, domain)
    print(devc_info_temperature)

    sensor_times = fds_data.index
    print(sensor_times)
    for t in range(1, 8):
        t_sensor = sensor_times[t]
        fig, axs = plt.subplots(2, 1)
        for i in [20, 30, 40, 50, 60, 70, 80]:
            t_artss, t_revert = get_time_step_artss(t_sensor, os.path.join(artss_data_path, str(i)))
            field_reader = FieldReader(t_artss, path=os.path.join(artss_data_path, str(i)))
            fields = field_reader.read_field_data()
            diff = calc_single_differences(devc_info_temperature, fds_data, field_reader, t_artss, t_sensor)
            x = []
            y = []
            y_fds = []
            y_artss = []
            tmp_dict = {}
            for k in devc_info_temperature:
                x.append(devc_info_temperature[k]['XYZ'][0])
                y.append(diff[k])
                y_fds.append(fds_data[k][t_sensor])
                if x[-1] in tmp_dict:
                    tmp_dict[x[-1]].append(y[-1])
                else:
                    tmp_dict[x[-1]] = [y[-1]]
                index = domain.get_index(*domain.get_ijk_from_xyz(devc_info_temperature[k]['XYZ'][0], devc_info_temperature[k]['XYZ'][2], devc_info_temperature[k]['XYZ'][1]))
                y_artss.append(kelvin_to_celsius(fields['T'][index]))
            line_x = []
            line_y = []
            for k in tmp_dict:
                line_x.append(k)
                line_y.append(sum(tmp_dict[k]) / len(tmp_dict[k]))
            axs[0].plot(line_x, line_y, label=f'{i} avg: {sum(y) / len(y):.2e}')
            axs[0].scatter(x, y, label=i)
            axs[1].plot(x, y_artss, label=f'{i}')
        axs[0].set_title('difference')
        axs[1].set_title('sensor data')
        axs[1].plot(x, y_fds, label='fds')
        plt.xlabel('position of temperature sensor')
        #axs[0].ylabel('difference [C]')
        #axs[1].ylabel('T [C]')
        axs[0].legend()
        axs[1].legend()
        plt.savefig(f'diff_{t_artss}.pdf')
        plt.close()

def plot_comparison_da():
    cwd = os.getcwd()
    a_path = os.path.join(artss_data_path,'with_da')
    xml = XML(FieldReader.get_xml_file_name(path=a_path), path=a_path)
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    devc_info_temperature, devc_info_thermocouple, fds_data = read_fds_data(fds_data_path, fds_input_file_name, domain)
    print(devc_info_temperature)

    sensor_times = fds_data.index[:100]
    print(sensor_times)
    print('artss path', artss_data_path)
    artss_times = FieldReader.get_all_time_steps(a_path)[:3051]
    print('artss times',artss_times)
    print('sensor times', sensor_times)

    artss_data = {}
    for path in ['with_da', 'without_da']:
        artss_data[path] = {}
        for sensor in devc_info_temperature:
            artss_data[path][sensor] = []
        for t_artss in artss_times:
            field_reader = FieldReader(t_artss, path=os.path.join(artss_data_path, str(x_pos)))
            fields_sim = field_reader.read_field_data()
            for sensor in devc_info_temperature:
                x_pos: float = devc_info_thermocouple[sensor]['XYZ'][0]
                index_sensor: int = devc_info_thermocouple[sensor]['index']
                value_artss = kelvin_to_celsius(fields_sim['T'][index_sensor])
                artss_data[path][sensor].append(value_artss)

    for sensor in devc_info_temperature:
        for path in artss_data:
            plt.plot(artss_times, artss_data[path][sensor], label=path.replace('_', ' '))
        
        fds_values = []
        for t in sensor_times:
            fds_values.append(fds_data[key][t])
        plt.plot(sensor_times, fds_values, label='fds thermocouple')

        fds_values = []
        for t in sensor_times:
            fds_values.append(fds_data[key_temp][t])
        plt.plot(sensor_times, fds_values, label='fds temperature')
        plt.legend() 
        plt.savefig(f'verlauf_sensor_{sensor}.pdf')
        plt.close()

def plot_sensor_data(fds_data_path: str, fds_input_file_name: str, artss_data_path: str):
    cwd = os.getcwd()
    a_path = os.path.join(artss_data_path, '30')
    xml = XML(FieldReader.get_xml_file_name(path=a_path), path=a_path)
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    devc_info_temperature, devc_info_thermocouple, fds_data = read_fds_data(fds_data_path, fds_input_file_name, domain)
    print(devc_info_temperature)

    sensor_times = fds_data.index[:13]
    print(sensor_times)
    print('artss path', artss_data_path)
    artss_times = FieldReader.get_all_time_steps(a_path)
    print('artss times', artss_times)
    print('sensor times', sensor_times)

    tmp_dict = {}
    axis_name = set()
    for d in devc_info_temperature:
        xpos = devc_info_temperature[d]['XYZ'][0]
        axis_name.add(xpos)
        if xpos in tmp_dict:
            tmp_dict[xpos].append(d)
        else:
            tmp_dict[xpos] = [d]

    for key in tmp_dict:
        print(key, tmp_dict[key])

    axis_name = list(axis_name)
    axis_name.sort()

    artss_data = {}
    for x_pos in [20, 30, 40, 50, 60, 70, 80]:
        artss_data[x_pos] = {}
        for sensor_x in axis_name:
            artss_data[x_pos][sensor_x] = [[], []]
        for t_artss in artss_times:
            field_reader = FieldReader(t_artss, path=os.path.join(artss_data_path, str(x_pos)))
            fields_sim = field_reader.read_field_data()
            for sensor_x in axis_name:
                for i in range(len(tmp_dict[sensor_x])):
                    index_sensor: int = devc_info_temperature[tmp_dict[sensor_x][i]]['index']
                    value_artss = kelvin_to_celsius(fields_sim['T'][index_sensor])
                    artss_data[x_pos][sensor_x][i].append(value_artss)

    for index, sensor_x in enumerate(axis_name):
        fig, axs = plt.subplots(2, 1)
        for i in range(len(tmp_dict[sensor_x])):
            key = tmp_dict[sensor_x][i]
            key_temp = key.replace('Thermocouple', 'Temperature')
            for x_pos in [30, 40, 50, 60, 70, 80]:
                axs[i].plot(artss_times, artss_data[x_pos][sensor_x][i], label=f'artss xpos {x_pos}')

            fds_values = []
            for t in sensor_times:
                fds_values.append(fds_data[key][t])
            axs[i].plot(sensor_times, fds_values, label='fds thermocouple')

            fds_values = []
            for t in sensor_times:
                fds_values.append(fds_data[key_temp][t])
            axs[i].plot(sensor_times, fds_values, label='fds temperature')

        plt.xlabel('t [s]')
        plt.ylabel('T [C]')
        plt.legend()
        plt.savefig(f'verlauf_sensor_{sensor_x}.pdf')
        plt.close()
