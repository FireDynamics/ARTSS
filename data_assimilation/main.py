import multiprocessing
import os
import time
from pprint import pprint
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
from numpy import ndarray

import ARTSS
import TCP_client
import data_assimilation
import downhill_simplex
import fds_utility
import gradient_based_optimisation
import utility
from ARTSS import XML, Domain
from data_assimilation import FieldReader, DAFile
from obstacle import Obstacle
from utility import log


def change_something(domain: Domain, field: list) -> list:
    new_field = field.copy()
    i = domain.domain_param['Nx'] // 2
    for j in range(domain.domain_param['Ny']):
        for k in range(domain.domain_param['Nz']):
            new_field[domain.calculate_index(i, j, k)] = 513
    return new_field


def create_gradient_field(Nx: int, Ny: int, Nz: int) -> ndarray:
    field = np.zeros(Nx * Ny * Nz)
    counter = 0
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                field[counter] = i
                counter += 1
    return field


def main(dry_run=False):
    cwd = os.getcwd()
    print(cwd)

    if dry_run:
        xml = XML('da.xml', path='example')
    else:
        xml = XML(FieldReader.get_xml_file_name(path='example'), path='example')

    xml.read_xml()
    domain = Domain(domain_param=xml.domain, obstacles=xml.obstacles, enable_computational_domain=xml.computational_domain)
    domain.print_info()
    domain.print_debug()

    source = xml.get_temperature_source()

    if not dry_run:
        client = TCP_client.TCPClient()
        client.connect()

    for index, t in enumerate([0.2, 0.5, 1.0], start=1):
        if dry_run:
            t_cur = t
            fields = {'T': create_gradient_field(domain.domain_param['Nx'],
                                                 domain.domain_param['Ny'],
                                                 domain.domain_param['Nz'])}
        else:
            t_cur = FieldReader.get_t_current()
            while t_cur < t:
                time.sleep(5)
                t_cur = FieldReader.get_t_current()

            reader = FieldReader(t)
            reader.print_header()
            fields = reader.read_field_data()

        field = change_something(domain, fields['T'])
        fields['T'] = field
        field_file_name = f'T_{t_cur}.dat'
        if dry_run:
            for k in fields:
                print((k, len(fields[k])))
            FieldReader.write_field_data_keys(file_name=field_file_name, data=fields,
                                              field_keys=['u', 'v', 'w', 'p', 'T', 'C'])
        else:
            reader.write_field_data(file_name=field_file_name, data=fields)

        source_type, temperature_source, random = \
            data_assimilation.change_heat_source(*source,
                                                 changes={'source_type': {},
                                                          'temperature_source': {
                                                              'x0': float(source[1]['x0']) + 10 * index},
                                                          'random': {}})
        da = DAFile()
        da.create_field_changes({'T': True}, field_file_name)
        da.create_temperature_source_changes(
            source_type=source_type,
            temperature_source=temperature_source,
            random=random)
        config_file_name = os.path.join(cwd, f'config_{t}.xml')
        da.write_xml(config_file_name, pretty_print=dry_run)

        if not dry_run:
            client.send_message(data_assimilation.create_message(t, config_file_name))
            time.sleep(10)


def tmp():
    obstacle = Obstacle("cube")
    obstacle.add_geometry([1, 2, 2, 3, 3, 4])
    obstacle.add_boundary(fields=['u', 'v', 'w'], patches=['front', 'back', 'left', 'right', 'bottom', 'top'],
                          boundary_condition='dirichlet', value=0)
    obstacle.add_boundary(fields=['p'], patches=['front', 'right', 'bottom', 'top'], boundary_condition='neumann',
                          value=1)
    obstacle.add_boundary(fields=['p'], patches=['back'], boundary_condition='neumann', value=10)
    obstacle.add_boundary(fields=['p'], patches=['left'], boundary_condition='dirichlet', value=11)
    da = DAFile()
    da.create_obstacle_changes([obstacle], True)
    da.write_xml('change_obstacle.xml', pretty_print=True)


def initial_guess_heat_source_position(time: float, fds_data: pd.DataFrame, temp_devices: Dict[str, any], threshold=22):
    """
    determine position of heat source based on the average of all sensors which values are higher than threshold
    :param time:
    :param fds_data:
    :param temp_devices:
    :param threshold:
    :return:
    """
    positions = {'x': [], 'y': [], 'z': []}
    for sensor in temp_devices:
        if fds_data[sensor][time] > threshold:
            x, y, z = temp_devices[sensor]['XYZ']
            positions['x'].append(x)
            positions['y'].append(z)  # in ARTSS and FDS y and z are swapped
            positions['z'].append(y)

    x = sum(positions['x']) / len(positions['x'])
    y = sum(positions['y']) / len(positions['y'])
    z = sum(positions['z']) / len(positions['z'])
    return x, y, z


def start_simulation(fds_data_path: str, fds_input_file_name: str, artss_data_path: str, artss_input_file_name: str, artss_root_path: str, port=7777, time_difference: float = 0):
    domain, xml = parse_artss_input_file(artss_data_path, artss_input_file_name=artss_input_file_name)

    fds_data, devices, sensor_times = parse_fds_data(fds_data_path=fds_data_path, fds_input_file_name=fds_input_file_name, time_difference=time_difference, domain=domain)
    pos_heat_source = initial_guess_heat_source_position(sensor_times[0], fds_data, devices['temperature'])

    file_da = open(os.path.join(artss_data_path, 'da_details.csv'), 'w')
    file_debug = open(os.path.join(artss_data_path, 'da_debug_details.dat'), 'w')

    heat_source = xml.get_temperature_source()
    delta = {
        'HRR': float(heat_source['temperature_source']['HRR']) * 0.05,
        'x0': 1.5,
        'y0': domain.domain_param['dy'] * 1,
        'z0': 0.5
    }

    cur = {
        'HRR': float(heat_source['temperature_source']['HRR']),
        'x0': pos_heat_source[0],
        'y0': pos_heat_source[1],
        'z0': pos_heat_source[2]
    }

    start_file_name = artss_input_file_name[:-4] + '_initial_guess.xml'
    ARTSS.change_xml(change={'x0': cur['x0'], 'z0': cur['z0']},
                     input_file=os.path.join(artss_data_path, artss_input_file_name), output_file=start_file_name,
                     artss_data_path=artss_data_path, artss_root_path=artss_root_path,
                     file_debug=file_debug)
    job = multiprocessing.Process(target=ARTSS.start_new_instance, args=(start_file_name, artss_data_path, os.path.join(artss_root_path, 'build', 'bin'), 'artss_data_assimilation_serial'))
    job.start()
    client = set_up_client(port=port)

    keys = ['HRR', 'x0', 'z0']

    log(f'cur: {cur}', file_debug)
    log(f'delta: {delta}', file_debug)
    log(f'keys: {keys}', file_debug)

    devc_info_temperature = devices['temperature']
    print('devices:')
    pprint(devc_info_temperature)

    downhill_simplex.opt_scipy(client=client, file_da=file_da, file_debug=file_debug,
                               sensor_times=sensor_times,
                               devc_info=devc_info_temperature, fds_data=fds_data,
                               artss_data_path=artss_data_path,
                               domain=domain, heat_source=heat_source,
                               cur=cur, delta=delta, keys=keys, n_iterations=5, artss=xml)


def set_up_client(ip_address: str = 'localhost', port: int = 7777) -> TCP_client:
    """
    set up TCP client and connect to ARTSS
    :param ip_address: IP address for connecting to the running simulation (default: localhost)
    :param port: port for connecting to the running simulation, can be found in input file (default: 7777)
    :return: connected TCP client
    """
    client = TCP_client.TCPClient()
    client.set_server_address(ip_address, port)
    client.connect()
    return client


def parse_artss_input_file(artss_data_path: str, artss_input_file_name: str = None, quiet=False) -> Tuple[Domain, XML]:
    artss_data_path = os.path.abspath(artss_data_path)
    if artss_input_file_name is None:
        artss_input_file_name = FieldReader.get_xml_file_name(artss_data_path)
    xml = XML(artss_input_file_name, path=artss_data_path)
    xml.read_xml()
    domain = Domain(domain_param=xml.domain, obstacles=xml.obstacles,
                    enable_computational_domain=xml.computational_domain)
    if not quiet:
        domain.print_info()
        domain.print_debug()
    return domain, xml


def parse_fds_data(fds_data_path: str, fds_input_file_name: str, domain: Domain, time_difference: float = 2, threshold=22, keys=['temperature']) -> [pd.DataFrame, Dict[str, Dict[str, any]],
                                                                                                                                                     List[float]]:
    devices, fds_data = fds_utility.read_fds_data(fds_data_path, fds_input_file_name, domain)
    starting_time = fds_utility.get_starting_time(fds_data, threshold, keys)
    starting_time_index = list(fds_data.index).index(starting_time)
    sensor_times = fds_data.index[starting_time_index:]
    res = [sensor_times[0]]
    for i in sensor_times[1:]:
        if i > res[-1] + time_difference:
            res.append(i)
    return fds_data, devices, res


def set_up(fds_data_path: str, fds_input_file_name: str, artss_data_path: str, port: int = 7777, time_difference: float = 0) \
        -> [TCP_client, Domain, XML, List[float], Dict[str, Dict[str, any]], pd.DataFrame, float]:
    """
    set up necessary objects for data assimilation
        :param fds_data_path: path to directory where the FDS and devc file are
    :param fds_input_file_name: name of FDS file (without ending .fds)
    :param artss_data_path: path to where the .vis directory of the simulation is
    :param port: port defined in the XML used by the simulation
    :param time_difference: restrict sensor times, two times have to be at least time_difference apart (in seconds)
    :return: TCPClient: connecting to ARTSS
             Domain: information about running simulation
             XML: parsed ARTSS input file
             sensor times
             devices: data about devices from FDS input file combined with location data from ARTSS
             fds_data: content of FDS device file
             offset: time when the first (temperature) sensor in FDS rises about a certain threshold (22C)
    """
    client = set_up_client(port=port)
    domain, xml = parse_artss_input_file(artss_data_path)
    fds_data, devices, sensor_times = parse_fds_data(fds_data_path=fds_data_path, fds_input_file_name=fds_input_file_name, time_difference=time_difference, domain=domain)
    print('sensor times:', fds_data.index)
    return client, domain, xml, sensor_times, devices, fds_data


if __name__ == '__main__':
    start_simulation(fds_data_path='example/fds_data', fds_input_file_name='tunnel',
                     artss_data_path='../thesis', artss_input_file_name='tunnel.xml', artss_root_path='..',
                     time_difference=2)
    # gradient_based_optimisation.start(artss_data_path='example',
    #                                   fds_data_path='example/fds_data', fds_input_file_name='tunnel',
    #                                   artss_path=os.path.join(os.getcwd(), '..'),
    #                                   parallel=True)
    # downhill_simplex.start(artss_data_path='../thesis/run3',
    #                        fds_data_path='example/fds_data', fds_input_file_name='tunnel',
    #                        parallel=True,
    #                        port=7777)
    # main(dry_run=False)
