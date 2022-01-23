import time
import os
import matplotlib.pyplot as plt

import numpy as np
from numpy import ndarray

import TCP_client
import data_assimilation
from ARTSS import XML, Domain, DAFile
from data_assimilation import FieldReader


def create_message(t_cur: float, config_file_name: str) -> bin:
    string_msg = str(t_cur) + ',' + config_file_name
    return string_msg.encode('utf-8')


def change_something(domain: Domain, field: list) -> list:
    new_field = field.copy()
    i = domain.domain_param['Nx'] // 2
    for j in range(domain.domain_param['Ny']):
        for k in range(domain.domain_param['Nz']):
            new_field[domain.calculate_index(i, j, k)] = 513
    return new_field


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
        xml = XML('da.xml')
    else:
        reader = FieldReader()
        reader.print_header()

        xml = XML(reader.get_xml_file_name())
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
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
            t_cur = reader.get_t_current()
            while t_cur < t:
                time.sleep(5)
                t_cur = reader.get_t_current()
            fields = reader.read_field_data(t_cur)

        field = change_something(domain, fields['T'])
        fields['T'] = field
        field_file_name = f'T_{t_cur}.dat'
        if dry_run:
            data_assimilation.write_field_data(file_name=field_file_name, data=fields,
                                               field_keys=['u', 'v', 'w', 'p', 'T', 'C'])
        else:
            reader.write_field_data(file_name=field_file_name, data=fields)

        source_type, temperature_source, random = \
            change_heat_source(*source,
                               changes={'source_type': {},
                                        'temperature_source': {'x0': float(source[1]['x0']) + 10 * index},
                                        'random': {}})
        da = DAFile()
        da.create_config({'u': False, 'v': False, 'w': False, 'p': False, 'T': True, 'C': False}, field_file_name)
        da.create_temperature_source_changes(
            source_type=source_type,
            temperature_source=temperature_source,
            random=random)
        config_file_name = os.path.join(cwd, f'config_{t}.xml')
        da.write_xml(config_file_name, pretty_print=dry_run)

        if not dry_run:
            client.send_message(create_message(t, config_file_name))

def gradient_tmp():

    reader = FieldReader()
    reader.print_header()
    
    xml = XML(reader.get_xml_file_name())
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()
    
    dt = reader.dt

    t_cur = reader.get_t_current()
    
    n = int(t_cur/dt)
    i = 300
    j = 15
    k = 16
    sensor_data= []
    print('iter')
    f = open('visualisation.dat', 'r')
    for i in range(6):
        f.readline()
    for i in range(1,34):
        print(i, i * dt)
        fields = []
        for i in range(6):
            fields.append(np.fromstring(f.readline(), dtype=np.float, sep=';'))
        sensor_data.append(fields[4][domain.calculate_index(i, j, k)])
        f.readline()
    f.close()
    print("plot")
    f = open('tmp.tmp', 'w')
    for i in sensor_data:
        f.write(str(i) + "\n")
    f.close()
    #f = open('tmp.tmp', 'r')
    #for i in f:
    #    sensor_data.append(float(i))
    plt.plot(sensor_data)
    plt.show()



if __name__ == '__main__':
    #gradient_tmp()
    main(dry_run=False)
