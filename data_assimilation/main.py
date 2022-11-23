import os
import time

import numpy as np
from numpy import ndarray

import TCP_client
import data_assimilation
import gradient_based_optimisation
from ARTSS import XML, Domain, DAFile, DAPackage
from data_assimilation import FieldReader
from obstacle import Obstacle


def create_message(t_cur: float, config_file_name: str) -> bin:
    package = DAPackage(t_cur, config_file_name)
    return package.pack()


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
            t_cur = FieldReader.get_t_current()
            while t_cur < t:
                time.sleep(5)
                t_cur = FieldReader.get_t_current()

            reader = FieldReader(t)
            reader.print_header()
            fields = reader.read_field_data(t)

        field = change_something(domain, fields['T'])
        fields['T'] = field
        field_file_name = f'T_{t_cur}.dat'
        if dry_run:
            for k in fields:
                print((k, len(fields[k])))
            data_assimilation.write_field_data_keys(file_name=field_file_name, data=fields,
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
        da.create_config({'u': False, 'v': False, 'w': False, 'p': False, 'T': True, 'C': False}, field_file_name)
        da.create_temperature_source_changes(
            source_type=source_type,
            temperature_source=temperature_source,
            random=random)
        config_file_name = os.path.join(cwd, f'config_{t}.xml')
        da.write_xml(config_file_name, pretty_print=dry_run)

        if not dry_run:
            client.send_message(create_message(t, config_file_name))
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
    da.create_config({'u': False, 'v': False, 'w': False, 'p': False, 'T': False, 'C': False})
    da.create_obstacle_changes([obstacle], True)
    da.write_xml('change_obstacle.xml', pretty_print=True)


if __name__ == '__main__':
    gradient_based_optimisation.start(artss_data_path='example',
                                      fds_data_path='example/fds_data', fds_input_file_name='tunnel',
                                      artss_path=os.path.join(os.getcwd(), '..'),
                                      parallel=True)
    # main(dry_run=False)
