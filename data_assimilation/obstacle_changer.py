#!/usr/bin/env python3
import os
from tqdm import tqdm
import time
from typing import List, Dict, Tuple

import numpy as np

import TCP_client
from ARTSS import XML, Domain
from data_assimilation import FieldReader, create_message, DAFile
from gradient_based_optimisation import get_time_step_artss
from obstacle import Obstacle


def create_obstacles() -> List[Obstacle]:
    obstacles = []

    obstacle_sw = Obstacle("cube SW")
    obstacle_sw.add_geometry([-6, -4.8, -6, -4.8, -1.2, 1.2])
    obstacle_sw.add_boundary(fields=['u', 'v', 'w'], patches=['front', 'back', 'left', 'right', 'bottom', 'top'],
                             boundary_condition='dirichlet', value=0)
    obstacle_sw.add_boundary(fields=['p'], patches=['front', 'back', 'left', 'right', 'bottom', 'top'],
                             boundary_condition='neumann', value=0)
    obstacles.append(obstacle_sw)

    obstacle_nw = Obstacle("cube NW")
    obstacle_nw.add_geometry([-6, -4.4, 4.4, 6, -1.2, 1.2])
    obstacle_nw.add_boundary(fields=['u', 'v', 'w'], patches=['front', 'back', 'left', 'right', 'bottom', 'top'],
                             boundary_condition='dirichlet', value=0)
    obstacle_nw.add_boundary(fields=['p'], patches=['front', 'back', 'left', 'right', 'bottom', 'top'],
                             boundary_condition='neumann', value=0)
    obstacles.append(obstacle_nw)

    obstacle_no = Obstacle("cube NO")
    obstacle_no.add_geometry([4.8, 6, 4.8, 6, -1.2, 1.2])
    obstacle_no.add_boundary(fields=['u', 'v', 'w'], patches=['front', 'back', 'left', 'right', 'bottom', 'top'],
                             boundary_condition='dirichlet', value=0)
    obstacle_no.add_boundary(fields=['p'], patches=['front', 'back', 'left', 'right', 'bottom', 'top'],
                             boundary_condition='neumann', value=0)
    obstacles.append(obstacle_no)

    obstacle_so = Obstacle("cube SO")
    obstacle_so.add_geometry([4.4, 6, -6, -4.4, -1.2, 1.2])
    obstacle_so.add_boundary(fields=['u', 'v', 'w'], patches=['front', 'back', 'left', 'right', 'bottom', 'top'],
                             boundary_condition='dirichlet', value=0)
    obstacle_so.add_boundary(fields=['p'], patches=['front', 'back', 'left', 'right', 'bottom', 'top'],
                             boundary_condition='neumann', value=0)
    obstacles.append(obstacle_so)
    return obstacles


def obstacle_wonder(artss_data_path: str):
    artss_data_path = os.path.abspath(artss_data_path)
    client = TCP_client.TCPClient()
    client.connect()

    xml = XML(FieldReader.get_xml_file_name(artss_data_path), path=artss_data_path)
    xml.read_xml()

    time_back = 1
    sensor_times = (np.array(range(20)) + 1) / 10 + time_back * xml.get_dt()
    print(sensor_times)
    obstacles = create_obstacles()

    counter = 1
    for t_sensor in sensor_times[:8]:
        wait_artss(t_sensor, artss_data_path)

        t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=xml.get_dt(),
                                                time_back=time_back * xml.get_dt())

        da = DAFile()
        da.create_obstacle_changes([obstacles[counter]], True)
        counter = (counter + 1) % 4
        config_file_name = f'change_obstacle_{int(t_sensor * 10)}.xml'
        config_file_path = os.path.join(artss_data_path, config_file_name)
        da.write_xml(config_file_path, pretty_print=True)
        client.send_message(create_message(t_revert, config_file_path))

    counter = 0
    for t_sensor in sensor_times[8:]:
        wait_artss(t_sensor, artss_data_path)

        t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=xml.get_dt(),
                                                time_back=time_back * xml.get_dt())

        da = DAFile()
        da.create_obstacle_changes([obstacles[counter], obstacles[counter + 2]], True)
        counter = (counter + 1) % 2
        config_file_name = f'change_obstacle_{int(t_sensor * 10)}.xml'
        config_file_path = os.path.join(artss_data_path, config_file_name)
        da.write_xml(config_file_path, pretty_print=True)
        client.send_message(create_message(t_revert, config_file_path))


def full_corridor_doors(artss_data_path: str):
    client = TCP_client.TCPClient()
    client.connect()

    xml = XML(FieldReader.get_xml_file_name(artss_data_path), path=artss_data_path)
    xml.read_xml()

    obstacles = xml.obstacles
    door1, neighbours1 = create_door(obstacles, door_identifier='Room1', name='door1')
    door2, neighbours2 = create_door(obstacles, door_identifier='Room2', name='door2')
    door3, neighbours3 = create_door(obstacles, door_identifier='Room3', name='door3')
    door4, neighbours4 = create_door(obstacles, door_identifier='Room4', name='door4')
    door5, neighbours5 = create_door(obstacles, door_identifier='Room5', name='door5')
    door6, neighbours6 = create_door(obstacles, door_identifier='Room6', name='door6')

    doors = [door1, door2, door3, door4, door5, door6]
    neighbours = [neighbours1, neighbours2, neighbours3, neighbours4, neighbours5, neighbours6]

    time_back = 1
    domain = Domain(xml.domain, xml.computational_domain, obstacles)
    for counter in range(len(doors)):
        print('door', counter)
        # close door
        t_sensor = 11 * xml.get_dt() + time_back * xml.get_dt() + counter * xml.get_dt() * 4
        wait_artss(t_sensor, artss_data_path)
        t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=xml.get_dt(),
                                                time_back=time_back * xml.get_dt())
        domain.add_obstacle(doors[counter])
        domain.print_info()
        da = DAFile()
        reader = FieldReader(t_artss, path=artss_data_path)
        field_changes, fields = set_value(fields=reader.read_field_data(), domain=domain,
                                          obstacle_names=[doors[counter].name],
                                          neighbouring_obstacle_patches={doors[counter].name: neighbours[counter]})
        field_key_changes: List[str] = merge_field_keys(field_changes, field_key_changes=[])
        field_file_name = f'fields_{t_artss:.5e}.hdf'
        FieldReader.write_field_data_keys(file_name=field_file_name, data=fields, field_keys=field_key_changes,
                                          path=artss_data_path)
        da.create_field_changes(field_changes, field_file_name=field_file_name)
        da.write_obstacle_changes(domain.get_obstacles())
        config_file_name = f'full_corridor_obstacle_{int(t_sensor * 10)}.xml'
        config_file_path = os.path.join(artss_data_path, config_file_name)
        da.write_xml(config_file_path, pretty_print=True)
        client.send_message(create_message(t_revert, config_file_name))

        # open door
        t_sensor = 13 * xml.get_dt() + time_back * xml.get_dt() + counter * xml.get_dt() * 4
        wait_artss(t_sensor, artss_data_path)
        t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=xml.get_dt(),
                                                time_back=time_back * xml.get_dt())
        door = domain.remove_obstacle(doors[counter].name)
        domain.print_info()
        da = DAFile()
        reader = FieldReader(t_artss, path=artss_data_path)
        field_changes, fields = set_ambient_temperature(domain=domain, obstacles=[door], value=299.14,
                                                        fields=reader.read_field_data())
        field_key_changes = merge_field_keys(field_changes, field_key_changes=[])
        field_file_name = f'fields_{t_artss:.5e}.hdf'
        FieldReader.write_field_data_keys(file_name=field_file_name, data=fields, field_keys=field_key_changes,
                                          path=artss_data_path)
        da.create_field_changes(field_changes, field_file_name=field_file_name)
        da.write_obstacle_changes(domain.get_obstacles() + [door])
        config_file_name = f'full_corridor_obstacle_{int(t_sensor * 10)}.xml'
        config_file_path = os.path.join(artss_data_path, config_file_name)
        da.write_xml(config_file_path, pretty_print=True)
        client.send_message(create_message(t_revert, config_file_name))


def merge_field_keys(field_changes: Dict[str, bool], field_key_changes: List[str]):
    for f in field_changes:
        if field_changes[f] and f not in field_key_changes:
            field_key_changes.append(f)
    return field_key_changes


def full_corridor_rooms(artss_data_path: str):
    client = TCP_client.TCPClient()
    client.connect()

    xml = XML(FieldReader.get_xml_file_name(artss_data_path), path=artss_data_path)
    xml.read_xml()

    obstacles = xml.obstacles

    domain = Domain(xml.domain, xml.computational_domain, obstacles)
    rooms = [replace_room2(domain)]

    time_back = 1
    for counter in range(len(rooms)):
        print('room', counter + 2)
        # close door
        t_sensor = 11 * xml.get_dt() + time_back * xml.get_dt() + counter * xml.get_dt() * 4
        wait_artss(t_sensor, artss_data_path)
        t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=xml.get_dt(),
                                                time_back=time_back * xml.get_dt())
        removed = []
        for o in rooms[counter][0]:
            domain.add_obstacle(o)
        for o_name in rooms[counter][1]:
            removed.append(domain.remove_obstacle(o_name))

        o_names = [x.name for x in rooms[counter][0]]
        domain.print_info()
        reader = FieldReader(t_artss, path=artss_data_path)
        field_changes, fields = set_value(fields=reader.read_field_data(), domain=domain,
                                          obstacle_names=o_names,
                                          neighbouring_obstacle_patches={})
        field_key_changes: List[str] = merge_field_keys(field_changes, field_key_changes=[])
        field_file_name = f'fields_{t_artss:.5e}.hdf'
        FieldReader.write_field_data_keys(file_name=field_file_name, data=fields, field_keys=field_key_changes,
                                          path=artss_data_path)
        da = DAFile()
        da.create_field_changes(field_changes, field_file_name=field_file_name)
        da.write_obstacle_changes(domain.get_obstacles() + removed)
        config_file_name = f'full_corridor_obstacle_{int(t_sensor * 10)}.xml'
        config_file_path = os.path.join(artss_data_path, config_file_name)
        da.write_xml(config_file_path, pretty_print=True)
        client.send_message(create_message(t_revert, config_file_name))

        # open door
        t_sensor = 13 * xml.get_dt() + time_back * xml.get_dt() + counter * xml.get_dt() * 4
        wait_artss(t_sensor, artss_data_path)
        t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=xml.get_dt(),
                                                time_back=time_back * xml.get_dt())
        for o in removed:
            domain.add_obstacle(o)
        o_names = [x.name for x in removed]
        removed = []
        for o in rooms[counter][0]:
            removed.append(domain.remove_obstacle(o.name))
        domain.print_info()
        da = DAFile()
        reader = FieldReader(t_artss, path=artss_data_path)
        field_changes, fields = set_ambient_temperature(fields=reader.read_field_data(), domain=domain,
                                                        obstacles=removed,
                                                        value=293.15)
        field_key_changes = merge_field_keys(field_changes, field_key_changes=[])
        field_changes, fields = set_value(fields=fields, domain=domain,
                                          obstacle_names=o_names,
                                          neighbouring_obstacle_patches={})
        field_key_changes = merge_field_keys(field_changes, field_key_changes=field_key_changes)

        field_file_name = f'fields_{t_artss:.5e}.hdf'
        FieldReader.write_field_data_keys(file_name=field_file_name, data=fields, field_keys=field_key_changes,
                                          path=artss_data_path)
        da.create_field_changes(field_changes, field_file_name=field_file_name)
        da.write_obstacle_changes(domain.get_obstacles() + removed)
        config_file_name = f'change_obstacle_{int(t_sensor * 10)}.xml'
        config_file_path = os.path.join(artss_data_path, config_file_name)
        da.write_xml(config_file_path, pretty_print=True)
        client.send_message(create_message(t_revert, config_file_name))


def steckler_door(artss_data_path: str):
    """
    experimental setup to close and open the door
    :param artss_data_path: path to where ARTSS was started
    """
    # artss_data_path = os.path.abspath(artss_data_path)
    client = TCP_client.TCPClient()
    client.connect()

    xml = XML(FieldReader.get_xml_file_name(artss_data_path), path=artss_data_path)
    xml.read_xml()

    obstacles = xml.obstacles
    # for o in obstacles:
    #     o.state = 'unmodified'
    door, neighbours = create_door(obstacles)

    time_back = 1

    # close door
    t_sensor = 81 * xml.get_dt() + time_back * xml.get_dt()
    wait_artss(t_sensor, artss_data_path)
    t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=xml.get_dt(),
                                            time_back=time_back * xml.get_dt())
    domain = Domain(xml.domain, xml.computational_domain, obstacles)
    domain.add_obstacle(door)
    domain.print_debug()
    da = DAFile()
    reader = FieldReader(t_artss, path=artss_data_path)
    field_changes, fields = set_value(domain=domain, obstacle_names=[door.name],
                                      neighbouring_obstacle_patches={door.name: neighbours},
                                      fields=reader.read_field_data())
    field_key_changes: List[str] = []
    for f in field_changes:
        if field_changes[f] and f not in field_key_changes:
            field_key_changes.append(f)
    field_file_name = f'fields_{t_artss:.5e}.hdf'
    FieldReader.write_field_data_keys(file_name=field_file_name, data=fields, field_keys=field_key_changes,
                                      path=artss_data_path)
    da.create_field_changes(field_changes, field_file_name=field_file_name)
    da.write_obstacle_changes(domain.get_obstacles())
    config_file_name = f'change_obstacle_{int(t_sensor * 10)}.xml'
    config_file_path = os.path.join(artss_data_path, config_file_name)
    da.write_xml(config_file_path, pretty_print=True)
    client.send_message(create_message(t_revert, config_file_name))

    # open door
    t_sensor = 101 * xml.get_dt() + time_back * xml.get_dt()  # 101
    wait_artss(t_sensor, artss_data_path)
    t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=xml.get_dt(),
                                            time_back=time_back * xml.get_dt())
    door = domain.remove_obstacle(door.name)
    domain.print_debug()
    da = DAFile()
    [field_changes, fields] = set_ambient_temperature(domain=domain, obstacles=[door], value=299.14,
                                                      fields=reader.read_field_data())
    # [field_changes, field_file_name] = set_gradient_x(t_artss=t_artss, domain=domain, obstacle=door,
    #                                                   artss_data_path=artss_data_path)
    field_key_changes: List[str] = []
    for f in field_changes:
        if field_changes[f] and f not in field_key_changes:
            field_key_changes.append(f)
    field_file_name = f'fields_{t_artss:.5e}.hdf'
    FieldReader.write_field_data_keys(file_name=field_file_name, data=fields, field_keys=field_key_changes,
                                      path=artss_data_path)
    da.create_field_changes(field_changes, field_file_name=field_file_name)
    da.write_obstacle_changes(domain.get_obstacles() + [door])
    config_file_name = f'change_obstacle_{int(t_sensor * 10)}.xml'
    config_file_path = os.path.join(artss_data_path, config_file_name)
    da.write_xml(config_file_path, pretty_print=True)
    client.send_message(create_message(t_revert, config_file_name))


def set_ambient_temperature(fields: Dict[str, np.ndarray], obstacles: List[Obstacle], domain: Domain, value: float) -> [
    Dict[str, bool], str]:
    for obstacle in obstacles:
        for j in range(obstacle.index['y1'], obstacle.index['y2'] + 1):
            for k in range(obstacle.index['z1'], obstacle.index['z2'] + 1):
                for i in range(obstacle.index['x1'], obstacle.index['x2'] + 1):
                    index = domain.get_index(i, j, k)
                    fields['T'][index] = value

    return {'T': True}, fields


def set_gradient_x(t_artss: float, obstacle: Obstacle, domain: Domain, artss_data_path: str) -> [Dict[str, bool], str]:
    """
    set gradient in x-direction. interpolate cell values from left neighbour cell to right neighbour cell
    :param t_artss: time of field to be changed
    :param obstacle: obstacle which will be removed
    :param domain: domain of ARTSS
    :param artss_data_path: path to where ARTSS was executed
    :return: A dictionary containing the changed fields and the name of changed field file
    """
    reader = FieldReader(t_artss, path=artss_data_path)
    field_T = reader.read_field_data()['T']
    field_file_name = f'temperature_{t_artss:.5e}'

    for j in range(obstacle.index['y1'], obstacle.index['y2'] + 1):
        for k in range(obstacle.index['z1'], obstacle.index['z2'] + 1):
            start_val = field_T[domain.get_index(obstacle.index['x1'] - 1, j, k)]
            end_val = field_T[domain.get_index(obstacle.index['x2'] + 1, j, k)]
            delta = (end_val - start_val) / (obstacle.index['x2'] - obstacle.index['x1'])
            counter = 1
            print(
                f"start val: {field_T[domain.get_index(obstacle.index['x1'] - 1, j, k)]}, index: {domain.get_index(obstacle.index['x1'] - 1, j, k)}")
            print(
                f"end val: {field_T[domain.get_index(obstacle.index['x2'] + 1, j, k)]}, index: {domain.get_index(obstacle.index['x2'] + 1, j, k)}")
            print(f"delta: {(end_val - start_val) / (obstacle.index['x2'] - obstacle.index['x1'])}")
            for i in range(obstacle.index['x1'], obstacle.index['x2'] + 1):
                index = domain.get_index(i, j, k)
                field_T[index] = start_val + delta * counter
                counter += 1

    FieldReader.write_field_data_keys(file_name=field_file_name,
                                      data={'T': field_T},
                                      field_keys=['T'],
                                      path=artss_data_path)
    return {'T': True}, field_file_name


def set_value(fields: Dict[str, np.ndarray], obstacle_names: List[str], domain: Domain,
              neighbouring_obstacle_patches: Dict[str, Dict[str, List[str]]], value=0) -> [Dict[str, bool],
                                                                                           Dict[str, np.ndarray]]:
    """
    set all cells given obstacle to zero in order to see the obstacle in the simulation
    :param fields: fields to be changed
    :param obstacle_names: name of obstacles which cell values should be changed
    :param domain: domain of ARTSS
    :param neighbouring_obstacle_patches: specify additional cells which should be changed.
    E.g. the patch of a neighbour obstacle
    :return:
    """
    for obstacle_name in obstacle_names:
        for field in fields:
            domain.set_value_of_obstacle_cells(value=value, field=fields[field], obstacle_name=obstacle_name)
            if obstacle_name in neighbouring_obstacle_patches:
                for o_name in neighbouring_obstacle_patches[obstacle_name]:
                    for patch in neighbouring_obstacle_patches[obstacle_name][o_name]:
                        domain.set_value_of_obstacle_patch(value=value, field=fields[field], obstacle_name=o_name,
                                                           patch=patch)

    return dict(zip(fields.keys(), [True] * len(fields))), fields


def create_door(obstacles: List[Obstacle], door_identifier: str = None, ground_coordinate: float = 0, name='door') \
        -> [Obstacle, Dict[str, List[int]]]:
    """
    door measurements are calculated based on the obstacle names "left from door", "right from door", "above door" +
    door identifier. E.g. door identifier = "Room1" results in looking for obstacle names which include
    "left from doorRoom1", "right from doorRoom1" and "above doorRoom1".
    """
    names = ['left from door', 'right from door', 'above door']
    if door_identifier is not None:
        names = [n + door_identifier for n in names]
    door_coordinates = {}
    neighbours = {}
    for obstacle in obstacles:
        if names[2] in obstacle.name:
            obstacle_above = obstacle
            neighbours[obstacle.name] = ['bottom']
        elif names[1] in obstacle.name:
            obstacle_right = obstacle
        elif names[0] in obstacle.name:
            obstacle_left = obstacle

    if obstacle_left.geometry['ox2'] < obstacle_right.geometry['ox1']:
        # door faces z-direction
        door_coordinates['ox1'] = obstacle_left.geometry['ox2']
        door_coordinates['ox2'] = obstacle_right.geometry['ox1']
        door_coordinates['oy1'] = ground_coordinate
        door_coordinates['oy2'] = obstacle_above.geometry['oy1']
        door_coordinates['oz1'] = max(obstacle_left.geometry['oz1'], obstacle_right.geometry['oz1'])
        door_coordinates['oz2'] = min(obstacle_left.geometry['oz2'], obstacle_right.geometry['oz2'])
        neighbours[obstacle_left.name] = ['right']
        neighbours[obstacle_right.name] = ['left']
    else:
        # door faces x-direction
        door_coordinates['ox1'] = max(obstacle_left.geometry['ox1'], obstacle_right.geometry['ox1'])
        door_coordinates['ox2'] = min(obstacle_left.geometry['ox2'], obstacle_right.geometry['ox2'])
        door_coordinates['oy1'] = ground_coordinate
        door_coordinates['oy2'] = obstacle_above.geometry['oy1']
        door_coordinates['oz1'] = obstacle_left.geometry['oz2']
        door_coordinates['oz2'] = obstacle_right.geometry['oz1']
        neighbours[obstacle_left.name] = ['back']
        neighbours[obstacle_right.name] = ['front']

    door = Obstacle(name=name, state='new')
    door.geometry = door_coordinates
    door.boundary = obstacle_above.boundary
    return door, neighbours


def replace_room2(domain: Domain) -> Tuple[List[Obstacle], List[str]]:
    created: List[Obstacle] = []
    room_blocker = Obstacle(name='room2', state='new')
    room_blocker.add_geometry([domain.obstacles['wall between 1 and 2'].geometry['ox1'],
                               domain.obstacles['wall between 2 and 3'].geometry['ox2'],
                               domain.obstacles['wall between 1 and 2'].geometry['oy1'],
                               domain.obstacles['wall between 1 and 2'].geometry['oy2'],
                               domain.obstacles['wall between 1 and 2'].geometry['oz1'],
                               domain.obstacles['wall between 1 and 2'].geometry['oz2'],
                               ])
    created.append(room_blocker)

    deleted: List[str] = ['wall between 1 and 2',
                          'wall between 2 and 3',
                          'room 2 right from doorRoom2',
                          'room 2 above doorRoom2',
                          'room 2 left from doorRoom2']
    return created, deleted


def wait_artss(t_sensor: float, artss_data_path: str):
    t_cur = FieldReader.get_t_current(path=artss_data_path)
    pbar = tqdm(total=t_sensor)
    pbar.update(t_cur)
    while t_cur <= t_sensor:
        time.sleep(1)
        t_new = FieldReader.get_t_current(path=artss_data_path)
        pbar.update(t_new - t_cur)
        t_cur = t_new

    pbar.close()


if __name__ == '__main__':
    # obstacle_wonder(artss_data_path='example')
    # steckler_door(artss_data_path='example')
    # full_corridor_doors(artss_data_path='full_corridor')
    full_corridor_rooms(artss_data_path='full_corridor')
