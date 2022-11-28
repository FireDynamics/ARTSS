#!/usr/bin/env python3
import os
from tqdm import tqdm
import time
from typing import List

import numpy as np

import TCP_client
from ARTSS import XML
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
        da.create_config({'u': False, 'v': False, 'w': False, 'p': False, 'T': False, 'C': False})
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
        da.create_config({'u': False, 'v': False, 'w': False, 'p': False, 'T': False, 'C': False})
        da.create_obstacle_changes([obstacles[counter], obstacles[counter + 2]], True)
        counter = (counter + 1) % 2
        config_file_name = f'change_obstacle_{int(t_sensor * 10)}.xml'
        config_file_path = os.path.join(artss_data_path, config_file_name)
        da.write_xml(config_file_path, pretty_print=True)
        client.send_message(create_message(t_revert, config_file_path))


def steckler_door(artss_data_path: str):
    artss_data_path = os.path.abspath(artss_data_path)
    client = TCP_client.TCPClient()
    client.connect()

    xml = XML(FieldReader.get_xml_file_name(artss_data_path), path=artss_data_path)
    xml.read_xml()

    obstacles = xml.obstacles
    # for o in obstacles:
    #     o.state = 'unmodified'
    door = create_door(obstacles)

    time_back = 1
    t_sensor = 81 * xml.get_dt() + time_back * xml.get_dt()

    wait_artss(t_sensor, artss_data_path)

    t_artss, t_revert = get_time_step_artss(t_sensor, artss_data_path, dt=xml.get_dt(),
                                            time_back=time_back * xml.get_dt())

    da = DAFile()
    da.create_config({'u': False, 'v': False, 'w': False, 'p': False, 'T': False, 'C': False})
    da.write_obstacle_changes(obstacles + [door], True)
    # da.remove_obstacle(obstacles, door)
    config_file_name = f'change_obstacle_{int(t_sensor * 10)}.xml'
    config_file_path = os.path.join(artss_data_path, config_file_name)
    da.write_xml(config_file_path, pretty_print=True)
    client.send_message(create_message(t_revert, config_file_path))


def create_door(obstacles: List[Obstacle], door_identifier: str = None, ground_coordinate: float = 0) -> Obstacle:
    """
    door measurements are calculated based on the obstacle names "left from door", "right from door", "above door" +
    door identifier. E.g. door identifier = "Room1" results in looking for obstacle names which include
    "left from doorRoom1", "right from doorRoom1" and "above doorRoom1".
    """
    names = ['left from door', 'right from door', 'above door']
    if door_identifier is not None:
        names = [n + door_identifier for n in names]
    door_coordinates = {}
    for obstacle in obstacles:
        if names[2] in obstacle.name:
            obstacle_above = obstacle
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
    else:
        # door faces x-direction
        door_coordinates['ox1'] = max(obstacle_left.geometry['ox1'], obstacle_right.geometry['ox1'])
        door_coordinates['ox2'] = min(obstacle_left.geometry['ox2'], obstacle_right.geometry['ox2'])
        door_coordinates['oy1'] = ground_coordinate
        door_coordinates['oy2'] = obstacle_above.geometry['oy1']
        door_coordinates['oz1'] = obstacle_left.geometry['oz2']
        door_coordinates['oz2'] = obstacle_right.geometry['oz1']

    door = Obstacle(name="door", state='new')
    door.geometry = door_coordinates
    door.boundary = obstacle_above.boundary
    return door


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
    steckler_door(artss_data_path='example')
