import socket
from typing import Type

import pandas as pd

from ARTSS import Domain
from data_assimilation import FieldReader


def start_listener():
    # Set up a TCP/IP server
    tcp_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    # Bind the socket to server address and port 81
    server_address = ('localhost', 81)
    tcp_socket.bind(server_address)

    # Listen on port 81
    tcp_socket.listen(1)

    while True:
        print("Waiting for connection")
        connection, client = tcp_socket.accept()

        print("Connected to client IP: {}".format(client))
        message: str = ''
        # Receive and print data 1024 bytes at a time, as long as the client is sending something
        while True:
            data = connection.recv(1024)
            print("Received data: {}".format(data))

            if not data:
                break
            else:
                message += data

        if message == 'exit':
            connection.close()
        else:
            received_sensor_data(message)

def received_sensor_data(file_name: str):
    data = read_FDS_file(file_name)

def read_FDS_file(file_name: str) -> dict:
    


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
