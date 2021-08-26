import numpy as np
import util as util


def read_data(filename):
    data_dict = {}
    f = open(filename, 'r')
    f.readline()  # header
    line = f.readline()
    while line:
        array = line.split(',')
        name = array.pop(0)
        data_dict[name] = np.array(array).astype(np.float)
        line = f.readline()
    return data_dict


class Steckler:
    def __init__(self):
        self.filename_door_temp = "_door_temperature.csv"
        self.filename_door_vel = "_door_velocity.csv"
        self.filename_experiment_setting = "experiment_setting.csv"
        self.filename_info = "_info.txt"
        self.filename_room = "_room_temperature.csv"
        self.filename_window_temp = "_window_temperature.csv"
        self.filename_window_vel = "_window_velocity.csv"

        self.info_dict = {}
        self.opening_vel_dict = {}
        self.opening_temp_dict = {}
        self.room_temp_dict = {}
        self.experiment_dict = {}
        self.coordinates_dict = {}

    def read_test_case(self, test_case):
        self.read_info(test_case + self.filename_info)

        if self.info_dict['door']:
            opening_name_vel = self.filename_door_vel
            opening_name_temp = self.filename_door_temp
        else:
            opening_name_vel = self.filename_window_vel
            opening_name_temp = self.filename_window_temp
        self.opening_vel_dict = read_data(test_case + opening_name_vel)
        self.opening_temp_dict = read_data(test_case + opening_name_temp)

        self.room_temp_dict = read_data(test_case + self.filename_room)

        self.get_experiment_setting(self.filename_experiment_setting)
        self.calculate_coordinates()

    def read_info(self, filename):
        f = open(filename, 'r')
        lines = list(f.readlines())
        delta_x = float(lines[4].split(' ')[3])
        delta_z = float(lines[5].split(' ')[3])
        door_data = lines[12].strip().split(' ')
        door = door_data[-1].endswith('Door')
        window = not door
        self.info_dict = {
            'delta x': delta_x,
            'delta z': delta_z,
            'door': door,
            'window': window,
            'opening': door_data[-2]
        }

    def get_experiment_setting(self, filename):
        f = open(filename, 'r')
        lines = f.readlines()

        volume = lines[1].split(',')
        self.experiment_dict = {
            'height': float(volume[1]),
            'width': float(volume[2]),
            'depth': float(volume[3]),
            'sill': float(lines[4].split(',')[1]),
            'soffit': float(lines[3].split(',')[1])
        }

        # door data
        for i in range(7, 14):
            door = lines[i].split(',')
            label = door.pop(0)
            self.experiment_dict[label] = np.array(door).astype(np.float)

        # window data
        for i in range(16, 19):
            window = lines[i].split(',')
            label = window.pop(0)
            self.experiment_dict[label] = np.array(window).astype(np.float)

        # gas burner
        for i in range(24, 32):
            gas_burner = lines[i].split(',')
            label = gas_burner.pop(0)
            self.experiment_dict[label] = np.array(gas_burner).astype(np.float)

    def calculate_coordinates(self):
        # opening is in the middle of the right wall
        mid_z = self.experiment_dict['depth'] / 2.
        opening_dim = self.experiment_dict[self.info_dict['opening']]
        self.coordinates_dict = {
            'left jamb': mid_z - opening_dim[1] / 2.,
            'right jamb': mid_z + opening_dim[1] / 2.,
            'sill': self.experiment_dict['soffit'] - opening_dim[0],
            'soffit': self.experiment_dict['soffit']
        }
        self.coordinates_dict['sensors opening x'] = np.arange(self.coordinates_dict['left jamb'] + self.info_dict['delta x'], self.coordinates_dict['right jamb'], self.info_dict['delta x'])
        self.coordinates_dict['sensors opening y'] = [self.experiment_dict['width'] - 0.102 / 2]
        self.coordinates_dict['sensors opening z'] = np.arange(self.coordinates_dict['sill'] + 0.057, self.coordinates_dict['soffit'], self.info_dict['delta z'])

        self.coordinates_dict['sensors room x'] = [self.experiment_dict['width'] - 0.305]
        self.coordinates_dict['sensors room y'] = [self.experiment_dict['depth'] - 0.305]
        self.coordinates_dict['sensors room z'] = np.arange(self.coordinates_dict['sill'] + 0.057, self.experiment_dict['height'], self.info_dict['delta z'])

    def print_info(self):
        print("info")
        util.print_dict(self.info_dict)
        print("\ndoor_vel")
        util.print_dict(self.opening_vel_dict)
        print("\ndoor_temp")
        util.print_dict(self.opening_temp_dict)
        print("\nroom temp")
        util.print_dict(self.room_temp_dict)
        print("\nsetting")
        util.print_dict(self.experiment_dict)
        print("\ncoords")
        util.print_dict(self.coordinates_dict)
        print(
            "\n####### ATTENTION! #######\n In the steckler case the x axis is pointing backwards, the y axis is pointing to the left side, the z axis is pointing upward.\n##########################")
