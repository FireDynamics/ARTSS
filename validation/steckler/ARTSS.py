from get_cell_type import XML, Domain
import numpy as np


def to_celsius(data):
    return np.array(data) - 272.15


def to_dict(data):
    return {'i': data[0], 'j': data[1], 'k': data[2], 'idx': data[3], 'x': data[4], 'y': data[5], 'z': data[6],
            'vel_x': data[7], 'vel_y': data[8], 'vel_z': data[9], 'vel_div': data[10], 'pressure': data[11],
            'temperature': data[12], 'concentration': data[13], 'sight': data[14], 'turb_visc': data[15],
            'temp_source': data[16]}


class ARTSS:
    def __init__(self):
        self.filename = None
        self.data = None
        self.domain = None

    def add_data_file(self, filename):
        self.filename = filename
        self._extract_artss_data()

    def _extract_artss_data(self):
        f = open(self.filename, 'r')
        i = []
        j = []
        k = []
        idx = []
        x_centre = []
        y_centre = []
        z_centre = []
        temperature = []
        vel_x = []
        vel_y = []
        vel_z = []
        vel_div = []
        pressure = []
        concentration = []
        sight = []
        turb_visc = []
        temp_source = []
        data = [i, j, k, idx, x_centre, y_centre, z_centre, vel_x, vel_y, vel_z, vel_div, pressure, temperature,
                concentration, sight, turb_visc, temp_source]
        f.readline()
        for line in f:
            array = line.split(',')
            if len(array) < 2:
                continue
            for i in range(len(array)):
                data[i].append(float(array[i]))
        f.close()
        self.data = to_dict(data)

    def add_xml(self, filename):
        xml = XML(filename)
        xml.read_xml()
        self.domain = Domain(xml.domain, xml.obstacles)

    def filter_data_on_indices(self, indices_i, indices_j, indices_k, data_type):
        i = indices_i[0]
        filter_data = []
        for j in indices_j:
            tmp = []
            for k in indices_k:
                tmp.append(self.data[data_type][self.domain.calculate_index(i, j, k)])
            filter_data.append(tmp)
        return filter_data

    def clip_data(self, interval_i, interval_j, interval_k, data_type):
        indices_i = list(range(interval_i[0], interval_i[1] + 1))
        indices_j = list(range(interval_j[0], interval_j[1] + 1))
        indices_k = list(range(interval_k[0], interval_k[1] + 1))
        return self.filter_data_on_indices(indices_i, indices_j, indices_k, data_type=data_type)

    def get_indices_from_artss(self, coords_x, coords_y, coords_z, offset_x=0, offset_y=0, offset_z=0):
        artss_i = []
        artss_j = []
        artss_k = []
        for x in coords_x:
            artss_i.append(self.domain.get_i_from_x(offset_x + x))
        for y in coords_y:
            artss_j.append(self.domain.get_j_from_y(offset_y + y))
        for z in coords_z:
            artss_k.append(self.domain.get_k_from_z(offset_z + z))
        return artss_i, artss_j, artss_k

    def print_info(self):
        print("####### ATTENTION! #######\n In ARTSS the x axis is pointing to the left side, the y axis is pointing upwards, the z axis is pointing backwards.\n##########################")
        if self.filename is not None:
            print(f'File data: {self.filename}')
        if self.domain is not None:
            self.domain.print_info()
