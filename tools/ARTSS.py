import xml.etree.ElementTree as ET

class XML:
    def __init__(self, filename):
        self.filename = filename
        self.has_obstacles = False
        self.obstacles = []
        self.domain = {}

    def read_xml(self):
        xml_tree = ET.parse(self.filename)
        root = xml_tree.getroot()
        domain_param = root.find('domain_parameters')
        for child in domain_param:
            self.domain[child.tag] = float(child.text)

        obstacles = root.find('obstacles')
        if obstacles.attrib['enabled'] == 'Yes':
            self.has_obstacles = True
            for child in obstacles:
                geometry = child.find('geometry')
                geometry.attrib['name'] = child.attrib['name']
                self.obstacles.append(geometry.attrib)

    def write_config(config_file_name: str, fields: list, t_cur: float):
        # TODO write config file. format:
        # <ARTSS>
        #   <fields_changed u="No" v="No" w="No" p="No" T="Yes" concentration="No"/>
        # </ARTSS>



class Domain:
    def __init__(self, domain_param, obstacles):
        self.domain_param = domain_param
        self.obstacles = []
        self.domain_boundary_cells = {}
        self.calculate_domain()
        self.calculate_obstacles(obstacles)
        self.calculate_indices()

    def calculate_domain(self):
        # from Domain.cpp/h
        self.domain_param['Lx'] = abs(self.domain_param['X2'] - self.domain_param['X1'])
        self.domain_param['Ly'] = abs(self.domain_param['Y2'] - self.domain_param['Y1'])
        self.domain_param['Lz'] = abs(self.domain_param['Z2'] - self.domain_param['Z1'])

        self.domain_param['lx'] = abs(self.domain_param['x2'] - self.domain_param['x1'])
        self.domain_param['ly'] = abs(self.domain_param['y2'] - self.domain_param['y1'])
        self.domain_param['lz'] = abs(self.domain_param['z2'] - self.domain_param['z1'])

        self.domain_param['dx'] = self.domain_param['lx'] / float(self.domain_param['nx'])
        self.domain_param['dy'] = self.domain_param['ly'] / float(self.domain_param['ny'])
        self.domain_param['dz'] = self.domain_param['lz'] / float(self.domain_param['nz'])

        self.domain_param['nx'] = int(self.domain_param['nx'])
        self.domain_param['ny'] = int(self.domain_param['ny'])
        self.domain_param['nz'] = int(self.domain_param['nz'])
        self.domain_param['Nx'] = int(round(self.domain_param['Lx'] / self.domain_param['dx'] + 2))
        self.domain_param['Ny'] = int(round(self.domain_param['Ly'] / self.domain_param['dy'] + 2))
        self.domain_param['Nz'] = int(round(self.domain_param['Lz'] / self.domain_param['dz'] + 2))

        self.domain_param['start x index'] = 1
        self.domain_param['start y index'] = 1
        self.domain_param['start z index'] = 1
        self.domain_param['end x index'] = self.domain_param['Nx'] - 2
        self.domain_param['end y index'] = self.domain_param['Ny'] - 2
        self.domain_param['end z index'] = self.domain_param['Nz'] - 2

    def match_grid(self, obstacle_coordinate, direction):
        # from function get_matching_index in Obstacle.cpp
        return int(round((-self.domain_param[direction.upper() + '1'] + obstacle_coordinate) / self.domain_param['d' + direction.lower()]))

    def calculate_obstacles(self, obstacles):
        for o in obstacles:
            name = o.pop('name')
            dictionary = {'name': name}
            for key in o.keys():
                coord = self.match_grid(float(o[key]), direction=key[1])
                if key[2] == '1':
                    coord += 1
                dictionary[key] = coord
            dictionary['stride x'] = dictionary['ox2'] - dictionary['ox1'] + 1
            dictionary['stride y'] = dictionary['oy2'] - dictionary['oy1'] + 1
            dictionary['stride z'] = dictionary['oz2'] - dictionary['oz1'] + 1
            self.obstacles.append(dictionary)

    def calculate_indices(self):
        self.calculate_domain_boundaries()
        self.calculate_obstacle_indices()

    def calculate_obstacle_indices(self):
        for o in self.obstacles:
            self.calculate_obstacle_boundaries(o)
            self.calculate_obstacle_inner(o)

    def calculate_obstacle_inner(self, o):
        inner_cells = []
        o['inner cells'] = inner_cells
        for i in range(o['ox1'] + 1, o['ox2']):
            for j in range(o['oy1'] + 1, o['oy2']):
                for k in range(o['oz1'] + 1, o['oz2']):
                    inner_cells.append(self.calculate_index(i, j, k))

    def calculate_obstacle_boundaries(self, o):
        front = []
        o['front'] = front
        back = []
        o['back'] = back
        for i in range(o['ox1'], o['ox2'] + 1):
            for j in range(o['oy1'], o['oy2'] + 1):
                front.append(self.calculate_index(i, j, o['oz1']))
                back.append(self.calculate_index(i, j, o['oz2']))

        bottom = []
        o['bottom'] = bottom
        top = []
        o['top'] = top
        for i in range(o['ox1'], o['ox2'] + 1):
            for k in range(o['oz1'], o['oz2'] + 1):
                bottom.append(self.calculate_index(i, o['oy1'], k))
                top.append(self.calculate_index(i, o['oy2'], k))

        left = []
        o['left'] = left
        right = []
        o['right'] = right
        for j in range(o['oy1'], o['oy2'] + 1):
            for k in range(o['oz1'], o['oz2'] + 1):
                left.append(self.calculate_index(o['ox1'], j, k))
                right.append(self.calculate_index(o['ox2'], j, k))

    def calculate_domain_boundaries(self):
        front = []
        self.domain_boundary_cells['front'] = front
        back = []
        self.domain_boundary_cells['back'] = back
        for i in range(self.domain_param['nx']):
            for j in range(self.domain_param['ny']):
                front.append(self.calculate_index(i, j, self.domain_param['start z index']-1))
                back.append(self.calculate_index(i, j, self.domain_param['end z index']+1))

        bottom = []
        self.domain_boundary_cells['bottom'] = bottom
        top = []
        self.domain_boundary_cells['top'] = top
        for i in range(self.domain_param['nx']):
            for k in range(self.domain_param['nz']):
                bottom.append(self.calculate_index(i, self.domain_param['start y index']-1, k))
                top.append(self.calculate_index(i, self.domain_param['end y index']+1, k))

        left = []
        self.domain_boundary_cells['left'] = left
        right = []
        self.domain_boundary_cells['right'] = right
        for j in range(self.domain_param['ny']):
            for k in range(self.domain_param['nz']):
                left.append(self.calculate_index(self.domain_param['start x index']-1, j, k))
                right.append(self.calculate_index(self.domain_param['end x index']+1, j, k))

    def calculate_index(self, i, j, k):
        return i + j * self.domain_param['Nx'] + k * self.domain_param['Nx'] * self.domain_param['Ny']

    def print_info(self):
        # alike to info output in logger
        print(f"Domain size inner cells: {self.domain_param['nx']} {self.domain_param['ny']} {self.domain_param['nz']}\n"
              f"step size (x|y|z): ({self.domain_param['dx']}|{self.domain_param['dy']}|{self.domain_param['dz']})")
        for o in self.obstacles:
            print(f"-- Obstacle {o['name']}\n"
                  f"   strides (x y z): {o['stride x']} {o['stride y']} {o['stride z']}\n"
                  f"   size of slices (Front|Back Bottom|Top Left|Right): {len(o['front'])}|{len(o['back'])} {len(o['bottom'])}|{len(o['top'])} {len(o['left'])}|{len(o['right'])}\n"
                  f"   size of Obstacle: {o['stride x'] * o['stride y'] * o['stride z']}\n"
                  f"   coords (x y z): ({o['ox1']}|{o['ox2']}) ({o['oy1']}|{o['oy2']}) ({o['oz1']}|{o['oz2']})")

    def print_debug(self):
        # alike to debug output in logger
        print(f"Nx: {self.domain_param['Nx']}, Ny: {self.domain_param['Ny']}, Nz: {self.domain_param['Nz']}")
        print(f"nx: {self.domain_param['nx']}, ny: {self.domain_param['ny']}, nz: {self.domain_param['nz']}")
        print(f"X: ({self.domain_param['X1']}|{self.domain_param['X2']}) Y: ({self.domain_param['Y1']}|{self.domain_param['Y2']}) Z: ({self.domain_param['Z1']}|{self.domain_param['Z2']})")
        print(f"Lx: {self.domain_param['Lx']}, Ly: {self.domain_param['Ly']}, Lz: {self.domain_param['Lz']}")
        print(f"lx: {self.domain_param['lx']}, ly: {self.domain_param['ly']}, lz: {self.domain_param['lz']}")
        print(f"dx: {self.domain_param['dx']}, dy: {self.domain_param['dy']}, dz: {self.domain_param['dz']}")
        print(f"X: ({self.domain_param['start x index']}|{self.domain_param['end x index']}) Y: ({self.domain_param['start y index']}|{self.domain_param['end y index']}) Z: ({self.domain_param['start z index']}|{self.domain_param['end z index']})")

    def get_ijk_from_xyz(self, coord_x, coord_y, coord_z):
        return self.get_i_from_x(coord_x), self.get_j_from_y(coord_y), self.get_k_from_z(coord_z)

    def get_i_from_x(self, coord_x):
        return int((coord_x - self.domain_param['X1']) / self.domain_param['dx'])

    def get_j_from_y(self, coord_y):
        return int((coord_y - self.domain_param['Y1']) / self.domain_param['dy'])

    def get_k_from_z(self, coord_z):
        return int((coord_z - self.domain_param['Z1']) / self.domain_param['dz'])

    def get_i(self, index, k, j):
        return int(float(index) - k * self.domain_param['Nx'] * self.domain_param['Ny'] - j * self.domain_param['Nx'])

    def get_j(self, index, k):
        return int((float(index) - k * self.domain_param['Nx'] * self.domain_param['Ny']) / self.domain_param['Nx'])

    def get_k(self, index):
        return int(float(index) / (self.domain_param['Nx'] * self.domain_param['Ny']))

    def get_type(self, index):
        k = self.get_k(index)
        j = self.get_j(index, k)
        i = self.get_i(index, k, j) 
        matches = []
        patches = ['front', 'back', 'bottom', 'top', 'left', 'right']
        for p in patches:
            if index in self.domain_boundary_cells[p]:
                matches.append(f'domain boundary cell ({p})')
        for o in self.obstacles:
            for p in patches:
                if index in o[p]:
                    matches.append(f'obstacle boundary cell ({o["name"]}/{p})')
            if index in o['inner cells']:
                matches.append('obstacle inner cell')
        return f'cell {index} ({i}|{j}|{k}) was found in {matches}'
