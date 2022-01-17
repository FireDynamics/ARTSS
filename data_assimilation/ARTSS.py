import xml.etree.ElementTree as ET
import xml.dom.minidom as md


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


class DAFile:
    def __init__(self):
        self.xml_root = ET.Element('ARTSS')

    def create_temperature_source_changes(self, source_type: str, dir: str, dissipation: bool,
                                          temperature_source: dict,
                                          random: dict):
        # create temperature source part
        # <ARTSS>
        #   <source type="ExplicitEuler" dir="y" temp_fct="Gauss" dissipation="No" random="Yes">
        #     <HRR> 25000. </HRR> <!-- Total heat release rate( in kW) -->
        #     <cp> 1.023415823 </cp> <!-- specific heat capacity( in kJ / kgK)-->
        #     <x0> 40. </x0>
        #     <y0> -3. </y0>
        #     <z0> 0. </z0>
        #     <sigma_x> 1.0 </sigma_x>
        #     <sigma_y> 1.5 </sigma_y>
        #     <sigma_z> 1.0 </sigma_z>
        #     <tau> 5. </tau>
        #     <random absolute="Yes" custom_seed="Yes" custom_steps="Yes">
        #       <seed> 0 </seed>
        #       <step_size> 0.1 </step_size>
        #       <range> 1 </range>
        #     </random>
        #   </source>
        # </ARTSS>
        source = ET.SubElement(self.xml_root, 'source', type=source_type, dir=dir,
                               dissipation='Yes' if dissipation else 'No', random='Yes' if random['enabled'] else 'No')
        for key in temperature_source:
            ET.SubElement(source, key).text = str(temperature_source[key])

        if random['enabled']:
            random_tree = ET.SubElement(source, 'random', absolute='Yes' if random['absolute'] else 'No',
                                        custom_seed='Yes' if random['custom_seed'] else 'No',
                                        custom_steps='Yes' if random['custom_steps'] else 'No')
            ET.SubElement(random_tree, 'range').text = str(random['range'])
            if random['custom_steps']:
                ET.SubElement(random_tree, 'step_size').text = str(random['step_size'])
            if random['custom_seed']:
                ET.SubElement(random_tree, 'seed').text = str(random['seed'])

    def create_config(self, fields: dict, field_file_name=''):
        # create config file. format:
        # <ARTSS>
        #   <fields_changed u="No" v="No" w="No" p="No" T="Yes" concentration="No" filename="field_file_name"/>
        # </ARTSS>
        changed = False
        field_values = {}
        for key in fields.keys():
            if fields[key]:
                field_values[key] = 'Yes'
                changed = True
            else:
                field_values[key] = 'No'

        if changed:
            ET.SubElement(self.xml_root, 'fields_changed', u=field_values['u'], v=field_values['v'],
                          w=field_values['w'], p=field_values['p'],
                          T=field_values['T'], concentration=field_values['C'], filename=field_file_name)
        else:
            ET.SubElement(self.xml_root, 'fields_changed', u=field_values['u'], v=field_values['v'],
                          w=field_values['w'], p=field_values['p'],
                          T=field_values['T'], concentration=field_values['C'])

    # write xml to file
    def write_xml(self, config_file_name: str, pretty_print=False):
        tree = ET.ElementTree(self.xml_root)
        tree.write(config_file_name)
        print(config_file_name)
        if pretty_print:
            dom = md.parse(config_file_name)
            pretty_format = dom.toprettyxml()
            f = open(config_file_name, 'w')
            f.write(pretty_format)
            f.close()


class Domain:
    def __init__(self, domain_param, obstacles):
        self.domain_param = domain_param
        self.obstacles = []
        self.domain_boundary_cells = {}
        self.control_computational_domain()
        self.calculate_domain()
        self.calculate_obstacles(obstacles)
        self.calculate_indices()

    def control_computational_domain(self):
        if "x1" not in self.domain_param.keys():
            self.domain_param['x1'] = self.domain_param['X1']
            self.domain_param['y1'] = self.domain_param['Y1']
            self.domain_param['z1'] = self.domain_param['Z1']
            self.domain_param['x2'] = self.domain_param['X2']
            self.domain_param['y2'] = self.domain_param['Y2']
            self.domain_param['z2'] = self.domain_param['Z2']

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
        return int(round((-self.domain_param[direction.upper() + '1'] + obstacle_coordinate) / self.domain_param[
            'd' + direction.lower()]))

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
                front.append(self.calculate_index(i, j, self.domain_param['start z index'] - 1))
                back.append(self.calculate_index(i, j, self.domain_param['end z index'] + 1))

        bottom = []
        self.domain_boundary_cells['bottom'] = bottom
        top = []
        self.domain_boundary_cells['top'] = top
        for i in range(self.domain_param['nx']):
            for k in range(self.domain_param['nz']):
                bottom.append(self.calculate_index(i, self.domain_param['start y index'] - 1, k))
                top.append(self.calculate_index(i, self.domain_param['end y index'] + 1, k))

        left = []
        self.domain_boundary_cells['left'] = left
        right = []
        self.domain_boundary_cells['right'] = right
        for j in range(self.domain_param['ny']):
            for k in range(self.domain_param['nz']):
                left.append(self.calculate_index(self.domain_param['start x index'] - 1, j, k))
                right.append(self.calculate_index(self.domain_param['end x index'] + 1, j, k))

    def calculate_index(self, i, j, k):
        return i + j * self.domain_param['Nx'] + k * self.domain_param['Nx'] * self.domain_param['Ny']

    def print_info(self):
        # alike to info output in logger
        print(
            f"Domain size inner cells: {self.domain_param['nx']} {self.domain_param['ny']} {self.domain_param['nz']}\n"
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
        print(
            f"X: ({self.domain_param['X1']}|{self.domain_param['X2']}) Y: ({self.domain_param['Y1']}|{self.domain_param['Y2']}) Z: ({self.domain_param['Z1']}|{self.domain_param['Z2']})")
        print(f"Lx: {self.domain_param['Lx']}, Ly: {self.domain_param['Ly']}, Lz: {self.domain_param['Lz']}")
        print(f"lx: {self.domain_param['lx']}, ly: {self.domain_param['ly']}, lz: {self.domain_param['lz']}")
        print(f"dx: {self.domain_param['dx']}, dy: {self.domain_param['dy']}, dz: {self.domain_param['dz']}")
        print(
            f"X: ({self.domain_param['start x index']}|{self.domain_param['end x index']}) Y: ({self.domain_param['start y index']}|{self.domain_param['end y index']}) Z: ({self.domain_param['start z index']}|{self.domain_param['end z index']})")

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
