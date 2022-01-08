import wsgiref.validate
from datetime import datetime
import numpy as np


def get_date_now() -> str:
    return datetime.now().strftime('%a %b %d %H:%M:%S %Y')


class FieldReader:
    def __init__(self):
        self.file_name = 'visualisation.dat'
        self.dt: float
        self.xml_file_name: str
        self.grid_resolution: dict
        self.fields: list
        self.date: datetime
        self.header: str
        self.read_header()

    def read_header(self):
        fname = open(self.file_name, 'r')
        # first line
        line = fname.readline()
        self.dt = float(line.split(';')[3])
        # second line
        line = fname.readline().split(';')[1:]
        self.grid_resolution = {'Nx': int(line[0]), 'Ny': int(line[1]), 'Nz': int(line[2])}
        # third line
        line = fname.readline().strip()
        self.fields = line.split(';')[1:]
        # fourth line
        line = fname.readline().strip()
        self.date = datetime.strptime(line.split(';')[1], '%a %b %d %H:%M:%S %Y')
        # fifth line
        line = fname.readline().strip()
        self.xml_file_name = line.split(';')[1]
        fname.close()

    def print_header(self):
        print(f'dt: {self.dt}')
        print(f'grid resolution: {self.grid_resolution}')
        print(f'fields: {self.fields}')
        print(f'date: {self.date}')
        print(f'xml file name: {self.xml_file_name}')

    def get_line_from_file(self, line_number: int) -> str:
        return self.get_lines_from_file([line_number])[0]

    def get_lines_from_file(self, line_numbers: list) -> list:
        lines = []
        max_val = max(line_numbers)
        file = open(self.file_name, 'r')
        for i, line in enumerate(file):
            if i in line_numbers:
                lines.append(line)
            if i > max_val:
                break
        file.close()
        return lines

    def get_t_current(self) -> float:
        first_line = self.get_line_from_file(0)
        t_cur = first_line.split(';')[1]
        return float(t_cur)

    def get_xml_file_name(self) -> str:
        return self.xml_file_name

    def get_grid_resolution(self) -> dict:
        return self.grid_resolution

    def get_fields(self) -> list:
        return self.fields

    def read_field_data(self, time_step: float) -> dict:
        t_cur = self.get_t_current()
        if time_step > t_cur:
            print(f'cannot read time step {time_step} as the current time step is {t_cur}')
            return {}
        else:
            number_of_fields = len(self.fields)
            steps = int(time_step / self.dt)
            print("t_cur", t_cur, self.dt)
            print("steps", steps)

            starting_line = 5 + number_of_fields * steps
            print(f"starting line for {time_step}: {starting_line}")
            lines = self.get_lines_from_file(list(range(starting_line, starting_line + number_of_fields + 1)))
            fields = {}
            for i in range(number_of_fields):
                fields[self.fields[i]] = np.fromstring(lines[i], dtype=np.float, sep=';')
            return fields

    def write_field_data(self, file_name: str, data: dict, t_cur: float):
        number_of_fields = len(self.fields)
        if number_of_fields != len(data.keys()):
            print(f'wrong number of fields: expected {number_of_fields} got {len(data.keys())}')
        else:
            file = open(file_name, 'w')
            for i in range(number_of_fields):
                field = data[self.fields[i]]
                line = ''
                for number in field:
                    line += f'{number};'
                line += '\n'
                file.write(line)
            file.close()

