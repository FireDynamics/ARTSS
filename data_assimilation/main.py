import time
import os
import TCP_client
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

    for key in random_changes:
        new_random[key] = random_changes[key]

    return new_source_type, new_temperature_source, new_random


def main():
    cwd = os.getcwd()
    print(cwd)

    xml2 = XML('da.xml')
    source = xml2.get_temperature_source()
    source_type, temperature_source, random = \
        change_heat_source(*source,
                           changes={'source_type': {}, 'temperature_source': {'x0': float(source[1]['x0']) + 10},
                                    'random': {}})

    da = DAFile()
    da.create_config({'u': True, 'v': False, 'w': True, 'p': True, 'T': False, 'C': False}, 'test.xml')
    da.create_temperature_source_changes(
        source_type=source_type,
        temperature_source=temperature_source,
        random=random)
    da.write_xml('text.xml', pretty_print=True)

    reader = FieldReader()
    reader.print_header()

    # optional, provides additional information
    xml = XML(reader.get_xml_file_name())
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    t_cur = 0.5
    # fields = reader.read_field_data(t_cur)
    # if len(fields.keys()) > 0:
    #    field = change_something(domain, fields['T'])
    #    fields['T'] = field
    #    field_file_name = 'test.txt'
    #    reader.write_field_data(field_file_name, fields, t_cur)
    #    config_file_name = f'config_{t_cur}.xml'
    #    xml.write_config(config_file_name, ['T'], t_cur)

    client = TCP_client.TCPClient()
    client.connect()

    for t in [0.2, 0.5, 1.0]:
        t_cur = reader.get_t_current()
        while t_cur < t:
            time.sleep(5)
            t_cur = reader.get_t_current()
        config_file_name = os.path.join(cwd, f'config_{t}.xml')
        client.send_message(create_message(t, config_file_name))


if __name__ == '__main__':
    main()
