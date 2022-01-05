import time
import TCP_client
from ARTSS import XML, Domain
from data_assimilation import FieldReader


def create_message(config_file_name: str, field_file_name: str) -> bin:
    string_msg = config_file_name + ',' + field_file_name
    return string_msg.encode('utf-8')


def change_something(domain: Domain, field: list) -> list:
    new_field = field.copy()
    i = domain.domain_param['Nx'] // 2
    for j in range(domain.domain_param['Ny']):
        for k in range(domain.domain_param['Nz']):
            new_field[domain.calculate_index(i, j, k)] = 513
    return new_field


if __name__ == '__main__':
    reader = FieldReader()
    reader.print_header()

    # optional, provides additional information
    xml = XML(reader.get_xml_file_name())
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    t_cur = 0.08
    fields = reader.read_field_data(t_cur)
    if len(fields.keys()) > 0:
        field = change_something(domain, fields['T'])
        fields['T'] = field
        field_file_name = 'test.txt'
        reader.write_field_data(field_file_name, fields, t_cur)
        config_file_name = 'config.xml'
        xml.write_config(config_file_name, ['T'], t_cur)

        client = TCP_client.TCPClient()
        client.connect()
        client.send_message(create_message(config_file_name, field_file_name))
