import time
import TCP_client
from ARTSS import XML, Domain
from data_assimilation import FieldReader


def create_message(t_cur: float, file_name: str) -> bin:
    string_msg = str(t_cur) + "|" + file_name
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

    xml = XML(reader.get_xml_file_name())
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    domain.print_debug()

    t_cur = 0.3
    fields = reader.read_field_data(t_cur)
    field = change_something(domain, fields['T'])
    fields['T'] = field
    file_name = 'test.txt'
    reader.write_field_data(file_name, fields, t_cur)

    client = TCP_client.TCPClient()
    client.connect()
    client.send_message(create_message(t_cur, file_name))
