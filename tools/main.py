import time
import TCP_client
from ARTSS import XML, Domain
from data_assimilation import FieldReader


def create_message(t_cur, file_name):
    string_msg = str(t_cur) + "|" + file_name
    return string_msg.encode('utf-8')


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    reader = FieldReader()
    reader.print_header()

    xml = XML(reader.get_xml_file_name())
    xml.read_xml()
    domain = Domain(xml.domain, xml.obstacles)
    domain.print_info()
    # domain.print_debug()

    print(reader.read_field_data(0.3))

    # client = TCP_client.TCPClient()
    # client.connect()
    # client.send_message(create_message(0.1, "test.txt"))
    # client.send_message(create_message(0.2, "test.txt"))
    # time.sleep(5)
    # client.send_message(create_message(0.3, "test.txt"))
