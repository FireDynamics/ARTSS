import socket
import sys
import time

def config_tcp():
    # Create a TCP/IP socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    # Connect the socket to the port where the server is listening
    server_address = ('localhost', 7777)
    sock.connect(server_address)
    return sock


def create_message(t_cur, file_name):
    string_msg = str(t_cur) + "|" + file_name
    return string_msg.encode('utf-8')


def send_message(sock, message):
    try:
        # Send data
        sock.sendall(message)

        # Look for the response
        amount_received = 0
        amount_expected = len("message was received")

        while amount_received < amount_expected:
            data = sock.recv(16)
            amount_received += len(data)
    finally:
        print("message was sent successfully")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    socket = config_tcp()
    time.sleep(5)
    send_message(socket, create_message(0.1, "test.txt"))
    send_message(socket, create_message(0.2, "test.txt"))
    time.sleep(5)
    send_message(socket, create_message(0.3, "test.txt"))
