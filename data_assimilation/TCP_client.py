import datetime
import socket
import struct


class TCPClient:
    def __init__(self):
        # Create a TCP/IP socket
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        timeval = struct.pack('ll', 20, 0)
        self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_SNDTIMEO, timeval)
        # Connect the socket to the port where the server is listening
        self.server_address = ('localhost', 7777)

    def set_server_address(self, address: str, port: int):
        self.server_address = (address, port)

    def connect(self):
        self.socket.connect(self.server_address)

    def send_message(self, message: bin):
        # Send data
        self.socket.sendall(message)
        print('message was sent. Waiting for response.')

        # Look for the response
        expected_response = "message was received: Rollback done"
        response = ''
        self.socket.settimeout(30)
        print('time', datetime.datetime.now())
        while len(response) < len(expected_response):
            try:
                reply = self.socket.recv(16)
            except Exception as e:
                err = e.args[0]
                if err == 'timed out':
                    print("an error occurred during the connection:", response, e)
                    break
                else:
                    print('is there an else?', e)
            else:
                response += reply.decode('utf-8')
        if response == expected_response:
            print("message was sent successfully. reply:", response)
        else:
            print("an error occurred during the connection:", response)
