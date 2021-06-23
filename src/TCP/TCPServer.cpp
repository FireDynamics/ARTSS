#include "TCPServer.h"

TCPServer::TCPServer(std::function<void(int, std::string)> onError) : BaseSocket(onError, TCP) {
    int opt = 1;
    setsockopt(this->sock, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(int));
    setsockopt(this->sock, SOL_SOCKET, SO_REUSEPORT, &opt, sizeof(int));
}

void TCPServer::bind_port(int port, std::function<void(int, std::string)> onError) {
    this->bind_port("0.0.0.0", port, onError);
}

void TCPServer::bind_port(const char *address, uint16_t port, std::function<void(int, std::string)> onError) {
    if (inet_pton(AF_INET, address, &this->address.sin_addr) <= 0) {
        onError(errno, "Invalid address. Address type not supported.");
        return;
    }

    this->address.sin_family = AF_INET;
    this->address.sin_port = htons(port);

    if (bind(this->sock, (const sockaddr *) &this->address, sizeof(this->address)) < 0) {
        onError(errno, "Cannot bind the socket.");
        return;
    }
}

void TCPServer::start_listening(std::function<void(int, std::string)> onError) {
    if (listen(this->sock, 10) < 0) {
        onError(errno, "Error: Server can't start_listening the socket.");
        return;
    }

    std::thread accept_thread(accept_connection, this, onError);
    accept_thread.detach();
}

void TCPServer::close_socket() {
    shutdown(this->sock, SHUT_RDWR);

    BaseSocket::close_socket();
}

void TCPServer::accept_connection(TCPServer *server, std::function<void(int, std::string)> onError) {
    sockaddr_in new_socket_info;
    socklen_t new_socket_info_length = sizeof(new_socket_info);

    int new_socket;
    while (!server->is_closed) {
        while ((new_socket = accept(server->sock, (sockaddr *) &new_socket_info, &new_socket_info_length)) < 0) {
            if (errno == EBADF || errno == EINVAL) return;

            onError(errno, "Error while accepting a new connection.");
            return;
        }

        if (!server->is_closed && new_socket >= 0) {
            auto *newSocket = new TCPSocket([](int e, std::string er) {
                FDR_UNUSED(e);
                FDR_UNUSED(er);
            }, new_socket);
            newSocket->set_address_struct(new_socket_info);

            server->on_new_connection(newSocket);
            newSocket->start_listening();
        }
    }
}
