/// \file      BaseSocket.cpp
/// \brief     Based on GitHub Project Asynchronous Sockets for C++ (https://github.com/eminfedar/async-sockets-cpp)
/// \date      June 23, 2021
/// \author    My Linh Wuerzburger

#include "BaseSocket.h"

#include <cerrno>

std::string BaseSocket::ip_to_string(sockaddr_in addr) {
    char ip[INET_ADDRSTRLEN];
    inet_ntop(AF_INET, &(addr.sin_addr), ip, INET_ADDRSTRLEN);

    return std::string(ip);
}

BaseSocket::BaseSocket(std::function<void(int, std::string)> onError,
                       SocketType socket_type, int socket_id) :
        sock(socket_id) {
    if (socket_id < 0) {
        if ((this->sock = socket(AF_INET, socket_type, 0)) < 0) {
            onError(errno, "Socket creating error.");
        }
    }
}

void BaseSocket::close_socket() {
    if (is_closed) {
        return;
    }

    is_closed = true;
    close(this->sock);
}

std::string BaseSocket::remote_address() {
    return ip_to_string(this->address);
}

int BaseSocket::remote_port() {
    return ntohs(this->address.sin_port);
}
