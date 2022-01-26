/// \file      TCPClient.cpp
/// \brief     Based on GitHub Project Asynchronous Sockets for C++ (https://github.com/eminfedar/async-sockets-cpp)
/// \date      June 23, 2021
/// \author    My Linh Wuerzburger

#include "TCPClient.h"
#include "../utility/Utility.h"
#include <cstring>

TCPSocket::TCPSocket(std::function<void(int, std::string)> on_error, int socket_id) :
        BaseSocket(on_error, TCP, socket_id) {
}

int TCPSocket::send_message(const std::string &message) {
    return this->send_message(message.c_str(), message.length());
}

int TCPSocket::send_message(const char *bytes, size_t bytes_length) {
    if (this->is_closed) {
        return -1;
    }
    ssize_t sent;
    if ((sent = send(this->sock, bytes, bytes_length, 0)) < 0) {
        auto logger = Utility::create_logger(typeid(this).name());
        logger->error("send");
    }
    return static_cast<int>(sent);
}

void TCPSocket::initiate_connection(const std::string &host, uint16_t port,
                                    const std::function<void()> &on_connected,
                                    std::function<void(int, std::string)> on_error) {
    struct addrinfo hints{}, *res, *it;
    memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_INET;
    hints.ai_socktype = SOCK_STREAM;

    int status;
    if ((status = getaddrinfo(host.c_str(), nullptr, &hints, &res)) != 0) {
        on_error(errno, "Invalid address." + std::string(gai_strerror(status)));
        return;
    }

    for (it = res; it != nullptr; it = it->ai_next) {
        if (it->ai_family == AF_INET) { // IPv4
            memcpy((void *) (&this->address), (void *) it->ai_addr, sizeof(sockaddr_in));
            break;  // for now, just get first ip (ipv4).
        }
    }

    freeaddrinfo(res);

    this->initiate_connection((uint32_t) this->address.sin_addr.s_addr, port, on_connected, on_error);
}

void TCPSocket::initiate_connection(uint32_t ipv4, uint16_t port,
                                    const std::function<void()> &on_connected,
                                    std::function<void(int, std::string)> onError) {
    this->address.sin_family = AF_INET;
    this->address.sin_port = htons(port);
    this->address.sin_addr.s_addr = ipv4;

    this->set_timeout(5);

    // Try to connect.
    if (connect(this->sock, (const sockaddr *) &this->address, sizeof(sockaddr_in)) < 0) {
        onError(errno, "Connection failed to the host.");
        return;
    }

    this->set_timeout(0);

    // Connected to the server, fire the event.
    on_connected();

    // Start listening from server:
    this->start_listening();
}

void TCPSocket::start_listening() {
    // Start listening the socket from thread.
    std::thread receive_listening(receive, this);
    receive_listening.detach();
}

void TCPSocket::set_address_struct(sockaddr_in addr) {
    this->address = addr;
}

sockaddr_in TCPSocket::get_address_struct() const {
    return this->address;
}

void TCPSocket::set_timeout(int seconds) {
    struct timeval tv{};
    tv.tv_sec = seconds;
    tv.tv_usec = 0;

    setsockopt(this->sock, SOL_SOCKET, SO_RCVTIMEO, (char *) &tv, sizeof(tv));
    setsockopt(this->sock, SOL_SOCKET, SO_SNDTIMEO, (char *) &tv, sizeof(tv));
}

void TCPSocket::receive(TCPSocket *socket) {
    char temp_buffer[socket->BUFFER_SIZE];
    ssize_t message_length;

    while ((message_length = recv(socket->sock, temp_buffer, socket->BUFFER_SIZE, 0)) > 0) {
        temp_buffer[message_length] = '\0';
        if (socket->on_message_received) {
            socket->on_message_received(std::string(temp_buffer).substr(0, message_length));
        }
        if (socket->on_raw_message_received) {
            socket->on_raw_message_received(temp_buffer, message_length);
        }
    }

    socket->close_socket();
    if (socket->on_socket_closed) {
        socket->on_socket_closed(errno);
    }
}
