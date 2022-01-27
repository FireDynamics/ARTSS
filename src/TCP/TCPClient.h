/// \file      TCPClient.h
/// \brief     Based on GitHub Project Asynchronous Sockets for C++ (https://github.com/eminfedar/async-sockets-cpp)
/// \date      June 23, 2021
/// \author    My Linh Wuerzburger

#ifndef FDR_TCPSOCKET_H
#define FDR_TCPSOCKET_H

#include "DllHelper.h"

#include "BaseSocket.h"
#include <string>
#include <functional>
#include <thread>

class EASYSOCKET_API TCPSocket : public BaseSocket {
 public:
    // Event Listeners:
    std::function<void(std::string)> on_message_received;
    std::function<void(const char*, ssize_t)> on_raw_message_received;
    std::function<void(int)> on_socket_closed;

    explicit TCPSocket(FDR_ON_ERROR, int socket_id = -1);

    int send_message(const std::string& message);
    int send_message(const char *bytes, size_t bytes_length);

    void initiate_connection(const std::string& host,
                             uint16_t port,
                             const std::function<void()>& on_connected = []() { },
                             FDR_ON_ERROR);
    void initiate_connection(uint32_t ipv4, uint16_t port,
                             const std::function<void()>& on_connected = []() { },
                             FDR_ON_ERROR);

    void start_listening();

    void set_address_struct(sockaddr_in addr);
    sockaddr_in get_address_struct() const;

 private:
    static void receive(TCPSocket *socket);

    void set_timeout(int seconds);
};

#endif
