/// \file      TCPServer.h
/// \brief     Based on GitHub Project Asynchronous Sockets for C++ (https://github.com/eminfedar/async-sockets-cpp)
/// \date      June 23, 2021
/// \author    My Linh Wuerzburger

#ifndef FDR_TCPSERVER_H
#define FDR_TCPSERVER_H

#include "DllHelper.h"

#include "TCPClient.h"
#include <string>
#include <functional>
#include <thread>

class EASYSOCKET_API TCPServer : public BaseSocket {
  public:
    // Event Listeners:
    std::function<void(TCPSocket *)> on_new_connection = [](TCPSocket *sock){FDR_UNUSED(sock)};

    explicit TCPServer(FDR_ON_ERROR);

    // Binding the server.
    void bind_port(int port, FDR_ON_ERROR);
    void bind_port(const char *address, uint16_t port, FDR_ON_ERROR);

    // Start listening the server.
    void start_listening(FDR_ON_ERROR);

    // Overriding close to add shutdown():
    void close_socket() override;

  private:
    static void accept_connection(TCPServer *server, FDR_ON_ERROR);
};

#endif
