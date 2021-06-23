/// \file      BaseSocket.h
/// \brief     Based on GitHub Project Asynchronous Sockets for C++ (https://github.com/eminfedar/async-sockets-cpp)
/// \date      June 23, 2021
/// \author    My Linh Wuerzburger

#ifndef FDR_BASESOCKET_H
#define FDR_BASESOCKET_H

#include "DllHelper.h"

#if defined(__linux__) || defined(__APPLE__)

#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>
#include <netdb.h>

#elif _WIN32
#include <winsock32.h>
#endif

#include <string>
#include <functional>

#define FDR_UNUSED(expr){ (void)(expr); }
#define FDR_ON_ERROR \
std::function<void(int, std::string)> onError = [](int errorCode, const std::string &errorMessage) { \
    { (void) (errorCode); }; \
    { (void) (errorMessage); } \
}

class EASYSOCKET_API BaseSocket {
  protected:
    int sock = 0;

    static std::string ip_to_string(sockaddr_in addr);

  public:
    const uint16_t BUFFER_SIZE = 0xFFFF;
    enum EASYSOCKET_API SocketType {
        TCP = SOCK_STREAM,
        UDP [[maybe_unused]] = SOCK_DGRAM
    };

    sockaddr_in address;
    bool is_closed = false;

    explicit BaseSocket(FDR_ON_ERROR, SocketType socket_type = TCP, int socket_id = -1);

    virtual void close();

    std::string remote_address();

    int remote_port();

    int file_descriptor() const { return this->sock; }
};

#endif
