/// \file       main.cpp
/// \brief      Hybrid (CPU and GPU) Version for solving Navier-Stokes equations in 3D
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <unistd.h>
#include <sys/socket.h>
#include <stdlib.h>
#include <netinet/in.h>

#include <thread>
#include <iostream>

#include "TimeIntegration.h"
#include "utility/tinyxml2.h"
#include "utility/Parameters.h"
#include "solver/SolverController.h"

#ifdef _OPENACC
    #include <openacc.h>
#endif

#define PORT 8080
void server(TimeIntegration &ti);

int main(int argc, char **argv) {
    // Initialisation
    // Parameters
    std::string XML_filename;
    auto params = Parameters::getInstance();
    if (argc > 1) {
        XML_filename.assign(argv[1]);
        params->parse(XML_filename);
    } else {
        std::cerr << "XML file missing" << std::endl;
        std::exit(1);
    }

    SolverController *sc = new SolverController();

#ifdef _OPENACC
    // Initialise GPU
    acc_device_t dev_type = acc_get_device_type();
    acc_init( dev_type );
#endif
       

    // Integrate over time and solve numerically
    // Time integration
    TimeIntegration ti(sc);
    server(ti);

    // Clean up
    delete sc;
    return 0;
}

void server(TimeIntegration &ti) {
    int server_fd, new_socket;
    struct sockaddr_in address;
    int opt = 1;
    int addrlen = sizeof(address);
    char buffer[1024] = {0};

    server_fd = socket(AF_INET, SOCK_STREAM, 0);
    setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt));
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(PORT);

    bind(server_fd, (struct sockaddr *)&address, sizeof(address));
    listen(server_fd, 3);
    std::thread thr(&TimeIntegration::run, &ti, 0.0);

    while (true) {
        new_socket = accept(server_fd, (struct sockaddr *)&address, (socklen_t*)&addrlen);
        read(new_socket, buffer, 1024);
        real t = std::stod(buffer);
        ti.stop();
        thr.join();
        thr = std::thread(&TimeIntegration::run, &ti, t);
    }
}
