/// \file      main_data_assimilation.cpp
/// \brief
/// \date      Jul 12, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include <iostream>
#include "TimeIntegration.h"
#include "utility/tinyxml2.h"
#include "utility/Parameters.h"
#include "solver/SolverController.h"
#include "mpi.h"
#include "TCP/TCPServer.h"

#ifdef _OPENACC
#include <openacc.h>
#endif

#define PORT 7777
void server();

int main(int argc, char **argv) {
    MPI_Init(nullptr, nullptr);
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int comm_size, comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    if (num_procs < 2 && comm_rank == 0) {
        // TODO Error Handling
        std::exit(1);
    }

    if (comm_rank == 0) {
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

        auto *sc = new SolverController();

#ifdef _OPENACC
        // Initialise GPU
        acc_device_t dev_type = acc_get_device_type();
        acc_init( dev_type );
#endif

        // Integrate over time and solve numerically
        // Time integration
        TimeIntegration ti(sc);
        ti.run();

        // Clean up
        delete sc;
    }

    if (comm_rank == 1) {
        server();
    }
    return 0;
}

void server() {
    bool simulation_is_running = true;
    MPI_Status status;
    MPI_Request request = MPI_REQUEST_NULL;  // communication is finished
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> logger = Utility::create_logger("TCP Server");
#endif
    TCPServer tcp_server;
    tcp_server.on_new_connection = [&](TCPSocket *new_client) {
#ifndef BENCHMARKING
        logger->info("New client: [{}:{}]", new_client->remote_address(), new_client->remote_port());
#endif
        new_client->on_message_received = [new_client, logger, &request, &status](const std::string &message) {  // message from client
#ifndef BENCHMARKING
            logger->info("received message from client {}:{} => {}", new_client->remote_address(), new_client->remote_port(), message);
#endif
            int flag;
            MPI_Iprobe(1, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
            if (!flag) {
#ifndef BENCHMARKING
                logger->warn("message couldn't be processed as a message is already being processed");
#endif
                MPI_Wait(&request, &status);
            }
            new_client->send_message("message was received");  // send a message back (acknowledgment/error/whatever)
            MPI_Isend(message.c_str(), message.size(), MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
#ifndef BENCHMARKING
            logger->debug("message was sent to rank 0");
#endif
        };

        new_client->on_socket_closed = [new_client, logger](int error_code) {
            // on close, client as well server
#ifndef BENCHMARKING
            logger->info("Socket closed: {}:{} -> {}", new_client->remote_address(), new_client->remote_port(), error_code);
#endif
        };
    };

    // bind the server to a port.
#ifndef BENCHMARKING
    logger->debug("bind server to a port");
#endif
    tcp_server.bind_port(PORT, [logger](int error_code, const std::string &error_message) {
        // BINDING FAILED:
#ifndef BENCHMARKING
        logger->critical("Binding failed - {} : {}", error_code, error_message);
#endif
        // binding wasn't possible, therefore exit ARTSS because something's wrong with the environment
        std::exit(1);
        // TODO Error Handling
    });

    // Start Listening the server.
#ifndef BENCHMARKING
    logger->debug("server: start listening");
#endif
    tcp_server.start_listening([logger](int error_code, const std::string &error_message) {
        // LISTENING FAILED:
#ifndef BENCHMARKING
        logger->critical("Listening failed - {} : {}", error_code, error_message);
#endif
        // listening not possible, therefore exit ARTSS because TCP Server cannot be used
        std::exit(1);
        // TODO Error Handling
    });
#ifndef BENCHMARKING
    logger->debug("completed configuration of TCP server");
#endif

    // TODO: You should do an input loop so the program will not terminated immediately: loop until main thread of MPI says to exit
    while (simulation_is_running) {
    }

    // close the server before exiting the program.
    tcp_server.close_socket();
}
