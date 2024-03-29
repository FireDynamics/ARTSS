/// \file       main.cpp
/// \brief      Hybrid (CPU and GPU) Version for solving Navier-Stokes equations in 3D
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include "TimeIntegration.h"
#include "domain/DomainData.h"
#include "domain/DomainController.h"
#include "solver/SolverController.h"
#include "utility/settings/Settings.h"

#ifdef _OPENACC
#include <openacc.h>
#endif
#ifdef ASSIMILATION
#include <mpi.h>
#include "TCP/TCPServer.h"
void server();
class tcp_server_error : std::runtime_error { ;
 public:
    explicit tcp_server_error(const std::string &message) : std::runtime_error(message) {}
    explicit tcp_server_error(const char *message) : std::runtime_error(message) {}
};
#endif

int main(int argc, char **argv) {
    try {
#ifdef ASSIMILATION
        MPI_Init(nullptr, nullptr);
        int num_procs;
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

        int comm_size, comm_rank;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

        if (num_procs < 2 && comm_rank == 0) {
            // TODO Error Handling
            std::cout << "Data data_assimilation has to be run with at least 2 processes, "
                         "currently available: " << num_procs << " processes" << std::endl;
            std::exit(1);
        }

        if (comm_rank == 0) {
#endif
            // Initialisation
            // Parameters
            std::string XML_filename;
            if (argc <= 1) {
                std::cerr << "XML file missing" << std::endl;
                std::exit(1);
            }

#ifdef _OPENACC
            // Initialise GPU
            acc_device_t dev_type = acc_get_device_type();
            acc_init(dev_type);
#endif

            Settings::Settings settings = Settings::parse_settings(argv[1]);
#ifdef ASSIMILATION
            MPI_Request request;
            int port = settings.assimilation_parameters.port;
            MPI_Isend(&port, 1, MPI_INT, 1, 10, MPI_COMM_WORLD,&request);
            MPI_Isend(settings.logging_parameters.level.c_str(),
                      static_cast<int>(settings.logging_parameters.level.size()), MPI_CHAR, 1, 0, MPI_COMM_WORLD,
                      &request);
#endif
            size_t multigrid_level = 0;

            const auto &solver = settings.solver_parameters.description;
            if (solver.find("NS") != std::string::npos || solver == SolverTypes::PressureSolver) {
                multigrid_level = settings.solver_parameters.pressure.solver.n_level;
            }
            DomainData::init(settings.physical_parameters, settings.domain_parameters, multigrid_level);
            DomainController::init(settings);

            SolverController *sc = new SolverController(settings);

            // Integrate over time and solve numerically
            // Time integration
            TimeIntegration ti(settings, sc);
            ti.run();

            // Clean up
            delete sc;
#ifdef ASSIMILATION
        }
        if (comm_rank == 1) {
            server();
        }
        MPI_Finalize();
#endif

    } catch (const std::exception &ex) {
        std::cerr << ex.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception encountered" << std::endl;
    }
    return 0;
}

#ifdef ASSIMILATION
void server() {
    bool simulation_is_running = true;
    MPI_Status status;
    MPI_Request request = MPI_REQUEST_NULL;  // communication is finished

    int port;
    MPI_Recv(&port, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD, &status);

#ifndef BENCHMARKING
    MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    int msg_len;
    MPI_Get_count(&status, MPI_CHAR, &msg_len);
    std::vector<char> msg;
    msg.resize(msg_len);
    MPI_Recv(msg.data(), msg_len, MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

    std::string log_file = "tcp_server.log";
    Utility::create_logger(msg.data(), log_file);;
    std::shared_ptr<spdlog::logger> logger = Utility::create_logger("TCPServer");
#endif
    TCPServer tcp_server;
    tcp_server.on_new_connection = [&](TCPSocket *new_client) {
        std::cout << fmt::format("New client: [{}:{}]", new_client->remote_address(), new_client->remote_port())
                  << std::endl;
        logger->info("New client: [{}:{}]", new_client->remote_address(), new_client->remote_port());
        new_client->on_raw_message_received = [new_client, logger, &request, &status](
                const char *message, const int message_size) {  // message from client
            std::cout << fmt::format("received message from client {}:{}", new_client->remote_address(),
                                     reinterpret_cast<const double*>(message)[0]) << std::endl;
            std::cout << fmt::format("received message from client {}:{}", new_client->remote_address(),
                                     new_client->remote_port()) << std::endl;
            logger->info("received message from client {}:{}", new_client->remote_address(), new_client->remote_port());
            std::cout << "message will be sent to rank 0" << std::endl;
            MPI_Isend(message, message_size + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &request);
            logger->debug("message was sent to rank 0");
            std::cout << "message was sent to rank 0" << std::endl;
            int flag = -1;
            while (true) {
                MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
                if (flag) {
                    break;
                }
            }
            int msg_len;
            MPI_Get_count(&status, MPI_CHAR, &msg_len);
            std::vector<char> msg;

            msg.resize(msg_len);
            MPI_Recv(msg.data(), msg_len, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            std::cout << "received reply: " << msg.data() << std::endl;
            logger->debug("received reply: {}", msg.data());
            new_client->send_message(fmt::format("message was received: {}", msg.data()));  // send a message back (acknowledgment/error/whatever)
        };

        new_client->on_socket_closed = [new_client, logger](int error_code) {
            // on close, client as well server
            std::cout << fmt::format("Socket closed: {}:{} -> {}", new_client->remote_address(), new_client->remote_port(),
                         error_code) << std::endl;
            logger->info("Socket closed: {}:{} -> {}", new_client->remote_address(), new_client->remote_port(),
                         error_code);
        };
    };

    // bind the server to a port.
    logger->debug("bind server to a port");
    tcp_server.bind_port(port, [logger](int error_code, const std::string &error_message) {
        // BINDING FAILED:
        logger->critical("Binding failed - {} : {}", error_code, error_message);
        // binding wasn't possible, therefore exit ARTSS because something's wrong with the environment
        throw tcp_server_error(fmt::format("Binding failed. {} : {}", error_code, error_message));
    });

    // Start Listening the server.
    logger->debug("server: start listening");
    tcp_server.start_listening([logger](int error_code, const std::string &error_message) {
        // LISTENING FAILED:
        logger->critical("Listening failed - {} : {}", error_code, error_message);
        // listening not possible, therefore exit ARTSS because TCP Server cannot be used
        throw tcp_server_error(fmt::format("listening failed {} : {}", error_code, error_message));
    });
    logger->debug("completed configuration of TCP server");

    // TODO: You should do an input loop so the program will not terminated immediately: loop until main thread of MPI says to exit
    while (simulation_is_running) {
        MPI_Recv(&simulation_is_running, 1, MPI_LOGICAL, 0, 77, MPI_COMM_WORLD, &status);
    }

    // close the server before exiting the program.
    std::cout << "close server" << std::endl;
    tcp_server.close_socket();
}
#endif
