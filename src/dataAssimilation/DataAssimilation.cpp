/// \file      DataAssimilation.cpp
/// \brief
/// \date      May 18, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include "DataAssimilation.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#ifdef ASSIMILATION
#include "../TCP/TCPServer.h"
#include "mpi.h"
#endif

DataAssimilation::DataAssimilation(const FieldController &field_controller) : m_field_controller(field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    Parameters *params = Parameters::getInstance();
    m_assimilated = (params->get("data_assimilation/enabled") == "Yes");
    if (m_assimilated) {
        std::string init = params->get("data_assimilation/class_name");
        if (init == AssimilationMethods::Standard) {
#ifndef BENCHMARKING
            m_logger->debug("found data assimilation class {}", init);
#endif
            m_func = new FieldIO();
        } else {
#ifndef BENCHMARKING
            m_logger->critical("Data Assimilation class {} is not defined", init);
#endif
            std::exit(1);
            // TODO Error Handling
        }
#ifdef ASSIMILATION
        config_MPI();
#else
        m_logger->critical("Data Assimilation can be only executed with its respective executable");
        std::exit(1);
#endif
    }
}

void DataAssimilation::assimilate(std::string &file_name) {
    auto *u_current = new Field(FieldType::U);
    auto *v_current = new Field(FieldType::V);
    auto *w_current = new Field(FieldType::W);
    auto *p_current = new Field(FieldType::P);
    auto *T_current = new Field(FieldType::T);
    auto *C_current = new Field(FieldType::RHO);

    m_t_cur = m_func->read(file_name, u_current, v_current, w_current, p_current, T_current, C_current);
    m_field_controller.replace_data(u_current, v_current, w_current, p_current, T_current, C_current);

    delete u_current;
    delete v_current;
    delete w_current;
    delete p_current;
    delete T_current;
    delete C_current;
}

void DataAssimilation::save_data(real t_cur) {
    if (m_assimilated) {
        auto u = m_field_controller.get_field_u_data();
        auto v = m_field_controller.get_field_v_data();
        auto w = m_field_controller.get_field_w_data();
        auto p = m_field_controller.get_field_w_data();
        auto T = m_field_controller.get_field_T_data();
        auto C = m_field_controller.get_field_concentration_data();

        m_func->write(t_cur, u, v, w, p, T, C);
    }
}

real DataAssimilation::get_new_time_value() const {
    return m_t_cur;
}

void DataAssimilation::initiate_rollback() {
    m_rollback = false;
    std::string file_name = "test";
    assimilate(file_name);
}

#ifdef ASSIMILATION
void DataAssimilation::config_MPI() {
    m_logger->info("config MPI");
    MPI_Init(NULL, NULL);
    m_logger->info("config MPI Init");

    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    m_logger->info("num_procs {}", num_procs);

    if (num_procs < 2) {
#ifndef BENCHMARKING
        m_logger->critical("Data Assimilation needs at least 2 processes, current number of processes: {}", num_procs);
#endif
        // TODO Error Handling
        std::exit(1);
    }

    int comm_size, comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    if (comm_rank == 1) {
        TCPServer tcp_server;
        tcp_server.on_new_connection = [&](TCPSocket *new_client) {
            std::cout << "New client: [";
            std::cout << new_client->remote_address() << ":" << new_client->remote_port() << "]" << std::endl;
            // TODO messages to logger

            new_client->on_message_received = [new_client](const std::string& message) {  // message from client
                std::cout << new_client->remote_address() << ":" << new_client->remote_port() << " => " << message << std::endl;

                new_client->send_message("OK!");  // send a message back (acknowledgment/error/whatever)
                // TODO process message from client
            };

            new_client->on_socket_closed = [new_client](int error_code) {
                // on close, client as well server
                std::cout << "Socket closed:" << new_client->remote_address() << ":" << new_client->remote_port() << " -> " << error_code << std::endl;
                // TODO ?
            };
        };

        // bind the server to a port.
        tcp_server.bind_port(7777, [](int error_code, const std::string &error_message) {
            // BINDING FAILED:
            std::cout << error_code << " : " << error_message << std::endl;
            // TODO: binding wasn't possible, therefore exit ARTSS because something's wrong with the environment
        });

        // Start Listening the server.
        tcp_server.start_listening([](int error_code, const std::string &error_message) {
            // LISTENING FAILED:
            std::cout << error_code << " : " << error_message << std::endl;
            // TODO: listening not possible, therefore exit ARTSS because TCP Server cannot be used
        });

        // TODO: You should do an input loop so the program will not terminated immediately: loop until main thread of MPI says to exit
        std::string input;
        getline(std::cin, input);
        while (input != "exit") {
            getline(std::cin, input);
        }

        // close the server before exiting the program.
        tcp_server.close_socket();
    }

}
#endif