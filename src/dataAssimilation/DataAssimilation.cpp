/// \file      DataAssimilation.cpp
/// \brief
/// \date      May 18, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include "DataAssimilation.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../TCP/TCPServer.h"
#include <mpi.h>

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
        config_MPI();
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

void DataAssimilation::config_MPI() {
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(NULL, NULL);

    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);

    if (comm_rank == 1) {
        TCPServer tcp_server;
        tcp_server.onNewConnection = [&](TCPSocket *new_client) {
            std::cout << "New client: [";
            std::cout << new_client->remoteAddress() << ":" << new_client->remotePort() << "]" << std::endl;

            new_client->onMessageReceived = [new_client](std::string message) {
                std::cout << new_client->remoteAddress() << ":" << new_client->remotePort() << " => " << message << std::endl;
                new_client->Send("OK!");
            };

            new_client->onSocketClosed = [new_client](int errorCode) {
                std::cout << "Socket closed:" << new_client->remoteAddress() << ":" << new_client->remotePort() << " -> " << errorCode << std::endl;
            };
        };

        // Bind the server to a port.
        tcp_server.Bind(8888, [](int errorCode, std::string errorMessage) {
            // BINDING FAILED:
            std::cout << errorCode << " : " << errorMessage << std::endl;
        });

        // Start Listening the server.
        tcp_server.Listen([](int errorCode, std::string errorMessage) {
            // LISTENING FAILED:
            std::cout << errorCode << " : " << errorMessage << std::endl;
        });

        // You should do an input loop so the program will not terminated immediately:
        std::string input;
        getline(std::cin, input);
        while (input != "exit") {
            getline(std::cin, input);
        }

        // Close the server before exiting the program.
        tcp_server.Close();
    }

}
