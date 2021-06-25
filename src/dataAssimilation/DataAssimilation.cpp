/// \file      DataAssimilation.cpp
/// \brief
/// \date      May 18, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include "DataAssimilation.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../TCP/TCPServer.h"
#include "mpi.h"

bool DataAssimilation::simulation_is_running = true;

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
        m_new_field_u = new Field(FieldType::U);
        m_new_field_v = new Field(FieldType::V);
        m_new_field_w = new Field(FieldType::W);
        m_new_field_p = new Field(FieldType::P);
        m_new_field_T = new Field(FieldType::T);
        m_new_field_C = new Field(FieldType::RHO);
        config_MPI();
    }
}

void DataAssimilation::initiate_rollback() {
    m_rollback = false;
    m_field_controller.replace_data(m_new_field_u, m_new_field_v, m_new_field_w, m_new_field_p, m_new_field_T, m_new_field_C);
}
void DataAssimilation::read_new_data(std::string &file_name) {
#ifndef BENCHMARKING
    m_logger->debug("read new data from {}", file_name);
#endif
    m_t_cur = m_func->read(file_name, m_new_field_u, m_new_field_v, m_new_field_w, m_new_field_p, m_new_field_T, m_new_field_C);
}

void DataAssimilation::save_data(real t_cur) {
    real *u = m_field_controller.get_field_u_data();
    real *v = m_field_controller.get_field_v_data();
    real *w = m_field_controller.get_field_w_data();
    real *p = m_field_controller.get_field_w_data();
    real *T = m_field_controller.get_field_T_data();
    real *C = m_field_controller.get_field_concentration_data();

    m_func->write(t_cur, u, v, w, p, T, C);
}

real DataAssimilation::get_new_time_value() const {
    return m_t_cur;
}

void DataAssimilation::config_rollback(const std::string &message) {
    std::vector<std::string> splitted_string = Utility::split(message, '|');
    m_t_cur = std::stod(splitted_string[0]);
#ifndef BENCHMARKING
    m_logger->debug("set new time value to {}", m_t_cur);
#endif
    m_rollback = true;
    read_new_data(splitted_string[1]);
}

void DataAssimilation::config_MPI() {
#ifndef BENCHMARKING
    m_logger->debug("config MPI");
#endif
    MPI_Init(nullptr, nullptr);
#ifndef BENCHMARKING
    m_logger->debug("config MPI Init");
#endif

    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#ifndef BENCHMARKING
    m_logger->debug("number of processes {}", num_procs);
#endif

    int comm_size, comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    if (num_procs < 2 && comm_rank == 0) {
#ifndef BENCHMARKING
        m_logger->critical("Data Assimilation needs at least 2 processes, current number of processes: {}", num_procs);
#endif
        // TODO Error Handling
        std::exit(1);
    }

    if (comm_rank == 1) {
        TCPServer tcp_server;
        tcp_server.on_new_connection = [&](TCPSocket *new_client) {
#ifndef BENCHMARKING
            m_logger->info("New client: [{}:{}]", new_client->remote_address(), new_client->remote_port());
#endif

            new_client->on_message_received = [new_client, this](const std::string &message) {  // message from client
#ifndef BENCHMARKING
                this->m_logger->info("received message from client {}:{} => {}", new_client->remote_address(), new_client->remote_port(), message);
#endif
                new_client->send_message("OK!");  // send a message back (acknowledgment/error/whatever)
                this->config_rollback(message);
            };

            new_client->on_socket_closed = [new_client, this](int error_code) {
                // on close, client as well server
#ifndef BENCHMARKING
                this->m_logger->info("Socket closed: {}:{} -> {}", new_client->remote_address(), new_client->remote_port(), error_code);
#endif
            };
        };

        // bind the server to a port.
#ifndef BENCHMARKING
        m_logger->debug("bind server to a port");
#endif
        tcp_server.bind_port(7777, [this](int error_code, const std::string &error_message) {
            // BINDING FAILED:
#ifndef BENCHMARKING
            m_logger->critical("Binding failed - {} : {}", error_code, error_message);
#endif
            // binding wasn't possible, therefore exit ARTSS because something's wrong with the environment
            std::exit(1);
            // TODO Error Handling
        });

        // Start Listening the server.
#ifndef BENCHMARKING
        m_logger->debug("server: start listening");
#endif
        tcp_server.start_listening([this](int error_code, const std::string &error_message) {
            // LISTENING FAILED:
#ifndef BENCHMARKING
            m_logger->critical("Listening failed - {} : {}", error_code, error_message);
#endif
            // listening not possible, therefore exit ARTSS because TCP Server cannot be used
            std::exit(1);
            // TODO Error Handling
        });
#ifndef BENCHMARKING
        m_logger->debug("completed configuration of TCP server");
#endif

        // TODO: You should do an input loop so the program will not terminated immediately: loop until main thread of MPI says to exit
        while (simulation_is_running) {
        }

        // close the server before exiting the program.
        tcp_server.close_socket();
    }

}