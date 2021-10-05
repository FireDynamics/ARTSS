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
#include "HRRChanger.h"

DataAssimilation::DataAssimilation(const SolverController &solver_controller,
                                   const FieldController &field_controller) :
        m_field_controller(field_controller),
        m_solver_controller(solver_controller),
        m_new_field_u(Field(FieldType::U)),
        m_new_field_v(Field(FieldType::V)),
        m_new_field_w(Field(FieldType::W)),
        m_new_field_p(Field(FieldType::P)),
        m_new_field_T(Field(FieldType::T)),
        m_new_field_C(Field(FieldType::RHO)){
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    Parameters *params = Parameters::getInstance();
    m_assimilated = (params->get("data_assimilation/enabled") == "Yes");
}

void DataAssimilation::initiate_rollback() {
    m_field_controller.replace_data(m_new_field_u, m_new_field_v, m_new_field_w, m_new_field_p, m_new_field_T, m_new_field_C);
}
void DataAssimilation::read_new_data(std::string &file_name) {
#ifndef BENCHMARKING
    m_logger->debug("read new data from {}", file_name);
#endif
    m_reader->read(file_name, m_new_field_u, m_new_field_v, m_new_field_w, m_new_field_p, m_new_field_T, m_new_field_C);
}

void DataAssimilation::save_data(real t_cur) {
    Field &u = m_field_controller.get_field_u();
    Field &v = m_field_controller.get_field_v();
    Field &w = m_field_controller.get_field_w();
    Field &p = m_field_controller.get_field_w();
    Field &T = m_field_controller.get_field_T();
    Field &C = m_field_controller.get_field_concentration();

    m_reader->write(t_cur, u, v, w, p, T, C);
}

real DataAssimilation::get_new_time_value() const {
    return m_t_cur;
}

void DataAssimilation::config_rollback(const char *msg) {
    std::vector<std::string> splitted_string = Utility::split(msg, '|');
    m_t_cur = std::stod(splitted_string[0]);
#ifndef BENCHMARKING
    m_logger->debug("set new time value to {}", m_t_cur);
#endif
    read_new_data(splitted_string[1]);
}

bool DataAssimilation::requires_rollback() {
    MPI_Status status;
    int flag = -1;
    MPI_Iprobe(1, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    if (flag) {
        int msg_len;
        MPI_Get_count(&status, MPI_CHAR, &msg_len);
        char *msg = new char[msg_len];
        MPI_Recv(msg, msg_len, MPI_CHAR, 1, status.MPI_TAG, MPI_COMM_WORLD, &status);
        config_rollback(msg);  // TODO(MPI): could be done in a third process
    }
    return flag;
}
