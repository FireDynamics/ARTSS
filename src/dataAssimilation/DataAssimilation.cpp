/// \file      DataAssimilation.cpp
/// \brief
/// \date      May 18, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include "DataAssimilation.h"
#include "../domain/DomainData.h"
#include "../TCP/TCPServer.h"
#include "mpi.h"
#include "TemperatureSourceChanger.h"

DataAssimilation::DataAssimilation(const SolverController &solver_controller,
                                   FieldController *field_controller,
                                   const Settings::Settings &settings) :
        m_settings(settings),
        m_field_controller(field_controller),
        m_solver_controller(solver_controller),
        m_new_field_u(Field(FieldType::U)),
        m_new_field_v(Field(FieldType::V)),
        m_new_field_w(Field(FieldType::W)),
        m_new_field_p(Field(FieldType::P)),
        m_new_field_T(Field(FieldType::T)),
        m_new_field_C(Field(FieldType::RHO)) {
    m_logger = Utility::create_logger(typeid(this).name());
    m_field_IO_handler = new FieldIO(settings.filename);
    if (m_settings.assimilation_parameters.enabled) {
        if (m_settings.assimilation_parameters.class_name == AssimilationMethods::standard) {
            m_parameter_handler = new ParameterReader();
        } else if (m_settings.assimilation_parameters.class_name == AssimilationMethods::temperature_source) {
            m_parameter_handler = new TemperatureSourceChanger(m_solver_controller,
                                                               m_settings.solver_parameters.temperature.source);
        } else {
            m_logger->error("assimilation method {} not known", m_settings.assimilation_parameters.class_name);
            std::exit(1);
        }
    }
}

void DataAssimilation::initiate_rollback() {
    m_field_controller->replace_data(m_new_field_u, m_new_field_v, m_new_field_w,
                                     m_new_field_p, m_new_field_T, m_new_field_C);
}

void DataAssimilation::save_data(real t_cur) {
    Field &u = m_field_controller->get_field_u();
    Field &v = m_field_controller->get_field_v();
    Field &w = m_field_controller->get_field_w();
    Field &p = m_field_controller->get_field_w();
    Field &T = m_field_controller->get_field_T();
    Field &C = m_field_controller->get_field_concentration();

    m_field_IO_handler->write_fields(t_cur, u, v, w, p, T, C);
}

real DataAssimilation::get_new_time_value() const {
    return m_t_cur;
}

bool DataAssimilation::config_rollback(const char *msg) {
    std::vector<std::string> divided_string = Utility::split(msg, ',');
    auto new_time = std::stod(divided_string[0]);
    m_logger->debug("current time step {}, new time {}", m_t_cur, new_time);
    if (m_t_cur < new_time) {
        m_logger->warn("simulation is currently at {}. Cannot rollback to {}", m_t_cur, new_time);
        return false;
    } else if (std::fabs(m_t_cur - new_time) < 1e-10) {  // current time step, no need to load original data
        m_t_cur = new_time;
        m_logger->info("set new time value to {}", m_t_cur);
        m_logger->debug("read config data from {}", divided_string[1]);
        auto[changes, field_changes] = m_parameter_handler->read_config(divided_string[1]);
        if (changes && field_changes.changed) {
            m_field_IO_handler->read_fields(field_changes,
                                            m_new_field_u, m_new_field_v, m_new_field_w,
                                            m_new_field_p, m_new_field_T, m_new_field_C);
            return true;
        }
        return changes;
    } else {
        m_t_cur = new_time;
        m_logger->info("set new time value to {}", m_t_cur);
        m_logger->debug("read config data from {}", divided_string[1]);
        auto[changes, field_changes] = m_parameter_handler->read_config(divided_string[1]);
        if (changes && field_changes.changed) {
            m_logger->debug("read field data from {}", field_changes.file_name);
            m_field_IO_handler->read_fields(m_t_cur, field_changes,
                                            m_new_field_u, m_new_field_v, m_new_field_w,
                                            m_new_field_p, m_new_field_T, m_new_field_C);
            return true;
        }
        return changes;
    }
}

bool DataAssimilation::requires_rollback(const real t_cur) {
    if (!m_settings.assimilation_parameters.enabled) {
        return false;
    }
    m_t_cur = t_cur;
    MPI_Status status;
    int flag = -1;
    MPI_Iprobe(1, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    m_logger->debug("probe: {}", flag);
    if (flag) {
        int msg_len;
        MPI_Get_count(&status, MPI_CHAR, &msg_len);
        std::vector<char> msg;

        msg.resize(msg_len);
        m_logger->debug("preparing to receive message");
        MPI_Recv(msg.data(), msg_len, MPI_CHAR, 1, status.MPI_TAG, MPI_COMM_WORLD, &status);
        m_logger->debug("received message: {}", msg.data());
        return config_rollback(msg.data());
    }
    return flag;
}
