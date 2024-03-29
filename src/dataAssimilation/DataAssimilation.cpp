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
#include "ObstacleChanger.h"
#include "../domain/DomainController.h"

DataAssimilation::DataAssimilation(const SolverController &solver_controller,
                                   FieldController *field_controller,
                                   const Settings::Settings &settings) :
        m_logger(Utility::create_logger(typeid(this).name())),
        m_settings(settings),
        m_field_controller(field_controller),
        m_solver_controller(solver_controller),
        m_new_field_u(Field(FieldType::U)),
        m_new_field_v(Field(FieldType::V)),
        m_new_field_w(Field(FieldType::W)),
        m_new_field_p(Field(FieldType::P)),
        m_new_field_T(Field(FieldType::T)),
        m_new_field_C(Field(FieldType::RHO)) {
    m_time_interval_counter = 0;
    m_output_time_interval = m_settings.assimilation_parameters.output_time_interval;

    m_field_IO_handler = new FieldIO(settings.filename,
                                     m_settings.assimilation_parameters.output_dir);

    if (m_settings.assimilation_parameters.enabled) {
        if (m_settings.assimilation_parameters.class_name == AssimilationMethods::standard) {
            m_parameter_handler = new ParameterReader();
        } else if (m_settings.assimilation_parameters.class_name == AssimilationMethods::temperature_source) {
            m_parameter_handler = new TemperatureSourceChanger(m_solver_controller,
                                                               m_settings.solver_parameters.temperature.source);
        } else if (m_settings.assimilation_parameters.class_name == AssimilationMethods::obstacle_changer) {
            m_parameter_handler = new ObstacleChanger(m_solver_controller, m_settings.obstacles_parameters);
        } else {
            m_logger->error("assimilation method {} not known", m_settings.assimilation_parameters.class_name);
            std::exit(1);
        }
    }
}

void DataAssimilation::initiate_time_skip(const real t_cur) {
    m_time_interval_counter = static_cast<int>(t_cur / m_output_time_interval) + 1;
    m_logger->debug("change counter to {} for {}", m_time_interval_counter, t_cur);
    m_field_controller->replace_data(m_new_field_u, m_new_field_v, m_new_field_w,
                                     m_new_field_p, m_new_field_T, m_new_field_C);
}

void DataAssimilation::save_data(real t_cur) {
    if (t_cur >= m_output_time_interval * m_time_interval_counter) {
        m_logger->debug("save data for {} with interval {} counter {}",
                        t_cur, m_output_time_interval,
                        m_time_interval_counter);
        m_time_interval_counter++;
        Field &u = m_field_controller->get_field_u();
        Field &v = m_field_controller->get_field_v();
        Field &w = m_field_controller->get_field_w();
        Field &p = m_field_controller->get_field_p();
        Field &T = m_field_controller->get_field_T();
        Field &C = m_field_controller->get_field_concentration();

        m_field_IO_handler->write_fields(t_cur, u, v, w, p, T, C);
    }
    m_field_IO_handler->create_meta_file(t_cur);
}

real DataAssimilation::get_new_time_value() const {
    return m_t_cur;
}

bool DataAssimilation::config_rollback(const char *msg) {
    const DataAssimilationPackageHeader package = *reinterpret_cast<const DataAssimilationPackageHeader *>(msg);
    const std::string file_name(msg + 12, package.file_name_len);
    m_logger->info("current time step {}, new time {}", m_t_cur, package.time);
    m_logger->info("new config file {}", file_name);
    if (m_t_cur < package.time) {
        m_logger->warn("simulation is currently at {}. Cannot rollback to {}", m_t_cur, package.time);
        return false;
    } else if (std::fabs(m_t_cur - package.time) < 1e-10) {  // current time step, no need to load original data
        m_t_cur = package.time;
        m_logger->info("already at time step {}", m_t_cur);
        m_logger->debug("read config data from {}", file_name);
        auto [changes, field_changes] = m_parameter_handler->read_config(file_name);
        if (field_changes.changed) {
            m_logger->debug("read field data from {}", field_changes.file_name);
            m_field_IO_handler->read_changed_fields(field_changes,
                                                    m_new_field_u, m_new_field_v, m_new_field_w,
                                                    m_new_field_p, m_new_field_T, m_new_field_C);
            return true;
        }
        return changes;
    } else {
        m_t_cur = package.time;
        m_logger->info("set new time value to {}", m_t_cur);
        m_logger->debug("read config data from {}", file_name);
        auto [changes, field_changes] = m_parameter_handler->read_config(file_name);
        m_logger->debug("read field data from {}", field_changes.file_name);
        m_logger->debug("read original field data @t: {}", m_t_cur);
        m_field_IO_handler->read_fields(m_t_cur, field_changes,
                                        m_new_field_u, m_new_field_v, m_new_field_w,
                                        m_new_field_p, m_new_field_T, m_new_field_C);
        auto domain_controller = DomainController::getInstance();
        if (field_changes.T_changed) {
            domain_controller->apply_boundary(m_new_field_T);
        }
        //TODO
        return changes || field_changes.changed;
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
        m_logger->debug("received chars1: {}", msg_len);
        MPI_Get_count(&status, MPI_CHAR, &msg_len);
        std::vector<char> msg;

        msg.resize(msg_len);
        m_logger->debug("preparing to receive message");
        MPI_Recv(msg.data(), msg_len, MPI_CHAR, 1, status.MPI_TAG, MPI_COMM_WORLD, &status);
        m_logger->debug("received chars: {}", msg_len);
        bool result = config_rollback(msg.data());
        std::string response = "Rollback done";
        MPI_Request request;
        m_logger->debug("send response: {}", response);
        MPI_Isend(response.c_str(), response.size() + 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &request);
        m_logger->debug("Isend done", msg_len);
        return result;
    }
    return flag;
}

bool DataAssimilation::load_data() {
    bool load_data = m_settings.assimilation_parameters.load_data;
    if (!load_data) {
        return load_data;
    }
    m_t_cur = m_settings.assimilation_parameters.time;

    m_field_IO_handler->read_fields(m_settings.assimilation_parameters.file,
                                    m_t_cur,
                                    m_new_field_u, m_new_field_v, m_new_field_w,
                                    m_new_field_p, m_new_field_T, m_new_field_C);

    return load_data;
}
