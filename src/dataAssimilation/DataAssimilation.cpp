/// \file      DataAssimilation.cpp
/// \brief
/// \date      May 18, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include "DataAssimilation.h"
#include "../utility/Parameters.h"
#include "../Domain.h"

DataAssimilation::DataAssimilation(const FieldController &field_controller) : m_field_controller(field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_IO = new FieldIO(m_field_controller);

    Parameters *params = Parameters::getInstance();
    m_assimilated = (params->get("data_assimilation/enabled") == "Yes");
    if (m_assimilated) {
        std::string init = params->get("data_assimilation/class/name");
        if (init == AssimilationMethods::Test) {
#ifndef BENCHMARKING
            m_logger->critical("worked");
#endif
        } else {
#ifndef BENCHMARKING
            m_logger->critical("Data Assimilation class {} is not defined", init);
#endif
            std::exit(1);
            // TODO Error Handling
        }
    }
}

void DataAssimilation::assimilate(real t_old) {
    Field *u_current = new Field(FieldType::U);
    Field *v_current = new Field(FieldType::V);
    Field *w_current = new Field(FieldType::W);
    Field *p_current = new Field(FieldType::P);
    Field *T_current = new Field(FieldType::T);
    Field *C_current = new Field(FieldType::RHO);
    m_field_IO->read(t_old, u_current, v_current, w_current, p_current, T_current, C_current);
    if (func->control(u_current, v_current, w_current, p_current, T_current, C_current)) {
        // TODO stop simulation ?

        func->assimilate(u_current, v_current, w_current, p_current, T_current, C_current);
        m_field_controller.replace_data(u_current, v_current, w_current, p_current, T_current, C_current);

        m_t_cur = t_old;
        m_rollback = true;
    }
    delete u_current;
    delete v_current;
    delete w_current;
    delete p_current;
    delete T_current;
    delete C_current;
}

void DataAssimilation::save_data(real t_cur) {
    if (m_assimilated) {
        m_field_IO->write_out(t_cur);
    }
}

real DataAssimilation::get_new_time_value() {
    return m_t_cur;
}
