/// \file      DataAssimilation.cpp
/// \brief
/// \date      May 18, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include "DataAssimilation.h"
#include "../utility/Parameters.h"
#include "../Domain.h"

DataAssimilation::DataAssimilation(FieldController *field_controller) : m_field_controller(field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    real dt = Parameters::getInstance()->get_real("physical_parameters/dt");
    m_field_IO = new FieldIO(m_field_controller, dt);

    Parameters *params = Parameters::getInstance();
    m_assimilated = (params->get("data_assimilation/enabled") == "Yes");
    if (m_assimilated) {
        std::string init = params->get("data_assimilation/class/name");
        if (init == AssimilationMethods::Test) {
            //func = new Layers(m_field_controller);
#ifndef BENCHMARKING
            m_logger->critical("worked");
#endif
        } else {
#ifndef BENCHMARKING
            m_logger->critical("Type {} is not defined", init);
#endif
            std::exit(1);
            // TODO Error Handling
        }
    }
}

void DataAssimilation::assimilate(real t_cur) {
    Field *u_current = new Field(FieldType::U);
    Field *v_current = new Field(FieldType::V);
    Field *w_current = new Field(FieldType::W);
    Field *p_current = new Field(FieldType::P);
    Field *T_current = new Field(FieldType::T);
    Field *C_current = new Field(FieldType::RHO);
    m_field_IO->read(t_cur, u_current, v_current, w_current, p_current, T_current, C_current);
    if (func->control()) {

    }
}
