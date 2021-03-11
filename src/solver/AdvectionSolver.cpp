/// \file       AdvectionSolver.cpp
/// \brief      Defines the steps to solve the advection equation
/// \date       August 22, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "AdvectionSolver.h"
#include "../interfaces/IAdvection.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"

AdvectionSolver::AdvectionSolver(FieldController *field_controlller) :
    m_u_lin(FieldType::U,
                Parameters::getInstance()->get_real("initial_conditions/w_lin"),
                0,
                Domain::getInstance()->get_size()),
    m_v_lin(FieldType::V,
                Parameters::getInstance()->get_real("initial_conditions/w_lin"),
                0,
                Domain::getInstance()->get_size()),
    m_w_lin(FieldType::W,
                Parameters::getInstance()->get_real("initial_conditions/w_lin"),
                0,
                Domain::getInstance()->get_size()) {
#ifndef BENCHMARKING
     m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = field_controlller;

    auto params = Parameters::getInstance();
    std::string advectionType = params->get("solver/advection/type");
    SolverSelection::SetAdvectionSolver(&adv, params->get("solver/advection/type"));

    m_u_lin.copyin();
    m_v_lin.copyin();
    m_w_lin.copyin();

    control();
}

AdvectionSolver::~AdvectionSolver() {
    delete adv;
}

//====================================== do_step =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronous kernel launching (true, default: false)
// ***************************************************************************************
void AdvectionSolver::do_step(real, bool sync) {
  // local variables and parameters for GPU
    Field &u = m_field_controller->field_u;
    Field &v = m_field_controller->field_v;
    Field &w = m_field_controller->field_w;
    Field &u0 = m_field_controller->field_u0;
    Field &v0 = m_field_controller->field_v0;
    Field &w0 = m_field_controller->field_w0;

#pragma acc data present(u_lin, v_lin, w_lin, u, u0, v, v0, w, w0)
    {
// 1. Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect ...");
#endif
        adv->advect(u, u0, m_u_lin, m_v_lin, m_w_lin, sync);
        adv->advect(v, v0, m_u_lin, m_v_lin, m_w_lin, sync);
        adv->advect(w, w0, m_u_lin, m_v_lin, m_w_lin, sync);
    }
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void AdvectionSolver::control() {
    auto params = Parameters::getInstance();
    if (params->get("solver/advection/field") != "u,v,w") {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(typeid(AdvectionSolver).name());
        logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
}
