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

AdvectionSolver::AdvectionSolver(FieldController *field_controlller) {
#ifndef BENCHMARKING
     m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = field_controlller;

    auto params = Parameters::getInstance();
    std::string advectionType = params->get("solver/advection/type");
    SolverSelection::SetAdvectionSolver(&adv, params->get("solver/advection/type"));

    auto domain = Domain::getInstance();

    real d_u_linm = params->get_real("initial_conditions/u_lin");
    real d_v_linm = params->get_real("initial_conditions/v_lin");
    real d_w_linm = params->get_real("initial_conditions/w_lin");

    u_linm = new Field(FieldType::U, d_u_linm, 0, domain->get_size());
    v_linm = new Field(FieldType::V, d_v_linm, 0, domain->get_size());
    w_linm = new Field(FieldType::W, d_w_linm, 0, domain->get_size());

    auto u_lin = u_linm;
    auto v_lin = v_linm;
    auto w_lin = w_linm;

    u_lin->copyin();
    v_lin->copyin();
    w_lin->copyin();

    control();
}

AdvectionSolver::~AdvectionSolver() {
    delete adv;
    delete u_linm;
    delete v_linm;
    delete w_linm;
}

//====================================== do_step =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronous kernel launching (true, default: false)
// ***************************************************************************************
void AdvectionSolver::do_step(real, bool sync) {
  // local variables and parameters for GPU
    auto u = m_field_controller->field_u;
    auto v = m_field_controller->field_v;
    auto w = m_field_controller->field_w;
    auto u0 = m_field_controller->field_u0;
    auto v0 = m_field_controller->field_v0;
    auto w0 = m_field_controller->field_w0;

    auto u_lin = u_linm;
    auto v_lin = v_linm;
    auto w_lin = w_linm;

#pragma acc data present(u_lin, v_lin, w_lin, u, u0, v, v0, w, w0)
    {
// 1. Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect ...");
#endif
        adv->advect(u, u0, u_lin, v_lin, w_lin, sync);
        adv->advect(v, v0, u_lin, v_lin, w_lin, sync);
        adv->advect(w, w0, u_lin, v_lin, w_lin, sync);
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
