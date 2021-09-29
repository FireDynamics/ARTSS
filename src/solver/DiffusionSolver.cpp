/// \file       DiffusionSolver.cpp
/// \brief      Defines the steps to solve the diffusion equation
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "DiffusionSolver.h"
#include "../interfaces/IDiffusion.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"

DiffusionSolver::DiffusionSolver(FieldController *field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = field_controller;
    auto params = Parameters::getInstance();
    std::string diffusionType = params->get("solver/diffusion/type");
    SolverSelection::SetDiffusionSolver(&this->dif, diffusionType);

    m_nu = params->get_real("physical_parameters/nu");
    control();
}

DiffusionSolver::~DiffusionSolver() {
    delete dif;
}

//====================================== do_step =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronous kernel launching (true, default: false)
// ***************************************************************************************
void DiffusionSolver::do_step(real, bool sync) {
#ifndef BENCHMARKING
    m_logger->info("Diffuse ...");
#endif
    // 1. Solve diffusion equation
    Field &u = m_field_controller->get_field_u();
    Field &v = m_field_controller->get_field_v();
    Field &w = m_field_controller->get_field_w();
    Field &u0 = m_field_controller->get_field_u0();
    Field &v0 = m_field_controller->get_field_v0();
    Field &w0 = m_field_controller->get_field_w0();
    Field &u_tmp = m_field_controller->get_field_u_tmp();
    Field &v_tmp = m_field_controller->get_field_v_tmp();
    Field &w_tmp = m_field_controller->get_field_w_tmp();

#pragma acc data present(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp)
    {
        dif->diffuse(u, u0, u_tmp, m_nu, sync);
        dif->diffuse(v, v0, v_tmp, m_nu, sync);
        dif->diffuse(w, w0, w_tmp, m_nu, sync);
    }
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void DiffusionSolver::control() {
    auto params = Parameters::getInstance();
    if (params->get("solver/diffusion/field") != "u,v,w") {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(typeid(DiffusionSolver).name());
        logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
}
