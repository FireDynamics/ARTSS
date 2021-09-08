/// \file       DiffusionTurbSolver.h
/// \brief      Defines the steps to solve the turbulent diffusion equation
/// \date       Aug 18, 2016
/// \author     Suryanarayana Maddu
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <string>


#include "DiffusionTurbSolver.h"


DiffusionTurbSolver::DiffusionTurbSolver(FieldController *field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = field_controller;

    auto params = Parameters::getInstance();

    //Diffusion
    std::string diffusionType = params->get("solver/diffusion/type");
    SolverSelection::SetDiffusionSolver(&this->dif, diffusionType);

    m_nu = params->get_real("physical_parameters/nu");

    // Turbulent viscosity
    std::string turbluenceType = params->get("solver/turbulence/type");
    SolverSelection::SetTurbulenceSolver(&this->mu_tub, turbluenceType);
    control();
}


DiffusionTurbSolver::~DiffusionTurbSolver() {
    delete dif;
    delete mu_tub;
}

//====================================== do_step =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronous kernel launching (true, default: false)
// ***************************************************************************************
void DiffusionTurbSolver::do_step(real, bool sync) {
// 1. Solve diffusion equation
// local variables and parameters for GPU
    Field &u = m_field_controller->get_field_u();
    Field &v = m_field_controller->get_field_v();
    Field &w = m_field_controller->get_field_w();
    Field &u0 = m_field_controller->get_field_u0();
    Field &v0 = m_field_controller->get_field_v0();
    Field &w0 = m_field_controller->get_field_w0();
    Field &u_tmp = m_field_controller->get_field_u_tmp();
    Field &v_tmp = m_field_controller->get_field_v_tmp();
    Field &w_tmp = m_field_controller->get_field_w_tmp();
    Field &nu_t = m_field_controller->get_field_nu_t();  // Eddy Viscosity

    auto nu = m_nu;

#pragma acc data present(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, nu_t)
    {
#ifndef BENCHMARKING
        m_logger->info("Calculating Turbulent viscosity ...");
#endif
        mu_tub->CalcTurbViscosity(nu_t, u, v, w, true);
#ifndef BENCHMARKING
        m_logger->info("Diffuse ...");
#endif
        dif->diffuse(u, u0, u_tmp, nu, nu_t, sync);
        dif->diffuse(v, v0, v_tmp, nu, nu_t, sync);
        dif->diffuse(w, w0, w_tmp, nu, nu_t, sync);
    }
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void DiffusionTurbSolver::control() {
    auto params = Parameters::getInstance();
    if (params->get("solver/diffusion/field") != "u,v,w") {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(typeid(DiffusionTurbSolver).name());
        logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
}
