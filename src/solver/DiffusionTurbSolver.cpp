/// \file       DiffusionTurbSolver.h
/// \brief      Defines the steps to solve the turbulent diffusion equation
/// \date       Aug 18, 2016
/// \author     Suryanarayana Maddu
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <vector>
#include <algorithm>

#include "DiffusionTurbSolver.h"


DiffusionTurbSolver::DiffusionTurbSolver(const Settings::solver_parameters &solver_settings, FieldController *field_controller) :
        m_solver_settings(solver_settings), m_field_controller(field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    //Diffusion
    SolverSelection::set_diffusion_solver(m_solver_settings.diffusion, &dif);

    // Turbulent viscosity
    SolverSelection::set_turbulence_solver(m_solver_settings.turbulence, &mu_tub);
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

    auto domain_data = DomainData::getInstance();
    real nu = domain_data->get_physical_parameters().nu.value();

#pragma acc data present(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, nu_t)
    {
#ifndef BENCHMARKING
        m_logger->info("Calculating Turbulent viscosity ...");
#endif
        mu_tub->calc_turbulent_viscosity(nu_t, u, v, w, true);
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
    auto diff_fields = m_solver_settings.diffusion.fields;
    std::sort(diff_fields.begin(), diff_fields.end());
    if (diff_fields != std::vector<FieldType>({FieldType::U, FieldType::V, FieldType::W})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO(issue 6) Error handling
    }
}
