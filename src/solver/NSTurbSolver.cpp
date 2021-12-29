/// \file       NSTurbSolver.cpp
/// \brief      Defines the steps to solve advection, diffusion, pressure and add sources (with LES turbulence)
/// \date       Feb 15, 2017
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "NSTurbSolver.h"

#include <string>
#include <vector>
#include <algorithm>

#include "SolverSelection.h"
#include "../domain/DomainData.h"

NSTurbSolver::NSTurbSolver(const Settings::solver_parameters &solver_settings, Settings::Settings const &settings, FieldController *field_controller) :
        m_settings(settings),
        m_solver_settings(solver_settings),
        m_field_controller(field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    // Advection of velocity
    SolverSelection::set_advection_solver(m_solver_settings.advection, &adv_vel);

    // Diffusion of velocity
    SolverSelection::set_diffusion_solver(m_solver_settings.diffusion, &dif_vel);

    // Turbulent viscosity
    SolverSelection::SetTurbulenceSolver(m_settings, &mu_tub, m_settings.get("solver/turbulence/type"));

    //Pressure
    SolverSelection::set_pressure_solver(m_solver_settings.pressure, &pres);

    // Source
    SolverSelection::set_source_solver(m_solver_settings.source.type, &sou_vel, m_solver_settings.source.direction);

    m_force_function = m_settings.get("solver/source/force_fct");
    control();
}

NSTurbSolver::~NSTurbSolver() {
    delete adv_vel;
    delete dif_vel;
    delete mu_tub;
    delete pres;
    delete sou_vel;
}

//=========================================== do_step ====================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void NSTurbSolver::do_step(real t, bool sync) {
    Field &u = m_field_controller->get_field_u();
    Field &v = m_field_controller->get_field_v();
    Field &w = m_field_controller->get_field_w();
    Field &u0 = m_field_controller->get_field_u0();
    Field &v0 = m_field_controller->get_field_v0();
    Field &w0 = m_field_controller->get_field_w0();
    Field &u_tmp = m_field_controller->get_field_u_tmp();
    Field &v_tmp = m_field_controller->get_field_v_tmp();
    Field &w_tmp = m_field_controller->get_field_w_tmp();
    Field &p = m_field_controller->get_field_p();
    Field &rhs = m_field_controller->get_field_rhs();
    Field &f_x = m_field_controller->get_field_force_x();
    Field &f_y = m_field_controller->get_field_force_y();
    Field &f_z = m_field_controller->get_field_force_z();
    Field &nu_t = m_field_controller->get_field_nu_t();  // nu_t - Eddy Viscosity

    auto nu = m_settings.get_real("physical_parameters/nu");

#pragma acc data present(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, \
                         p, rhs, f_x, f_y, f_z, nu_t)
    {
        // 1. Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect ...");
#endif
        adv_vel->advect(u, u0, u0, v0, w0, sync);
        adv_vel->advect(v, v0, u0, v0, w0, sync);
        adv_vel->advect(w, w0, u0, v0, w0, sync);

        // Couple velocity to prepare for diffusion
        FieldController::couple_vector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);

        // 2. Solve turbulent diffusion equation
#ifndef BENCHMARKING
        m_logger->info("Calculating Turbulent viscosity ...");
#endif
        mu_tub->calc_turbulent_viscosity(nu_t, u, v, w, true);

#ifndef BENCHMARKING
        m_logger->info("Diffuse ...");
#endif
        dif_vel->diffuse(u, u0, u_tmp, nu, nu_t, sync);
        dif_vel->diffuse(v, v0, v_tmp, nu, nu_t, sync);
        dif_vel->diffuse(w, w0, w_tmp, nu, nu_t, sync);

        // Couple data to prepare for adding source
        FieldController::couple_vector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);

        // 3. Add force
        if (m_force_function != SourceMethods::Zero) {
#ifndef BENCHMARKING
            m_logger->info("Add source ...");
#endif
            sou_vel->add_source(u, v, w, f_x, f_y, f_z, sync);

            // Couple data
            FieldController::couple_vector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);
        }

        // 4. Solve pressure equation and project
        // Calculate divergence of u
        pres->divergence(rhs, u_tmp, v_tmp, w_tmp, sync);

        // Solve pressure equation
#ifndef BENCHMARKING
        m_logger->info("Pressure ...");
#endif
        pres->pressure(p, rhs, t, sync);

        // Correct
        pres->projection(u, v, w, u_tmp, v_tmp, w_tmp, p, sync);

        // 5. Sources updated in Solver::update_sources, TimeIntegration

        if (sync) {
#pragma acc wait
        }
    }
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void NSTurbSolver::control() {
    auto fields = Utility::split(m_settings.get("solver/advection/field"), ',');
    std::sort(fields.begin(), fields.end());
    if (fields != std::vector<std::string>({"u", "v", "w"})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }

    auto diff_fields = Utility::split(m_settings.get("solver/diffusion/field"), ',');
    std::sort(diff_fields.begin(), diff_fields.end());
    if (diff_fields != std::vector<std::string>({"u", "v", "w"})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }

    if (m_settings.get("solver/pressure/field") != Mapping::get_field_type_name(FieldType::P)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
}
