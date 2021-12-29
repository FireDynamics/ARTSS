/// \file       NSTempTurbConSolver.cpp
/// \brief      Defines the (fractional) steps to solve the incompressible Navier-Stokes equations with force f(T), turbulence and concentration C
/// \date       Oct 02, 2017
/// \author     Küsters
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "NSTempTurbConSolver.h"

#include <string>
#include <vector>
#include <algorithm>

#include "SolverSelection.h"
#include "../domain/DomainData.h"

NSTempTurbConSolver::NSTempTurbConSolver(const Settings::solver_parameters &solver_settings, Settings::Settings const &settings, FieldController *field_controller) :
        m_settings(settings),
        m_field_controller(field_controller),
        m_solver_settings(solver_settings) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    // Advection of velocity
    SolverSelection::set_advection_solver(m_solver_settings.advection, &adv_vel);

    // Advection of temperature
    SolverSelection::set_advection_solver(m_solver_settings.temperature.advection, &adv_temp);

    // Advection of concentration
    SolverSelection::set_advection_solver(m_solver_settings.concentration.advection, &adv_temp);

    // Diffusion of velocity
    SolverSelection::SetDiffusionSolver(m_settings, &dif_vel, m_settings.get("solver/diffusion/type"));

    // Turbulent viscosity for velocity diffusion
    SolverSelection::SetTurbulenceSolver(m_settings, &mu_tub, m_settings.get("solver/turbulence/type"));

    // Diffusion of temperature
    SolverSelection::SetDiffusionSolver(m_settings, &dif_temp, m_settings.get("solver/temperature/diffusion/type"));

    // Diffusion for concentration
    SolverSelection::SetDiffusionSolver(m_settings, &dif_con, m_settings.get("solver/concentration/diffusion/type"));

    // Pressure
    SolverSelection::SetPressureSolver(m_solver_settings.pressure, &pres);

    // Source of velocity
    SolverSelection::set_source_solver(m_solver_settings.source.type, &sou_vel, m_solver_settings.source.direction);

    // Source of temperature
    SolverSelection::set_source_solver(m_solver_settings.temperature.source.type, &sou_temp, m_solver_settings.temperature.source.dir);

    // Source of concentration
    SolverSelection::set_source_solver(m_solver_settings.concentration.source.type, &sou_con, m_solver_settings.concentration.source.dir);

    // Constants
    m_dir_vel = m_settings.get("solver/source/dir");
    m_has_turbulence_temperature = m_settings.get_bool("solver/temperature/turbulence/include");
    m_has_turbulence_concentration = m_settings.get_bool("solver/concentration/turbulence/include");
    m_has_dissipation = m_settings.get_bool("solver/temperature/source/dissipation");
    m_forceFct = m_settings.get("solver/source/force_fct");
    m_tempFct = m_settings.get("solver/temperature/source/temp_fct");
    m_conFct = m_settings.get("solver/concentration/source/con_fct");
    SolverSelection::set_temperature_source_function(m_solver_settings.temperature.source, &m_source_function_temperature);
    SolverSelection::set_concentration_source_function(m_solver_settings.concentration.source, &m_source_function_concentration);
    control();
}

NSTempTurbConSolver::~NSTempTurbConSolver() {
    delete adv_vel;
    delete dif_vel;
    delete adv_temp;
    delete dif_temp;
    delete adv_con;
    delete dif_con;
    delete mu_tub;
    delete pres;
    delete sou_vel;
    delete sou_temp;
    delete sou_con;
}

//====================================== do_step =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt          time step
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void NSTempTurbConSolver::do_step(real t, bool sync) {
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
    Field &T = m_field_controller->get_field_T();
    Field &T0 = m_field_controller->get_field_T0();
    Field &T_tmp = m_field_controller->get_field_T_tmp();
    Field &C = m_field_controller->get_field_concentration();
    Field &C0 = m_field_controller->get_field_concentration0();
    Field &C_tmp = m_field_controller->get_field_concentration_tmp();
    Field &f_x = m_field_controller->get_field_force_x();
    Field &f_y = m_field_controller->get_field_force_y();
    Field &f_z = m_field_controller->get_field_force_z();
    Field &S_T = m_field_controller->get_field_source_T();
    Field &S_C = m_field_controller->get_field_source_concentration();
    Field &nu_t = m_field_controller->get_field_nu_t();      // nu_t - Eddy Viscosity
    Field &kappa_t = m_field_controller->get_field_kappa();  // kappa_t - Eddy thermal diffusivity
    Field &gamma_t = m_field_controller->get_field_gamma();  // gamma_t - Eddy mass diffsusivity

#pragma acc data present(u, u0, u_tmp, v, v0, v_tmp, w, \
                         w0, w_tmp, p, rhs, T, T0, T_tmp, \
                         C, C0, C_tmp, f_x, f_y, f_z, S_T, S_C, \
                         nu_t, kappa_t, gamma_t)
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
        const real nu = m_settings.get_real("physical_parameters/nu");

        dif_vel->diffuse(u, u0, u_tmp, nu, nu_t, sync);
        dif_vel->diffuse(v, v0, v_tmp, nu, nu_t, sync);
        dif_vel->diffuse(w, w0, w_tmp, nu, nu_t, sync);

        // Couple data to prepare for adding source
        FieldController::couple_vector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);

// 3. Add force
        if (m_forceFct != SourceMethods::Zero) {
#ifndef BENCHMARKING
            m_logger->info("Add momentum source ...");
#endif
            sou_vel->add_source(u, v, w, f_x, f_y, f_z, sync);

            // Couple data to prepare for adding source
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

// 5. Solve Temperature and link back to force

        // Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect Temperature ...");
#endif
        adv_temp->advect(T, T0, u, v, w, sync);

        // Couple temperature to prepare for diffusion
        FieldController::couple_scalar(T, T0, T_tmp, sync);

        // Solve diffusion equation
        // turbulence
        real kappa = m_settings.get_real("physical_parameters/kappa");
        if (m_has_turbulence_temperature) {
            real Pr_T = m_settings.get_real("solver/temperature/turbulence/Pr_T");
            real rPr_T = 1. / Pr_T;

            size_t bsize = kappa_t.get_size();
#pragma acc parallel loop independent present(kappa_t, nu_t) async
            for (size_t i = 0; i < bsize; ++i) {
                kappa_t[i] = nu_t[i] * rPr_T;  // kappa_turb = nu_turb/Pr_turb
            }

#ifndef BENCHMARKING
            m_logger->info("Diffuse turbulent Temperature ...");
#endif
            dif_temp->diffuse(T, T0, T_tmp, kappa, kappa_t, sync);

            // Couple temperature to prepare for adding source
            FieldController::couple_scalar(T, T0, T_tmp, sync);
        } else {
            // no turbulence
            if (kappa != 0.) {

#ifndef BENCHMARKING
                m_logger->info("Diffuse Temperature ...");
#endif
                dif_temp->diffuse(T, T0, T_tmp, kappa, sync);

                // Couple temperature to prepare for adding source
                FieldController::couple_scalar(T, T0, T_tmp, sync);
            }
        }

        // Add dissipation
        if (m_has_dissipation) {

#ifndef BENCHMARKING
            m_logger->info("Add dissipation ...");
#endif
            sou_temp->dissipate(T, u, v, w, sync);

            // Couple temperature
            FieldController::couple_scalar(T, T0, T_tmp, sync);
        } else if (m_tempFct != SourceMethods::Zero) {  // Add source
#ifndef BENCHMARKING
            m_logger->info("Add temperature source ...");
#endif
            sou_temp->add_source(T, S_T, sync);

            // Couple temperature
            FieldController::couple_scalar(T, T0, T_tmp, sync);
        }

// 6. Solve Concentration

        // Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect Concentration ...");
#endif
        adv_con->advect(C, C0, u, v, w, sync);

        // Couple concentration to prepare for diffusion
        FieldController::couple_scalar(C, C0, C_tmp, sync);

        // Solve diffusion equation
        // turbulence
        real gamma = m_settings.get_real("solver/concentration/diffusion/gamma");
        if (m_has_turbulence_concentration) {
            real Sc_T = m_settings.get_real("solver/concentration/turbulence/Sc_T");
            real rSc_T = 1. / Sc_T;

            size_t bsize = gamma_t.get_size();
#pragma acc parallel loop independent present(gamma_t, nu_t) async
            for (size_t i = 0; i < bsize; ++i) {
                gamma_t[i] = nu_t[i] * rSc_T;  // gamma_turb = nu_turb/Sc_turb
            }

#ifndef BENCHMARKING
            m_logger->info("Diffuse turbulent Concentration ...");
#endif
            dif_con->diffuse(C, C0, C_tmp, gamma, gamma_t, sync);

            // Couple concentration to prepare for adding source
            FieldController::couple_scalar(C, C0, C_tmp, sync);
        } else {
            // no turbulence
            if (gamma != 0.) {
#ifndef BENCHMARKING
                m_logger->info("Diffuse Concentration ...");
#endif
                dif_con->diffuse(C, C0, C_tmp, gamma, sync);

                // Couple concentration to prepare for adding source
                FieldController::couple_scalar(C, C0, C_tmp, sync);
            }
        }

        // Add source
        if (m_conFct != SourceMethods::Zero) {
#ifndef BENCHMARKING
            m_logger->info("Add concentration source ...");
#endif
            sou_con->add_source(C, S_C, sync);

            // Couple concentration
            FieldController::couple_scalar(C, C0, C_tmp, sync);
        }

// 7. Sources updated in Solver::update_sources, TimeIntegration

        if (sync) {
#pragma acc wait
        }
    }
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void NSTempTurbConSolver::control() {
    auto adv_fields = m_solver_settings.advection.fields;
    std::sort(adv_fields.begin(), adv_fields.end());
    if (adv_fields != std::vector<FieldType>({FieldType::U, FieldType::V, FieldType::W})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }

    auto diff_fields = m_solver_settings.diffusion.fields;
    std::sort(diff_fields.begin(), diff_fields.end());
    if (adv_fields != std::vector<FieldType>({FieldType::U, FieldType::V, FieldType::W})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }

    if (m_solver_settings.temperature.advection.fields != std::vector<FieldType>{FieldType::T}) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (m_solver_settings.temperature.diffusion.fields != std::vector<FieldType>{FieldType::T}) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (m_solver_settings.concentration.advection.fields != std::vector<FieldType>{FieldType::RHO}) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (m_solver_settings.concentration.diffusion.fields != std::vector<FieldType>{FieldType::RHO}) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (m_solver_settings.pressure.field != FieldType::P) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
}

void NSTempTurbConSolver::update_source(real t_cur) {
    m_source_function_temperature->update_source(m_field_controller->get_field_source_T(), t_cur);
    m_source_function_concentration->update_source(m_field_controller->get_field_source_concentration(), t_cur);
}
