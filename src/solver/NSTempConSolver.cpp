/// \file       NSTempConSolver.cpp
/// \brief      Navier-Stokes Solver with force f(T)
/// \details    Defines the steps to solve advection, diffusion, pressure and add sources (dependent on T), solves temperature and concentration equation
/// \date       Sep 27, 2017
/// \author     Küsters
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "NSTempConSolver.h"

#include <string>
#include <vector>
#include <algorithm>

#include "../pressure/VCycleMG.h"
#include "../DomainData.h"
#include "../boundary/BoundaryData.h"
#include "SolverSelection.h"

NSTempConSolver::NSTempConSolver(Settings::Settings const &settings, FieldController *field_controller) :
        m_settings(settings) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = field_controller;

    // Advection of velocity
    SolverSelection::SetAdvectionSolver(m_settings, &adv_vel, settings.get("solver/advection/type"));

    // Advection of temperature
    SolverSelection::SetAdvectionSolver(m_settings, &adv_temp, settings.get("solver/temperature/advection/type"));

    // Advection of concentration
    SolverSelection::SetAdvectionSolver(m_settings, &adv_con, settings.get("solver/concentration/advection/type"));

    // Diffusion of velocity
    SolverSelection::SetDiffusionSolver(m_settings, &dif_vel, settings.get("solver/diffusion/type"));

    // Diffusion of temperature
    SolverSelection::SetDiffusionSolver(m_settings, &dif_temp, settings.get("solver/temperature/diffusion/type"));

    // Diffusion for concentration
    SolverSelection::SetDiffusionSolver(m_settings, &dif_con, settings.get("solver/concentration/diffusion/type"));

    // Pressure
    SolverSelection::SetPressureSolver(m_settings, &pres, settings.get("solver/pressure/type"));

    // Source of velocity
    SolverSelection::SetSourceSolver(m_settings, &sou_vel, settings.get("solver/source/type"));

    // Source of temperature
    SolverSelection::SetSourceSolver(m_settings, &sou_temp, settings.get("solver/temperature/source/type"));

    // Source of concentration
    SolverSelection::SetSourceSolver(m_settings, &sou_con, settings.get("solver/concentration/source/type"));

    // Constants
    m_dir_vel = settings.get("solver/source/dir");

    m_hasDissipation = settings.get_bool("solver/temperature/source/dissipation");
    m_forceFct = settings.get("solver/source/force_fct");
    m_tempFct = settings.get("solver/temperature/source/temp_fct");
    m_conFct = settings.get("solver/concentration/source/con_fct");
    control();
}

NSTempConSolver::~NSTempConSolver() {
    delete adv_vel;
    delete dif_vel;
    delete adv_temp;
    delete dif_temp;
    delete adv_con;
    delete dif_con;
    delete pres;
    delete sou_vel;
    delete sou_temp;
    delete sou_con;
}

//=========================================== do_step ====================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void NSTempConSolver::do_step(real t, bool sync) {
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

    auto nu = m_settings.get_real("physical_parameters/nu");
    auto kappa = m_settings.get_real("physical_parameters/kappa");
    auto gamma = m_settings.get_real("solver/concentration/diffusion/gamma");

    auto dir_vel = m_dir_vel;

#pragma acc data present(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, p, rhs, T, T0, T_tmp, C, C0, C_tmp, f_x, f_y, f_z, S_T, S_C)
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

// 2. Solve diffusion equation
        if (nu != 0.) {
#ifndef BENCHMARKING
            m_logger->info("Diffuse ...");
#endif
            dif_vel->diffuse(u, u0, u_tmp, nu, sync);
            dif_vel->diffuse(v, v0, v_tmp, nu, sync);
            dif_vel->diffuse(w, w0, w_tmp, nu, sync);

            // Couple data to prepare for adding source
            FieldController::couple_vector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);
        }

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
        pres->pressure(p, rhs, t, sync);        //only multigrid cycle, divergence and velocity update (in case of NS) need to be added

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
        if (kappa != 0.) {
#ifndef BENCHMARKING
            m_logger->info("Diffuse Temperature ...");
#endif
            dif_temp->diffuse(T, T0, T_tmp, kappa, sync);

            // Couple temperature to prepare for adding source
            FieldController::couple_scalar(T, T0, T_tmp, sync);
        }

        // Add dissipation
        if (m_hasDissipation) {
#ifndef BENCHMARKING
            m_logger->info("Add dissipation ...");
#endif
            sou_temp->dissipate(m_settings, T, u, v, w, sync);

            // Couple temperature
            FieldController::couple_scalar(T, T0, T_tmp, sync);
        }

        // Add source
        if (m_tempFct != SourceMethods::Zero) {
#ifndef BENCHMARKING
            m_logger->info("Add temperature source ...");
#endif
            sou_temp->add_source(T, S_T, sync);

            // Couple temperature
            FieldController::couple_scalar(T, T0, T_tmp, sync);
        }

// 6. Solve for concentration
        // Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect Concentration ...");
#endif
        adv_con->advect(C, C0, u, v, w, sync);

        // Couple concentration to prepare for diffusion
        FieldController::couple_scalar(C, C0, C_tmp, sync);

        // Solve diffusion equation
        if (gamma != 0.) {
#ifndef BENCHMARKING
            m_logger->info("Diffuse Concentration ...");
#endif
            dif_con->diffuse(C, C0, C_tmp, gamma, sync);

            // Couple concentration to prepare for adding source
            FieldController::couple_scalar(C, C0, C_tmp, sync);
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
void NSTempConSolver::control() {
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

    if (m_settings.get("solver/temperature/advection/field") != Field::get_field_type_name(FieldType::T)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (m_settings.get("solver/temperature/diffusion/field") != Field::get_field_type_name(FieldType::T)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (m_settings.get("solver/concentration/advection/field") != Field::get_field_type_name(FieldType::RHO)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (m_settings.get("solver/concentration/diffusion/field") != Field::get_field_type_name(FieldType::RHO)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (m_settings.get("solver/pressure/field") != Field::get_field_type_name(FieldType::P)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
}
