/// \file       NSTempSolver.cpp
/// \brief      Defines the steps to solve advection, diffusion, pressure and add sources (dependent on T), solves temperature equation
/// \date       Feb 15, 2017
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "NSTempSolver.h"

#include <string>
#include <vector>
#include <algorithm>

#include "SolverSelection.h"
#include "../domain/DomainData.h"

NSTempSolver::NSTempSolver(const Settings::solver_parameters &solver_settings, FieldController *field_controller) :
        m_solver_settings(solver_settings),
        m_field_controller(field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    // Advection of velocity
    SolverSelection::set_advection_solver(m_solver_settings.advection, &adv_vel);

    // Advection of temperature
    SolverSelection::set_advection_solver(m_solver_settings.temperature.advection, &adv_temp);

    // Diffusion of velocity
    SolverSelection::set_diffusion_solver(m_solver_settings.diffusion, &dif_vel);

    // Diffusion of temperature
    SolverSelection::set_diffusion_solver(m_solver_settings.diffusion, &dif_temp);

    // Pressure
    SolverSelection::set_pressure_solver(m_solver_settings.pressure, &pres);

    // Source of velocity
    SolverSelection::set_source_solver(m_solver_settings.source.type, &sou_vel, m_solver_settings.source.direction);

    // Source of temperature
    SolverSelection::set_source_solver(m_solver_settings.temperature.source.type, &sou_temp, m_solver_settings.temperature.source.dir);

    // temperature function
    SolverSelection::set_temperature_source_function(m_solver_settings.temperature.source, &m_source_function_temperature);
    m_add_source = m_solver_settings.source.force_fct != SourceMethods::Zero;
    m_add_temp_source = m_solver_settings.temperature.source.temp_fct != SourceMethods::Zero;
    control();
}

// ==================================== Destructor ====================================
// ***************************************************************************************
NSTempSolver::~NSTempSolver() {
    delete adv_vel;
    delete dif_vel;
    delete adv_temp;
    delete dif_temp;
    delete pres;
    delete sou_vel;
    delete sou_temp;
}

//====================================== do_step =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void NSTempSolver::do_step(real t, bool sync) {
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
    Field &f_x = m_field_controller->get_field_force_x();
    Field &f_y = m_field_controller->get_field_force_y();
    Field &f_z = m_field_controller->get_field_force_z();
    Field &S_T = m_field_controller->get_field_source_T();

    auto domain_data = DomainData::getInstance();
    real nu = domain_data->get_physical_parameters().nu.value();
    real kappa = domain_data->get_physical_parameters().kappa.value();

#pragma acc data present(u, u0, u_tmp, v, v0, v_tmp, w, \
                            w0, w_tmp, p, rhs, T, T0, T_tmp, \
                            f_x, f_y, f_z, S_T)
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
        if (m_add_source) {
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
        if (m_solver_settings.temperature.source.dissipation) {

#ifndef BENCHMARKING
            m_logger->info("Add dissipation ...");
#endif
            sou_temp->dissipate(T, u, v, w, sync);

            // Couple temperature
            FieldController::couple_scalar(T, T0, T_tmp, sync);
        }

        // Add source
        if (m_add_temp_source) {

#ifndef BENCHMARKING
            m_logger->info("Add temperature source ...");
#endif
            sou_temp->add_source(T, S_T, sync);

            // Couple temperature
            FieldController::couple_scalar(T, T0, T_tmp, sync);
        }

// 6. Sources updated in Solver::update_sources, TimeIntegration

        if (sync) {
#pragma acc wait
        }
    }//end data
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void NSTempSolver::control() {
    auto adv_fields = m_solver_settings.advection.fields;
    std::sort(adv_fields.begin(), adv_fields.end());
    if (adv_fields != std::vector<FieldType>({FieldType::U, FieldType::V, FieldType::W})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error Handling
    }

    auto diff_fields = m_solver_settings.diffusion.fields;
    std::sort(diff_fields.begin(), diff_fields.end());
    if (diff_fields != std::vector<FieldType>({FieldType::U, FieldType::V, FieldType::W})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error Handling
    }

    if (m_solver_settings.temperature.advection.fields != std::vector<FieldType>{FieldType::T}) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error Handling
    }
    if (m_solver_settings.temperature.diffusion.fields != std::vector<FieldType>{FieldType::T}) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error Handling
    }
    if (m_solver_settings.pressure.field != FieldType::P) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        // TODO Error Handling
        std::exit(1);
    }
}

void NSTempSolver::update_source(real t_cur) {
    m_source_function_temperature->update_source(m_field_controller->get_field_source_T(), t_cur);
}

void NSTempSolver::replace_heat_source(const Settings::solver::temperature_source &temperature_source) {
    SolverSelection::set_temperature_source_function(temperature_source, &m_source_function_temperature);
    m_add_temp_source = m_solver_settings.temperature.source.temp_fct != SourceMethods::Zero;
}

void NSTempSolver::update_obstacle_change() {
    m_source_function_temperature->update_obstacle_change();
}
