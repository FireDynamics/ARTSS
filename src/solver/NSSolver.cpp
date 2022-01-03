/// \file       NSSolver.cpp
/// \brief      Defines the (fractional) steps to solve the incompressible Navier-Stokes equations
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "NSSolver.h"

#include <string>
#include <vector>
#include <algorithm>

#include "SolverSelection.h"
#include "../domain/DomainData.h"

NSSolver::NSSolver(const Settings::solver_parameters &solver_settings, FieldController *field_controller) :
        m_solver_settings(solver_settings),
        m_field_controller(field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    //Advection of velocity
    SolverSelection::set_advection_solver(m_solver_settings.advection, &adv_vel);

    //Diffusion of velocity
    SolverSelection::set_diffusion_solver(m_solver_settings.diffusion, &dif_vel);

    SolverSelection::set_pressure_solver(m_solver_settings.pressure, &pres);

    //Source
    SolverSelection::set_source_solver(m_solver_settings.source.type, &sou_vel, m_solver_settings.source.direction);
    m_add_source = m_solver_settings.source.force_fct != SourceMethods::Zero;
    control();
}

NSSolver::~NSSolver() {
    delete adv_vel;
    delete dif_vel;
    delete pres;
    delete sou_vel;
}

//========================================== do_step ======================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void NSSolver::do_step(real t, bool sync) {
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

    real nu = DomainData::getInstance()->get_physical_parameters().nu.value();

#pragma acc data present(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, p, rhs, f_x, f_y, f_z)
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
    }//end data
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void NSSolver::control() {
    auto fields_adv = m_solver_settings.advection.fields;
    std::sort(fields_adv.begin(), fields_adv.end());
    if (fields_adv != std::vector<FieldType>({FieldType::U, FieldType::V, FieldType::W})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error Handling
    }

    auto fields_dif = m_solver_settings.diffusion.fields;
    std::sort(fields_dif.begin(), fields_dif.end());
    if (fields_dif != std::vector<FieldType>({FieldType::U, FieldType::V, FieldType::W})) {
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
        std::exit(1);
        // TODO Error Handling
    }
}
