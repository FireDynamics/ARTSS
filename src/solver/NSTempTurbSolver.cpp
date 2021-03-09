/// \file       NSTempTurbSolver.cpp
/// \brief      Defines the (fractional) steps to solve the turbulent incompressible Navier-Stokes equations with force f(T)
/// \date       Feb 15, 2017
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "NSTempTurbSolver.h"
#include "../pressure/VCycleMG.h"
#include "../utility/Parameters.h"
#include "../utility/Utility.h"
#include "../Domain.h"
#include "SolverSelection.h"
#include "../boundary/BoundaryData.h"

NSTempTurbSolver::NSTempTurbSolver(FieldController *field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = field_controller;

    auto params = Parameters::getInstance();

    // Advection of velocity
    SolverSelection::SetAdvectionSolver(&adv_vel, params->get("solver/advection/type"));

    // Advection of temperature
    SolverSelection::SetAdvectionSolver(&adv_temp, params->get("solver/temperature/advection/type"));

    // Diffusion of velocity
    SolverSelection::SetDiffusionSolver(&dif_vel, params->get("solver/diffusion/type"));

    m_nu = params->get_real("physical_parameters/nu");

    // Turbulent viscosity for velocity diffusion
    SolverSelection::SetTurbulenceSolver(&mu_tub, params->get("solver/turbulence/type"));

    // Diffusion of temperature
    SolverSelection::SetDiffusionSolver(&dif_temp, params->get("solver/temperature/diffusion/type"));

    m_kappa = params->get_real("physical_parameters/kappa");

    // Pressure
    SolverSelection::SetPressureSolver(&pres, params->get("solver/pressure/type"), m_field_controller->field_p, m_field_controller->field_rhs);

    // Source of velocity
    SolverSelection::SetSourceSolver(&sou_vel, params->get("solver/source/type"));

    // Source of temperature
    SolverSelection::SetSourceSolver(&sou_temp, params->get("solver/temperature/source/type"));

    // Constants
    m_dir_vel = params->get("solver/source/dir");
    m_hasTurbulence = (params->get("solver/temperature/turbulence/include") == "Yes");
    m_hasDissipation = (params->get("solver/temperature/source/dissipation") == "Yes");
    m_forceFct = params->get("solver/source/force_fct");
    m_tempFct = params->get("solver/temperature/source/temp_fct");
    control();
}

NSTempTurbSolver::~NSTempTurbSolver() {
    delete adv_vel;
    delete dif_vel;
    delete mu_tub;
    delete adv_temp;
    delete dif_temp;
    delete pres;
    delete sou_vel;
    delete sou_temp;
}

//========================================== do_step ======================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void NSTempTurbSolver::do_step(real t, bool sync) {

    auto params = Parameters::getInstance();

    // local variables and parameters for GPU
    auto u = m_field_controller->field_u;
    auto v = m_field_controller->field_v;
    auto w = m_field_controller->field_w;
    auto u0 = m_field_controller->field_u0;
    auto v0 = m_field_controller->field_v0;
    auto w0 = m_field_controller->field_w0;
    auto u_tmp = m_field_controller->field_u_tmp;
    auto v_tmp = m_field_controller->field_v_tmp;
    auto w_tmp = m_field_controller->field_w_tmp;
    auto p = m_field_controller->field_p;
    auto p0 = m_field_controller->field_p0;
    auto rhs = m_field_controller->field_rhs;
    auto T = m_field_controller->field_T;
    auto T0 = m_field_controller->field_T0;
    auto T_tmp = m_field_controller->field_T_tmp;
    auto f_x = m_field_controller->field_force_x;
    auto f_y = m_field_controller->field_force_y;
    auto f_z = m_field_controller->field_force_z;
    auto S_T = m_field_controller->field_source_T;
    auto nu_t = m_field_controller->field_nu_t;        // nu_t - Eddy Viscosity
    auto kappa_t = m_field_controller->field_kappa_t;  // kappa_t - Eddy thermal diffusivity

    auto d_nu_t = nu_t->data;
    auto d_kappa_t = kappa_t->data;

    size_t bsize = Domain::getInstance()->get_size(u->get_level());

    auto nu = m_nu;
    auto kappa = m_kappa;
    auto dir_vel = m_dir_vel;

#pragma acc data present(u, u0, u_tmp, v, v0, v_tmp, w, \
                         w0, w_tmp, p, p0, rhs, T, T0, T_tmp, \
                         fx, fy, fz, S_T, nu_t, kappa_t)
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
        mu_tub->CalcTurbViscosity(nu_t, u, v, w, true);


#ifndef BENCHMARKING
        m_logger->info("Diffuse ...");
#endif
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
        if (m_hasTurbulence) {
            real Pr_T = params->get_real("solver/temperature/turbulence/Pr_T");
            real rPr_T = 1. / Pr_T;

            // kappa_turb = nu_turb/Pr_turb
#pragma acc parallel loop independent present(d_kappa_t[:bsize], d_nu_t[:bsize]) async
            for (size_t i = 0; i < bsize; ++i) {
                d_kappa_t[i] = d_nu_t[i] * rPr_T;
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
        if (m_hasDissipation) {
#ifndef BENCHMARKING
            m_logger->info("Add dissipation ...");
#endif
            sou_temp->dissipate(T, u, v, w, sync);

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

// 6. Sources updated in Solver::update_sources, TimeIntegration

        if (sync) {
#pragma acc wait
        }
    }
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void NSTempTurbSolver::control() {
    auto params = Parameters::getInstance();
#ifndef BENCHMARKING
    auto logger = Utility::create_logger(typeid(NSTempTurbSolver).name());
#endif

    if (params->get("solver/advection/field") != "u,v,w") {
#ifndef BENCHMARKING
        logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (params->get("solver/diffusion/field") != "u,v,w") {
#ifndef BENCHMARKING
        logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (params->get("solver/temperature/advection/field") != BoundaryData::getFieldTypeName(FieldType::T)) {
#ifndef BENCHMARKING
        logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (params->get("solver/temperature/diffusion/field") != BoundaryData::getFieldTypeName(FieldType::T)) {
#ifndef BENCHMARKING
        logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
    if (params->get("solver/pressure/field") != BoundaryData::getFieldTypeName(FieldType::P)) {
#ifndef BENCHMARKING
        logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
}
