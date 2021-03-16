/// \file       NSSolver.cpp
/// \brief      Defines the (fractional) steps to solve the incompressible Navier-Stokes equations
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "NSSolver.h"
#include "../pressure/VCycleMG.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"
#include "../boundary/BoundaryData.h"

NSSolver::NSSolver(FieldController *field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = field_controller;

    auto params = Parameters::getInstance();

    //Advection of velocity
    std::string advectionType = params->get("solver/advection/type");
    SolverSelection::SetAdvectionSolver(&adv_vel, advectionType);

    //Diffusion of velocity
    std::string diffusionType = params->get("solver/diffusion/type");
    SolverSelection::SetDiffusionSolver(&dif_vel, diffusionType);

    m_nu = params->get_real("physical_parameters/nu");

    //Pressure
    std::string pressureType = params->get("solver/pressure/type");
    SolverSelection::SetPressureSolver(&pres, pressureType, m_field_controller->field_p, m_field_controller->field_rhs);

    //Source
    std::string sourceType = params->get("solver/source/type");
    SolverSelection::SetSourceSolver(&sou, sourceType);

    m_sourceFct = params->get("solver/source/force_fct");
    control();
}

NSSolver::~NSSolver() {
    delete adv_vel;
    delete dif_vel;
    delete pres;
    delete sou;
}

//========================================== do_step ======================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void NSSolver::do_step(real t, bool sync) {
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
    auto f_x = m_field_controller->field_force_x;
    auto f_y = m_field_controller->field_force_y;
    auto f_z = m_field_controller->field_force_z;

    auto d_u = u->data;
    auto d_v = v->data;
    auto d_w = w->data;
    auto d_u0 = u0->data;
    auto d_v0 = v0->data;
    auto d_w0 = w0->data;
    auto d_u_tmp = u_tmp->data;
    auto d_v_tmp = v_tmp->data;
    auto d_w_tmp = w_tmp->data;
    auto d_p = p->data;
    auto d_p0 = p0->data;
    auto d_rhs = rhs->data;
    auto d_fx = f_x->data;
    auto d_fy = f_y->data;
    auto d_fz = f_z->data;

    size_t bsize = Domain::getInstance()->get_size(u->get_level());

    auto nu = m_nu;

#pragma acc data present(d_u[:bsize], d_u0[:bsize], d_u_tmp[:bsize], d_v[:bsize], d_v0[:bsize], d_v_tmp[:bsize], d_w[:bsize], d_w0[:bsize], d_w_tmp[:bsize], d_p[:bsize], d_p0[:bsize], d_rhs[:bsize], d_fx[:bsize], d_fy[:bsize], d_fz[:bsize])
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
        if (m_sourceFct != SourceMethods::Zero) {
#ifndef BENCHMARKING
            m_logger->info("Add source ...");
#endif
            sou->add_source(u, v, w, f_x, f_y, f_z, sync);
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
#ifndef BENCHMARKING
    auto logger = Utility::create_logger(typeid(NSSolver).name());
#endif
    auto params = Parameters::getInstance();

    if (params->get("solver/advection/field") != "u,v,w") {
#ifndef BENCHMARKING
        logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error Handling
    }
    if (params->get("solver/diffusion/field") != "u,v,w") {
#ifndef BENCHMARKING
        logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error Handling
    }
    if (params->get("solver/pressure/field") != BoundaryData::get_field_type_name(FieldType::P)) {
#ifndef BENCHMARKING
        logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error Handling
    }
}
