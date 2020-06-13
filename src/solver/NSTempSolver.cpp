/// \file 		NSTempSolver.cpp
/// \brief 		Defines the (fractional) steps to solve the incompressible Navier-Stokes equations with force f(T)
/// \date 		Feb 15, 2017
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>

#include "NSTempSolver.h"
#include "../pressure/VCycleMG.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"
#include "../boundary/BoundaryData.h"

NSTempSolver::NSTempSolver() {

    auto params = Parameters::getInstance();

    // Advection of velocity
    SolverSelection::SetAdvectionSolver(&adv_vel, params->get("solver/advection/type"));

    // Advection of temperature
    SolverSelection::SetAdvectionSolver(&adv_temp, params->get("solver/temperature/advection/type"));

    // Diffusion of velocity
    SolverSelection::SetDiffusionSolver(&dif_vel, params->get("solver/diffusion/type"));

    m_nu = params->getReal("physical_parameters/nu");

    // Diffusion of temperature
    SolverSelection::SetDiffusionSolver(&dif_temp, params->get("solver/temperature/diffusion/type"));

    m_kappa = params->getReal("physical_parameters/kappa");

    // Pressure
    SolverSelection::SetPressureSolver(&pres, params->get("solver/pressure/type"), p, rhs);

    // Source of velocity
    SolverSelection::SetSourceSolver(&sou_vel, params->get("solver/source/type"));

    // Source of temperature
    SolverSelection::SetSourceSolver(&sou_temp, params->get("solver/temperature/source/type"));

    // Constants
    m_dir_vel = params->get("solver/source/dir");

    m_hasDissipation = (params->get("solver/temperature/source/dissipation") == "Yes");
    m_forceFct = params->get("solver/source/force_fct");
    m_tempFct = params->get("solver/temperature/source/temp_fct");
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

//====================================== DoStep =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param	dt			time step
/// \param  sync		synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void NSTempSolver::DoStep(real t, bool sync) {

    // local variables and parameters for GPU
    auto u = SolverI::u;
    auto v = SolverI::v;
    auto w = SolverI::w;
    auto u0 = SolverI::u0;
    auto v0 = SolverI::v0;
    auto w0 = SolverI::w0;
    auto u_tmp = SolverI::u_tmp;
    auto v_tmp = SolverI::v_tmp;
    auto w_tmp = SolverI::w_tmp;
    auto p = SolverI::p;
    auto p0 = SolverI::p0;
    auto rhs = SolverI::rhs;
    auto T = SolverI::T;
    auto T0 = SolverI::T0;
    auto T_tmp = SolverI::T_tmp;
    auto f_x = SolverI::f_x;
    auto f_y = SolverI::f_y;
    auto f_z = SolverI::f_z;
    auto S_T = SolverI::S_T;

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
    auto d_T = T->data;
    auto d_T0 = T0->data;
    auto d_T_tmp = T_tmp->data;
    auto d_fx = f_x->data;
    auto d_fy = f_y->data;
    auto d_fz = f_z->data;
    auto d_S_T = S_T->data;

    size_t bsize = Domain::getInstance()->GetSize(u->GetLevel());

    auto nu = m_nu;
    auto kappa = m_kappa;
    auto dir_vel = m_dir_vel;

#pragma acc data present(    d_u[:bsize], d_u0[:bsize], d_u_tmp[:bsize], d_v[:bsize], d_v0[:bsize], d_v_tmp[:bsize], d_w[:bsize], \
                            d_w0[:bsize], d_w_tmp[:bsize], d_p[:bsize], d_p0[:bsize], d_rhs[:bsize], d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], \
                            d_fx[:bsize], d_fy[:bsize], d_fz[:bsize], d_S_T[:bsize])
    {
// 1. Solve advection equation
#ifndef BENCHMARKING
        std::cout << "Advect ..." << std::endl;
        //TODO Logger
#endif
        adv_vel->advect(u, u0, u0, v0, w0, sync);
        adv_vel->advect(v, v0, u0, v0, w0, sync);
        adv_vel->advect(w, w0, u0, v0, w0, sync);


        // Couple velocity to prepare for diffusion
        SolverI::CoupleVector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);

// 2. Solve diffusion equation
        if (nu != 0.) {
#ifndef BENCHMARKING
            std::cout << "Diffuse ..." << std::endl;
            //TODO Logger
#endif
            dif_vel->diffuse(u, u0, u_tmp, nu, sync);
            dif_vel->diffuse(v, v0, v_tmp, nu, sync);
            dif_vel->diffuse(w, w0, w_tmp, nu, sync);

            // Couple data to prepare for adding source
            SolverI::CoupleVector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);
        }

// 3. Add force
        if (m_forceFct != SourceMethods::Zero) {
#ifndef BENCHMARKING
            std::cout << "Add momentum source ..." << std::endl;
            //TODO Logger
#endif
            sou_vel->addSource(u, v, w, f_x, f_y, f_z, sync);

            // Couple data to prepare for adding source
            SolverI::CoupleVector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);
        }

// 4. Solve pressure equation and project
        // Calculate divergence of u
        pres->Divergence(rhs, u_tmp, v_tmp, w_tmp, sync);

        // Solve pressure equation
#ifndef BENCHMARKING
        std::cout << "Pressure ..." << std::endl;
        //TODO Logger
#endif
        pres->pressure(p, rhs, t, sync);        //only multigrid cycle, divergence and velocity update (in case of NS) need to be added

        // Correct
        pres->Project(u, v, w, u_tmp, v_tmp, w_tmp, p, sync);

// 5. Solve Temperature and link back to force
        // Solve advection equation
#ifndef BENCHMARKING
        std::cout << "Advect Temperature ..." << std::endl;
        //TODO Logger
#endif
        adv_temp->advect(T, T0, u, v, w, sync);

        // Couple temperature to prepare for diffusion
        SolverI::CoupleScalar(T, T0, T_tmp, sync);

        // Solve diffusion equation
        if (kappa != 0.) {

#ifndef BENCHMARKING
            std::cout << "Diffuse Temperature ..." << std::endl;
            //TODO Logger
#endif
            dif_temp->diffuse(T, T0, T_tmp, kappa, sync);

            // Couple temperature to prepare for adding source
            SolverI::CoupleScalar(T, T0, T_tmp, sync);
        }

        // Add dissipation
        if (m_hasDissipation) {

#ifndef BENCHMARKING
            std::cout << "Add dissipation ..." << std::endl;
            //TODO Logger
#endif
            sou_temp->Dissipate(T, u, v, w, sync);

            // Couple temperature
            SolverI::CoupleScalar(T, T0, T_tmp, sync);
        }

        // Add source
        if (m_tempFct != SourceMethods::Zero) {

#ifndef BENCHMARKING
            std::cout << "Add temperature source ..." << std::endl;
            //TODO Logger
#endif
            sou_temp->addSource(T, S_T, sync);

            // Couple temperature
            SolverI::CoupleScalar(T, T0, T_tmp, sync);
        }

// 6. Sources updated in Solver::UpdateSources, TimeIntegration

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
    auto params = Parameters::getInstance();
    if (params->get("solver/advection/field") != "u,v,w") {
        std::cout << "Fields not specified correctly!" << std::endl;
        //TODO Logger + Error Handling
        std::flush(std::cout);
        std::exit(1);
    }
    if (params->get("solver/diffusion/field") != "u,v,w") {
        std::cout << "Fields not specified correctly!" << std::endl;
        //TODO Logger + Error Handling
        std::flush(std::cout);
        std::exit(1);
    }
    if (params->get("solver/temperature/advection/field") != BoundaryData::getFieldTypeName(FieldType::T)) {
        std::cout << "Fields not specified correctly!" << std::endl;
        //TODO Logger + Error Handling
        std::flush(std::cout);
        std::exit(1);
    }
    if (params->get("solver/temperature/diffusion/field") != BoundaryData::getFieldTypeName(FieldType::T)) {
        std::cout << "Fields not specified correctly!" << std::endl;
        //TODO Logger + Error Handling
        std::flush(std::cout);
        std::exit(1);
    }
    if (params->get("solver/pressure/field") != BoundaryData::getFieldTypeName(FieldType::P)) {
        std::cout << "Fields not specified correctly!" << std::endl;
        //TODO Logger + Error Handling
        std::flush(std::cout);
        std::exit(1);
    }
}
