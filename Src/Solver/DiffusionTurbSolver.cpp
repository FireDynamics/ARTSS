/// \file 		DiffusionTurbSolver.h
/// \brief 		Defines the steps to solve the turbulent diffusion equation
/// \date 		August 18, 2016
/// \author 	Suryanarayana Maddu
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>

#include "DiffusionTurbSolver.h"
#include "../Utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"

DiffusionTurbSolver::DiffusionTurbSolver() {

    auto params = Parameters::getInstance();

    //Diffusion
    std::string diffusionType = params->get("solver/diffusion/type");
    SolverSelection::SetDiffusionSolver(&this->dif, diffusionType);

    m_nu = params->getReal("physical_parameters/nu");

    // Turbulent viscosity
    std::string turbluenceType = params->get("solver/turbulence/type");
    SolverSelection::SetTurbulenceSolver(&this->mu_tub, turbluenceType);
    control();
}


DiffusionTurbSolver::~DiffusionTurbSolver() {
    delete dif;
    delete mu_tub;
}

//====================================== DoStep =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param	dt			time step
/// \param	sync		synchronous kernel launching (true, default: false)
// ***************************************************************************************
void DiffusionTurbSolver::DoStep(real t, bool sync) {

// 1. Solve diffusion equation
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
    auto nu_t = SolverI::nu_t;     //Eddy Viscosity

    auto d_u = u->data;
    auto d_v = v->data;
    auto d_w = w->data;
    auto d_u0 = u0->data;
    auto d_v0 = v0->data;
    auto d_w0 = w0->data;
    auto d_u_tmp = u_tmp->data;
    auto d_v_tmp = v_tmp->data;
    auto d_w_tmp = w_tmp->data;
    auto d_nu_t = nu_t->data;

    size_t bsize = Domain::getInstance()->GetSize(u->GetLevel());

    auto nu = m_nu;

#pragma acc data present(d_u[:bsize], d_u0[:bsize], d_u_tmp[:bsize], d_v[:bsize], d_v0[:bsize], d_v_tmp[:bsize], d_w[:bsize], d_w0[:bsize], d_w_tmp[:bsize], d_nu_t[:bsize]) //EV
    {
#ifndef PROFILING
        std::cout << "Calculating Turbulent viscosity ..." << std::endl;
        //TODO Logger
#endif
        mu_tub->CalcTurbViscosity(nu_t, u, v, w, true);
#ifndef PROFILING
        std::cout << "Diffuse ..." << std::endl;
        //TODO Logger
#endif
        dif->diffuse(u, u0, u_tmp, nu, nu_t, sync);
        dif->diffuse(v, v0, v_tmp, nu, nu_t, sync);
        dif->diffuse(w, w0, w_tmp, nu, nu_t, sync);
    }//end data
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void DiffusionTurbSolver::control() {
    auto params = Parameters::getInstance();
    if (params->get("solver/diffusion/field") != "u,v,w") {
        std::cout << "Fields not specified correctly!" << std::endl;
        std::flush(std::cout);
        std::exit(1);
        //TODO Error handling + Logger
    }
}
