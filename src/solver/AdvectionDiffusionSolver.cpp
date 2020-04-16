/// \file 		AdvectionDiffusionSolver.h
/// \brief 		Defines the steps to solve the advection and diffusion equation
/// \date 		May 20, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <spdlog/spdlog.h>

#include "AdvectionDiffusionSolver.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"
#include "../utility/Utility.h"

AdvectionDiffusionSolver::AdvectionDiffusionSolver() {
#ifndef PROFILING
    m_logger = Utility::createLogger(typeid(this).name());
#endif
    auto params = Parameters::getInstance();
    std::string advectionType = params->get("solver/advection/type");
    SolverSelection::SetAdvectionSolver(&this->adv, advectionType);

    std::string diffusionType = params->get("solver/diffusion/type");
    SolverSelection::SetDiffusionSolver(&this->dif, diffusionType);

    m_nu = params->getReal("physical_parameters/nu");

    control();
}

AdvectionDiffusionSolver::~AdvectionDiffusionSolver() {
    delete adv;
    delete dif;
}


//====================================== DoStep =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param	dt			time step
/// \param	sync		synchronous kernel launching (true, default: false)
// ***************************************************************************************

void AdvectionDiffusionSolver::DoStep(real t, bool sync) {

// local variables and parameters
    auto u = SolverI::u;
    auto v = SolverI::v;
    auto w = SolverI::w;
    auto u0 = SolverI::u0;
    auto v0 = SolverI::v0;
    auto w0 = SolverI::w0;
    auto u_tmp = SolverI::u_tmp;
    auto v_tmp = SolverI::v_tmp;
    auto w_tmp = SolverI::w_tmp;

    auto d_u = u->data;
    auto d_v = v->data;
    auto d_w = w->data;
    auto d_u0 = u0->data;
    auto d_v0 = v0->data;
    auto d_w0 = w0->data;
    auto d_u_tmp = u_tmp->data;
    auto d_v_tmp = v_tmp->data;
    auto d_w_tmp = w_tmp->data;

    size_t bsize = Domain::getInstance()->GetSize(u->GetLevel());

    auto nu = m_nu;

#pragma acc data present(d_u[:bsize], d_u0[:bsize], d_u_tmp[:bsize], d_v[:bsize], d_v0[:bsize], d_v_tmp[:bsize], d_w[:bsize], d_w0[:bsize], d_w_tmp[:bsize])
    {
// 1. Solve advection equation
#ifndef PROFILING
        spdlog::info("Advect ...");
#endif
        adv->advect(u, u0, u0, v0, w0, sync);
        adv->advect(v, v0, u0, v0, w0, sync);
        adv->advect(w, w0, u0, v0, w0, sync);

// 2. Couple data to prepare for diffusion
        SolverI::CoupleVector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);

// 3. Solve diffusion equation
        if (nu != 0.) {
#ifndef PROFILING
            spdlog::info("Diffuse ...");
#endif
            dif->diffuse(u, u0, u_tmp, nu, sync);
            dif->diffuse(v, v0, v_tmp, nu, sync);
            dif->diffuse(w, w0, w_tmp, nu, sync);
        }

// if sync kernel launching, wait for kernel to be finished before starting next kernel launch
        if (sync) {
#pragma acc wait
        }
    }//end data
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void AdvectionDiffusionSolver::control() {
    auto params = Parameters::getInstance();
    if (params->get("solver/advection/field") != "u,v,w") {
        spdlog::error("Fields not specified correctly!");
        std::exit(1);
        //TODO Error handling
    }
    if (params->get("solver/diffusion/field") != "u,v,w") {
        spdlog::error("Fields not specified correctly!");
        std::exit(1);
        //TODO Error handling
    }
}
