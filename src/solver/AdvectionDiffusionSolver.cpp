/// \file     AdvectionDiffusionSolver.h
/// \brief    Defines the steps to solve the advection and diffusion equation
/// \date     May 20, 2016
/// \author   Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "AdvectionDiffusionSolver.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"
#include "../utility/Utility.h"


AdvectionDiffusionSolver::AdvectionDiffusionSolver(FieldController *field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = field_controller;

    auto params = Parameters::getInstance();
    std::string advectionType = params->get("solver/advection/type");
    SolverSelection::SetAdvectionSolver(&this->adv, advectionType);

    std::string diffusionType = params->get("solver/diffusion/type");
    SolverSelection::SetDiffusionSolver(&this->dif, diffusionType);

    m_nu = params->get_real("physical_parameters/nu");

    control();
}

AdvectionDiffusionSolver::~AdvectionDiffusionSolver() {
    delete adv;
    delete dif;
}


//================================= DoStep ============================
// *******************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt          time step
/// \param  sync        synchronous kernel launching (true, default: false)
// *******************************************************************
void AdvectionDiffusionSolver::do_step(real t, bool sync) {
// local variables and parameters
    auto u = m_field_controller->field_u;
    auto v = m_field_controller->field_v;
    auto w = m_field_controller->field_w;
    auto u0 = m_field_controller->field_u0;
    auto v0 = m_field_controller->field_v0;
    auto w0 = m_field_controller->field_w0;
    auto u_tmp = m_field_controller->field_u_tmp;
    auto v_tmp = m_field_controller->field_v_tmp;
    auto w_tmp = m_field_controller->field_w_tmp;

    auto d_u = u->data;
    auto d_v = v->data;
    auto d_w = w->data;
    auto d_u0 = u0->data;
    auto d_v0 = v0->data;
    auto d_w0 = w0->data;
    auto d_u_tmp = u_tmp->data;
    auto d_v_tmp = v_tmp->data;
    auto d_w_tmp = w_tmp->data;

    size_t bsize = Domain::getInstance()->get_size(u->get_level());

    auto nu = m_nu;

#pragma acc data present(d_u[:bsize], d_u0[:bsize], d_u_tmp[:bsize], d_v[:bsize], d_v0[:bsize], d_v_tmp[:bsize], d_w[:bsize], d_w0[:bsize], d_w_tmp[:bsize])
    {
// 1. Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect ...");
#endif
        adv->advect(*u, *u0, *u0, *v0, *w0, sync);
        adv->advect(*v, *v0, *u0, *v0, *w0, sync);
        adv->advect(*w, *w0, *u0, *v0, *w0, sync);

// 2. Couple data to prepare for diffusion
        FieldController::couple_vector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);

// 3. Solve diffusion equation
        if (nu != 0.) {
#ifndef BENCHMARKING
            m_logger->info("Diffuse ...");
#endif
            dif->diffuse(*u, *u0, *u_tmp, nu, sync);
            dif->diffuse(*v, *v0, *v_tmp, nu, sync);
            dif->diffuse(*w, *w0, *w_tmp, nu, sync);
        }

// if sync kernel launching, wait for kernel to be finished before starting next kernel launch
        if (sync) {
#pragma acc wait
        }
    }//end data
}

//================================== Check data =============================
// ************************************************************************
/// \brief  Checks if field specified correctly
// ************************************************************************
void AdvectionDiffusionSolver::control() {
    auto params = Parameters::getInstance();
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(typeid(AdvectionDiffusionSolver).name());
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
}
