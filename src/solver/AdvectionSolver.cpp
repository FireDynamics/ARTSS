/// \file       AdvectionSolver.cpp
/// \brief      Defines the steps to solve the advection equation
/// \date       August 22, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "AdvectionSolver.h"
#include "../interfaces/IAdvection.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"

AdvectionSolver::AdvectionSolver() {

    auto params = Parameters::getInstance();
    std::string advectionType = params->get("solver/advection/type");
    SolverSelection::SetAdvectionSolver(&adv, params->get("solver/advection/type"));

    real d_u_linm = params->get_real("initial_conditions/u_lin");
    real d_v_linm = params->get_real("initial_conditions/v_lin");
    real d_w_linm = params->get_real("initial_conditions/w_lin");

    u_linm = new Field(FieldType::U, d_u_linm);
    v_linm = new Field(FieldType::V, d_v_linm);
    w_linm = new Field(FieldType::W, d_w_linm);

    auto u_lin = u_linm;
    auto v_lin = v_linm;
    auto w_lin = w_linm;

    auto d_u_lin = u_lin->data;
    auto d_v_lin = v_lin->data;
    auto d_w_lin = w_lin->data;

    size_t bsize = Domain::getInstance()->get_size(u_linm->GetLevel());

#pragma acc enter data copyin(d_u_lin[:bsize], d_v_lin[:bsize], d_w_lin[:bsize])
    control();
}

AdvectionSolver::~AdvectionSolver() {
    delete adv;

    auto d_u_lin = u_linm->data;
    auto d_v_lin = v_linm->data;
    auto d_w_lin = w_linm->data;

    size_t bsize = Domain::getInstance()->get_size(u_linm->GetLevel());

#pragma acc exit data delete(d_u_lin[:bsize], d_v_lin[:bsize], d_w_lin[:bsize])

    delete u_linm;
    delete v_linm;
    delete w_linm;
}

//====================================== do_step =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronous kernel launching (true, default: false)
// ***************************************************************************************
void AdvectionSolver::do_step(real t, bool sync) {
  // local variables and parameters for GPU
    auto u = ISolver::u;
    auto v = ISolver::v;
    auto w = ISolver::w;
    auto u0 = ISolver::u0;
    auto v0 = ISolver::v0;
    auto w0 = ISolver::w0;

    auto d_u = u->data;
    auto d_v = v->data;
    auto d_w = w->data;
    auto d_u0 = u0->data;
    auto d_v0 = v0->data;
    auto d_w0 = w0->data;

    auto u_lin = u_linm;
    auto v_lin = v_linm;
    auto w_lin = w_linm;

    auto d_u_lin = u_lin->data;
    auto d_v_lin = v_lin->data;
    auto d_w_lin = w_lin->data;

    size_t bsize = Domain::getInstance()->get_size(u->GetLevel());

#pragma acc data present(d_u_lin[:bsize], d_v_lin[:bsize], d_w_lin[:bsize], d_u[:bsize], d_u0[:bsize], d_v[:bsize], d_v0[:bsize], d_w[:bsize], d_w0[:bsize])
    {
// 1. Solve advection equation
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(AdvectionSolver).name());
        m_logger->info("Advect ...");
#endif
        adv->advect(u, u0, u_lin, v_lin, w_lin, sync);
        adv->advect(v, v0, u_lin, v_lin, w_lin, sync);
        adv->advect(w, w0, u_lin, v_lin, w_lin, sync);
    }//end data
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void AdvectionSolver::control() {
    auto params = Parameters::getInstance();
    if (params->get("solver/advection/field") != "u,v,w") {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(AdvectionSolver).name());
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
}
