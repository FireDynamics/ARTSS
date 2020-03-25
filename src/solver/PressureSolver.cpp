/// \file 		PressureSolver.cpp
/// \brief 		Defines the steps to solve the pressure Poisson equation
/// \date 		Sep 13, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include <spdlog/spdlog.h>

#include "PressureSolver.h"
#include "../pressure/VCycleMG.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"
#include "../boundary/BoundaryData.h"

PressureSolver::PressureSolver() {
    auto params = Parameters::getInstance();
    SolverSelection::SetPressureSolver(&this->pres, params->get("solver/pressure/type"), p, rhs);
    control();
}

PressureSolver::~PressureSolver() {
    delete pres;
}

//========================================= DoStep =======================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param	dt			time step
/// \param  sync		synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void PressureSolver::DoStep(real t, bool sync) {

#ifndef PROFILING
    spdlog::info("Pressure ...");
#endif

// 1. Solve pressure Poisson equation

    // local variables and parameters for GPU
    auto p = SolverI::p;
    auto rhs = SolverI::rhs;
    auto d_p = p->data;
    auto d_rhs = rhs->data;

    size_t bsize = Domain::getInstance()->GetSize(p->GetLevel());

#pragma acc data present(d_p[:bsize], d_rhs[:bsize])
    {
        pres->pressure(p, rhs, t, sync);
    }//end data
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void PressureSolver::control() {
    auto params = Parameters::getInstance();
    if (params->get("solver/pressure/field") != BoundaryData::getFieldTypeName(FieldType::P)) {
        spdlog::error("Fields not specified correctly!");
        std::exit(1);
        //TODO Error handling
    }
}
