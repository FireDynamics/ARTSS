/// \file       PressureSolver.cpp
/// \brief      Defines the steps to solve the pressure Poisson equation
/// \date       Sep 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>

#include "PressureSolver.h"
#include "../pressure/VCycleMG.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"
#include "../boundary/BoundaryData.h"

PressureSolver::PressureSolver(FieldController *field_controller) {
    m_field_controller = field_controller;
    auto params = Parameters::getInstance();
    SolverSelection::SetPressureSolver(&this->pres, params->get("solver/pressure/type"), m_field_controller->field_p, m_field_controller->field_rhs);
    control();
}

PressureSolver::~PressureSolver() {
    delete pres;
}

//========================================= do_step =======================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void PressureSolver::do_step(real t, bool sync) {

#ifndef BENCHMARKING
    std::cout << "Pressure ..." << std::endl;
    //TODO Logger
#endif

// 1. Solve pressure Poisson equation

    // local variables and parameters for GPU
    auto p = m_field_controller->field_p;
    auto rhs = m_field_controller->field_rhs;
    auto d_p = p->data;
    auto d_rhs = rhs->data;

    size_t bsize = Domain::getInstance()->get_size(p->get_level());

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
        std::cout << "Fields not specified correctly!" << std::endl;
        std::flush(std::cout);
        std::exit(1);
        //TODO Error handling + Logger
    }
}
