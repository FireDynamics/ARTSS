/// \file       PressureSolver.cpp
/// \brief      Defines the steps to solve the pressure Poisson equation
/// \date       Sep 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#include "PressureSolver.h"


PressureSolver::PressureSolver() {
    auto params = Parameters::getInstance();
    auto p_type = params->get("solver/pressure/type");
    SolverSelection::SetPressureSolver(&this->pres, p_type, p, rhs);
    control();
}

PressureSolver::~PressureSolver() {
    delete pres;
}

//==================================== DoStep ==================================
// *****************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void PressureSolver::do_step(real t, bool sync) {
#ifndef BENCHMARKING
    auto m_logger = Utility::create_logger(typeid(PressureSolver).name());
    m_logger->info("Pressure ...");
#endif

// 1. Solve pressure Poisson equation

    // local variables and parameters for GPU
    auto p = ISolver::p;
    auto rhs = ISolver::rhs;
    auto d_p = p->data;
    auto d_rhs = rhs->data;

    size_t bsize = Domain::getInstance()->get_size(p->GetLevel());

#pragma acc data present(d_p[:bsize], d_rhs[:bsize])
    {
        pres->pressure(p, rhs, t, sync);
    }  // end data
}

//================================== Check data =============================
// *****************************************************************************
/// \brief  Checks if field specified correctly
// *****************************************************************************
void PressureSolver::control() {
    auto params = Parameters::getInstance();
    auto p_field = params->get("solver/pressure/field");
    if (p_field != BoundaryData::getFieldTypeName(FieldType::P)) {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(PressureSolver).name());
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
}
