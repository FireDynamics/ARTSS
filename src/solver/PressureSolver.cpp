/// \file       PressureSolver.cpp
/// \brief      Defines the steps to solve the pressure Poisson equation
/// \date       Sep 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "PressureSolver.h"


PressureSolver::PressureSolver(Settings::Settings const &settings, FieldController *field_controller) :
        m_settings(settings) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(PressureSolver).name());
#endif
    m_field_controller = field_controller;
    auto p_type = m_settings.get("solver/pressure/type");
    SolverSelection::SetPressureSolver(m_settings, &this->pres, p_type);
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
    m_logger->info("Pressure ...");
#endif

    // 1. Solve pressure Poisson equation
    // local variables and parameters for GPU
    Field &p = m_field_controller->get_field_p();
    Field &rhs = m_field_controller->get_field_rhs();
#pragma acc data present(p, rhs)
    {
        pres->pressure(p, rhs, t, sync);
    }
}

//================================== Check data =============================
// *****************************************************************************
/// \brief  Checks if field specified correctly
// *****************************************************************************
void PressureSolver::control() {
    auto p_field = m_settings.get("solver/pressure/field");
    if (p_field != Field::get_field_type_name(FieldType::P)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO(issue 6) Error handling
    }
}
