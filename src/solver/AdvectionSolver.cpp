/// \file       AdvectionSolver.cpp
/// \brief      Defines the steps to solve the advection equation
/// \date       August 22, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <string>
#include <vector>
#include <algorithm>

#include "AdvectionSolver.h"
#include "../interfaces/IAdvection.h"
#include "../Domain.h"
#include "SolverSelection.h"

AdvectionSolver::AdvectionSolver(
        Settings::Settings const &settings,
        FieldController *field_controller,
        real u_lin, real v_lin, real w_lin) :
    m_settings(settings),
    m_field_controller(field_controller),
    m_u_lin(FieldType::U, u_lin),
    m_v_lin(FieldType::V, v_lin),
    m_w_lin(FieldType::W, w_lin) {
#ifndef BENCHMARKING
     m_logger = Utility::create_logger(settings, typeid(this).name());
#endif
    std::string advectionType = m_settings.get("solver/advection/type");
    SolverSelection::SetAdvectionSolver(m_settings, &adv, m_settings.get("solver/advection/type"));

    m_u_lin.copyin();
    m_v_lin.copyin();
    m_w_lin.copyin();

    control();
}

AdvectionSolver::AdvectionSolver(Settings::Settings const &settings, FieldController *field_controller) :
    AdvectionSolver(
            settings,
            field_controller,
            settings.get_real("initial_conditions/u_lin"),
            settings.get_real("initial_conditions/v_lin"),
            settings.get_real("initial_conditions/w_lin")) {
}

AdvectionSolver::~AdvectionSolver() {
    delete adv;
}

//====================================== do_step =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronous kernel launching (true, default: false)
// ***************************************************************************************
void AdvectionSolver::do_step(real, bool sync) {
    Field &u = m_field_controller->get_field_u();
    Field &v = m_field_controller->get_field_v();
    Field &w = m_field_controller->get_field_w();
    Field &u0 = m_field_controller->get_field_u0();
    Field &v0 = m_field_controller->get_field_v0();
    Field &w0 = m_field_controller->get_field_w0();

#pragma acc data present(u_lin, v_lin, w_lin, u, u0, v, v0, w, w0)
    {
// 1. Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect ...");
#endif
        adv->advect(u, u0, m_u_lin, m_v_lin, m_w_lin, sync);
        adv->advect(v, v0, m_u_lin, m_v_lin, m_w_lin, sync);
        adv->advect(w, w0, m_u_lin, m_v_lin, m_w_lin, sync);
    }
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void AdvectionSolver::control() {
    auto fields = Utility::split(m_settings.get("solver/advection/field"), ',');
    std::sort(fields.begin(), fields.end());

    if (fields != std::vector<std::string>({"u", "v", "w"})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
}
