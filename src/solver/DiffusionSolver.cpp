/// \file       DiffusionSolver.cpp
/// \brief      Defines the steps to solve the diffusion equation
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <string>
#include <vector>
#include <algorithm>

#include "DiffusionSolver.h"
#include "../interfaces/IDiffusion.h"
#include "../Domain.h"
#include "SolverSelection.h"

DiffusionSolver::DiffusionSolver(Settings::Settings const &settings, FieldController *field_controller) :
        m_settings(settings) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(m_settings, typeid(this).name());
#endif
    m_field_controller = field_controller;

    std::string diffusionType = m_settings.get("solver/diffusion/type");
    SolverSelection::SetDiffusionSolver(m_settings, &this->dif, diffusionType);

    control();
}

DiffusionSolver::~DiffusionSolver() {
    delete dif;
}

//====================================== do_step =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param  dt      time step
/// \param  sync    synchronous kernel launching (true, default: false)
// ***************************************************************************************
void DiffusionSolver::do_step(real, bool sync) {
#ifndef BENCHMARKING
    m_logger->info("Diffuse ...");
#endif
    // 1. Solve diffusion equation
    Field &u = m_field_controller->get_field_u();
    Field &v = m_field_controller->get_field_v();
    Field &w = m_field_controller->get_field_w();
    Field &u0 = m_field_controller->get_field_u0();
    Field &v0 = m_field_controller->get_field_v0();
    Field &w0 = m_field_controller->get_field_w0();
    Field &u_tmp = m_field_controller->get_field_u_tmp();
    Field &v_tmp = m_field_controller->get_field_v_tmp();
    Field &w_tmp = m_field_controller->get_field_w_tmp();

#pragma acc data present(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp)
    {
        real nu = m_settings.get_real("physical_parameters/nu");
        dif->diffuse(u, u0, u_tmp, nu, sync);
        dif->diffuse(v, v0, v_tmp, nu, sync);
        dif->diffuse(w, w0, w_tmp, nu, sync);
    }
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void DiffusionSolver::control() {
    auto fields = Utility::split(m_settings.get("solver/diffusion/field"), ',');
    std::sort(fields.begin(), fields.end());
    if (fields != std::vector<std::string>({"u", "v", "w"})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
}
