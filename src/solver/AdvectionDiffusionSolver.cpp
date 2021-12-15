/// \file     AdvectionDiffusionSolver.h
/// \brief    Defines the steps to solve the advection and diffusion equation
/// \date     May 20, 2016
/// \author   Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "AdvectionDiffusionSolver.h"

#include <string>
#include <vector>
#include <algorithm>

#include "SolverSelection.h"
#include "../domain/DomainData.h"


AdvectionDiffusionSolver::AdvectionDiffusionSolver(Settings::Settings const &settings, FieldController *field_controller) :
        m_settings(settings) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = field_controller;

    std::string advectionType = m_settings.get("solver/advection/type");
    SolverSelection::SetAdvectionSolver(m_settings, &this->adv, advectionType);

    std::string diffusionType = m_settings.get("solver/diffusion/type");
    SolverSelection::SetDiffusionSolver(m_settings, &this->dif, diffusionType);

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
void AdvectionDiffusionSolver::do_step(real, bool sync) {
// local variables and parameters
    Field &u = m_field_controller->get_field_u();
    Field &v = m_field_controller->get_field_v();
    Field &w = m_field_controller->get_field_w();
    Field &u0 = m_field_controller->get_field_u0();
    Field &v0 = m_field_controller->get_field_v0();
    Field &w0 = m_field_controller->get_field_w0();
    Field &u_tmp = m_field_controller->get_field_u_tmp();
    Field &v_tmp = m_field_controller->get_field_v_tmp();
    Field &w_tmp = m_field_controller->get_field_w_tmp();

    auto nu = m_settings.get_real("physical_parameters/nu");

#pragma acc data present(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp)
    {
// 1. Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect ...");
#endif
        adv->advect(u, u0, u0, v0, w0, sync);
        adv->advect(v, v0, u0, v0, w0, sync);
        adv->advect(w, w0, u0, v0, w0, sync);

// 2. Couple data to prepare for diffusion
        FieldController::couple_vector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);

// 3. Solve diffusion equation
        if (nu != 0.) {
#ifndef BENCHMARKING
            m_logger->info("Diffuse ...");
#endif
            dif->diffuse(u, u0, u_tmp, nu, sync);
            dif->diffuse(v, v0, v_tmp, nu, sync);
            dif->diffuse(w, w0, w_tmp, nu, sync);
        }

// if sync kernel launching, wait for kernel to be finished before starting next kernel launch
        if (sync) {
#pragma acc wait
        }
    }
}

//================================== Check data =============================
// ************************************************************************
/// \brief  Checks if field specified correctly
// ************************************************************************
void AdvectionDiffusionSolver::control() {
    auto fields = Utility::split(m_settings.get("solver/advection/field"), ',');
    std::sort(fields.begin(), fields.end());
    if (fields != std::vector<std::string>({"u", "v", "w"})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }

    auto diff_fields = Utility::split(m_settings.get("solver/diffusion/field"), ',');
    std::sort(diff_fields.begin(), diff_fields.end());
    if (diff_fields != std::vector<std::string>({"u", "v", "w"})) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        // TODO Error handling
    }
}
