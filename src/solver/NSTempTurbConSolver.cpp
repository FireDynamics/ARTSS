/// \file 		NSTempTurbConSolver.cpp
/// \brief 		Defines the (fractional) steps to solve the incompressible Navier-Stokes equations with force f(T), turbulence and concentration C
/// \date 		Oct 02, 2017
/// \author 	Küsters
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "NSTempTurbConSolver.h"
#include "../pressure/VCycleMG.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"
#include "../boundary/BoundaryData.h"

NSTempTurbConSolver::NSTempTurbConSolver() {

    auto params = Parameters::getInstance();

    // Advection of velocity
    SolverSelection::SetAdvectionSolver(&adv_vel, params->get("solver/advection/type"));

    // Advection of temperature
    SolverSelection::SetAdvectionSolver(&adv_temp, params->get("solver/temperature/advection/type"));

    // Advection of concentration
    SolverSelection::SetAdvectionSolver(&adv_con, params->get("solver/concentration/advection/type"));

    // Diffusion of velocity
    SolverSelection::SetDiffusionSolver(&dif_vel, params->get("solver/diffusion/type"));

    m_nu = params->get_real("physical_parameters/nu");

    // Turbulent viscosity for velocity diffusion
    SolverSelection::SetTurbulenceSolver(&mu_tub, params->get("solver/turbulence/type"));

    // Diffusion of temperature
    SolverSelection::SetDiffusionSolver(&dif_temp, params->get("solver/temperature/diffusion/type"));

    m_kappa = params->get_real("physical_parameters/kappa");

    // Diffusion for concentration
    SolverSelection::SetDiffusionSolver(&dif_con, params->get("solver/concentration/diffusion/type"));

    m_gamma = params->get_real("solver/concentration/diffusion/gamma");

    // Pressure
    SolverSelection::SetPressureSolver(&pres, params->get("solver/pressure/type"), p, rhs);

    // Source of velocity
    SolverSelection::SetSourceSolver(&sou_vel, params->get("solver/source/type"));

    // Source of temperature
    SolverSelection::SetSourceSolver(&sou_temp, params->get("solver/temperature/source/type"));

    // Source of concentration
    SolverSelection::SetSourceSolver(&sou_con, params->get("solver/concentration/source/type"));

    // Constants
    m_dir_vel = params->get("solver/source/dir");
    m_has_turbulence_temperature = (params->get("solver/temperature/turbulence/include") == "Yes");
    m_has_turbulence_concentration = (params->get("solver/concentration/turbulence/include") == "Yes");
    m_has_dissipation = (params->get("solver/temperature/source/dissipation") == "Yes");
    m_forceFct = params->get("solver/source/force_fct");
    m_tempFct = params->get("solver/temperature/source/temp_fct");
    m_conFct = params->get("solver/concentration/source/con_fct");
    control();
}

NSTempTurbConSolver::~NSTempTurbConSolver() {
    delete adv_vel;
    delete dif_vel;
    delete adv_temp;
    delete dif_temp;
    delete adv_con;
    delete dif_con;
    delete mu_tub;
    delete pres;
    delete sou_vel;
    delete sou_temp;
    delete sou_con;
}

//====================================== do_step =================================
// ***************************************************************************************
/// \brief  brings all calculation steps together into one function
/// \param	dt			time step
/// \param  sync		synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void NSTempTurbConSolver::do_step(real t, bool sync) {
    auto params = Parameters::getInstance();

    // local variables and parameters for GPU
    auto u = ISolver::u;
    auto v = ISolver::v;
    auto w = ISolver::w;
    auto u0 = ISolver::u0;
    auto v0 = ISolver::v0;
    auto w0 = ISolver::w0;
    auto u_tmp = ISolver::u_tmp;
    auto v_tmp = ISolver::v_tmp;
    auto w_tmp = ISolver::w_tmp;
    auto p = ISolver::p;
    auto p0 = ISolver::p0;
    auto rhs = ISolver::rhs;
    auto T = ISolver::T;
    auto T0 = ISolver::T0;
    auto T_tmp = ISolver::T_tmp;
    auto C = ISolver::concentration;
    auto C0 = ISolver::concentration0;
    auto C_tmp = ISolver::concentration_tmp;
    auto f_x = ISolver::f_x;
    auto f_y = ISolver::f_y;
    auto f_z = ISolver::f_z;
    auto S_T = ISolver::S_T;
    auto S_C = ISolver::S_concentration;
    auto nu_t = ISolver::nu_t;            //nu_t - Eddy Viscosity
    auto kappa_t = ISolver::kappa_t;        //kappa_t - Eddy thermal diffusivity
    auto gamma_t = ISolver::gamma_t;        //gamma_t - Eddy mass diffsusivity

    auto d_u = u->data;
    auto d_v = v->data;
    auto d_w = w->data;
    auto d_u0 = u0->data;
    auto d_v0 = v0->data;
    auto d_w0 = w0->data;
    auto d_u_tmp = u_tmp->data;
    auto d_v_tmp = v_tmp->data;
    auto d_w_tmp = w_tmp->data;
    auto d_p = p->data;
    auto d_p0 = p0->data;
    auto d_rhs = rhs->data;
    auto d_T = T->data;
    auto d_T0 = T0->data;
    auto d_T_tmp = T_tmp->data;
    auto d_C = C->data;
    auto d_C0 = C0->data;
    auto d_C_tmp = C_tmp->data;
    auto d_fx = f_x->data;
    auto d_fy = f_y->data;
    auto d_fz = f_z->data;
    auto d_S_T = S_T->data;
    auto d_S_C = S_C->data;
    auto d_nu_t = nu_t->data;
    auto d_kappa_t = kappa_t->data;
    auto d_gamma_t = gamma_t->data;

    size_t bsize = Domain::getInstance()->get_size(u->GetLevel());

    auto nu = m_nu;
    auto kappa = m_kappa;
    auto gamma = m_gamma;
    auto dir_vel = m_dir_vel;

#ifndef BENCHMARKING
    auto m_logger = Utility::create_logger(typeid(NSTempTurbConSolver).name());
#endif

#pragma acc data present(    d_u[:bsize], d_u0[:bsize], d_u_tmp[:bsize], d_v[:bsize], d_v0[:bsize], d_v_tmp[:bsize], d_w[:bsize], \
                            d_w0[:bsize], d_w_tmp[:bsize], d_p[:bsize], d_p0[:bsize], d_rhs[:bsize], d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], \
                            d_C[:bsize], d_C0[:bsize], d_C_tmp[:bsize], d_fx[:bsize], d_fy[:bsize], d_fz[:bsize], d_S_T[:bsize], d_S_C[:bsize], \
                            d_nu_t[:bsize], d_kappa_t[:bsize], d_gamma_t[:bsize])
    {
// 1. Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect ...");
#endif
        adv_vel->advect(u, u0, u0, v0, w0, sync);
        adv_vel->advect(v, v0, u0, v0, w0, sync);
        adv_vel->advect(w, w0, u0, v0, w0, sync);

        // Couple velocity to prepare for diffusion
        ISolver::couple_vector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);

// 2. Solve turbulent diffusion equation
#ifndef BENCHMARKING
        m_logger->info("Calculating Turbulent viscosity ...");
#endif
        mu_tub->CalcTurbViscosity(nu_t, u, v, w, true);

#ifndef BENCHMARKING
        m_logger->info("Diffuse ...");
#endif
        dif_vel->diffuse(u, u0, u_tmp, nu, nu_t, sync);
        dif_vel->diffuse(v, v0, v_tmp, nu, nu_t, sync);
        dif_vel->diffuse(w, w0, w_tmp, nu, nu_t, sync);

        // Couple data to prepare for adding source
        ISolver::couple_vector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);

// 3. Add force
        if (m_forceFct != SourceMethods::Zero) {
#ifndef BENCHMARKING
            m_logger->info("Add momentum source ...");
#endif
            sou_vel->add_source(u, v, w, f_x, f_y, f_z, sync);

            // Couple data to prepare for adding source
            ISolver::couple_vector(u, u0, u_tmp, v, v0, v_tmp, w, w0, w_tmp, sync);
        }

// 4. Solve pressure equation and project
        // Calculate divergence of u
        pres->divergence(rhs, u_tmp, v_tmp, w_tmp, sync);

        // Solve pressure equation
#ifndef BENCHMARKING
        m_logger->info("Pressure ...");
#endif
        pres->pressure(p, rhs, t, sync);

        // Correct
        pres->projection(u, v, w, u_tmp, v_tmp, w_tmp, p, sync);

// 5. Solve Temperature and link back to force

        // Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect Temperature ...");
#endif
        adv_temp->advect(T, T0, u, v, w, sync);

        // Couple temperature to prepare for diffusion
        ISolver::couple_scalar(T, T0, T_tmp, sync);

        // Solve diffusion equation
        // turbulence
        if (m_has_turbulence_temperature) {
            real Pr_T = params->get_real("solver/temperature/turbulence/Pr_T");
            real rPr_T = 1. / Pr_T;

#pragma acc parallel loop independent present(d_kappa_t[:bsize], d_nu_t[:bsize]) async
            for (size_t i = 0; i < bsize; ++i) {
                d_kappa_t[i] = d_nu_t[i] * rPr_T; // kappa_turb = nu_turb/Pr_turb
            }

#ifndef BENCHMARKING
            m_logger->info("Diffuse turbulent Temperature ...");
#endif
            dif_temp->diffuse(T, T0, T_tmp, kappa, kappa_t, sync);

            // Couple temperature to prepare for adding source
            ISolver::couple_scalar(T, T0, T_tmp, sync);
        } else {
            // no turbulence
            if (kappa != 0.) {

#ifndef BENCHMARKING
                m_logger->info("Diffuse Temperature ...");
#endif
                dif_temp->diffuse(T, T0, T_tmp, kappa, sync);

                // Couple temperature to prepare for adding source
                ISolver::couple_scalar(T, T0, T_tmp, sync);
            }
        }

        // Add dissipation
        if (m_has_dissipation) {

#ifndef BENCHMARKING
            m_logger->info("Add dissipation ...");
#endif
            sou_temp->dissipate(T, u, v, w, sync);

            // Couple temperature
            ISolver::couple_scalar(T, T0, T_tmp, sync);
        }

            // Add source
        else if (m_tempFct != SourceMethods::Zero) {

#ifndef BENCHMARKING
            m_logger->info("Add temperature source ...");
#endif
            sou_temp->add_source(T, S_T, sync);

            // Couple temperature
            ISolver::couple_scalar(T, T0, T_tmp, sync);
        }

// 6. Solve Concentration

        // Solve advection equation
#ifndef BENCHMARKING
        m_logger->info("Advect Concentration ...");
#endif
        adv_con->advect(C, C0, u, v, w, sync);

        // Couple concentration to prepare for diffusion
        ISolver::couple_scalar(C, C0, C_tmp, sync);

        // Solve diffusion equation
        // turbulence
        if (m_has_turbulence_concentration) {
            real Sc_T = params->get_real("solver/concentration/turbulence/Sc_T");
            real rSc_T = 1. / Sc_T;

#pragma acc parallel loop independent present(d_gamma_t[:bsize], d_nu_t[:bsize]) async
            for (size_t i = 0; i < bsize; ++i) d_gamma_t[i] = d_nu_t[i] * rSc_T; // gamma_turb = nu_turb/Sc_turb

#ifndef BENCHMARKING
            m_logger->info("Diffuse turbulent Concentration ...");
#endif
            dif_con->diffuse(C, C0, C_tmp, gamma, gamma_t, sync);

            // Couple concentration to prepare for adding source
            ISolver::couple_scalar(C, C0, C_tmp, sync);
        } else {
            // no turbulence
            if (gamma != 0.) {
#ifndef BENCHMARKING
                m_logger->info("Diffuse Concentration ...");
#endif
                dif_con->diffuse(C, C0, C_tmp, gamma, sync);

                // Couple concentration to prepare for adding source
                ISolver::couple_scalar(C, C0, C_tmp, sync);
            }
        }

        // Add source
        if (m_conFct != SourceMethods::Zero) {

#ifndef BENCHMARKING
            m_logger->info("Add concentration source ...");
#endif
            sou_con->add_source(C, S_C, sync);

            // Couple concentration
            ISolver::couple_scalar(C, C0, C_tmp, sync);
        }

// 7. Sources updated in Solver::update_sources, TimeIntegration

        if (sync) {
#pragma acc wait
        }
    }//end data
}

//======================================= Check data ==================================
// ***************************************************************************************
/// \brief  Checks if field specified correctly
// ***************************************************************************************
void NSTempTurbConSolver::control() {
    auto params = Parameters::getInstance();
#ifndef BENCHMARKING
    auto m_logger = Utility::create_logger(typeid(NSTempTurbConSolver).name());
#endif
    if (params->get("solver/advection/field") != "u,v,w") {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
    if (params->get("solver/diffusion/field") != "u,v,w") {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
    if (params->get("solver/temperature/advection/field") != BoundaryData::getFieldTypeName(FieldType::T)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
    if (params->get("solver/concentration/advection/field") != BoundaryData::getFieldTypeName(FieldType::RHO)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
    if (params->get("solver/temperature/diffusion/field") != BoundaryData::getFieldTypeName(FieldType::T)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
    if (params->get("solver/concentration/diffusion/field") != BoundaryData::getFieldTypeName(FieldType::RHO)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
    if (params->get("solver/pressure/field") != BoundaryData::getFieldTypeName(FieldType::P)) {
#ifndef BENCHMARKING
        m_logger->error("Fields not specified correctly!");
#endif
        std::exit(1);
        //TODO Error handling
    }
}
