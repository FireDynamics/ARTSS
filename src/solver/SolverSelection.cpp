/// \file       SolverSelection.cpp
/// \brief      Selects the solver
/// \date       December 18, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "SolverSelection.h"
#include "../advection/SLAdvect.h"
#include "../diffusion/JacobiDiffuse.h"
#include "../diffusion/ColoredGaussSeidelDiffuse.h"
#include "../diffusion/ExplicitDiffuse.h"
#include "../pressure/VCycleMG.h"
#include "../source/ExplicitEulerSource.h"
#include "../turbulence/ConstSmagorinsky.h"
#include "../turbulence/DynamicSmagorinsky.h"


// =================== Set advection solver ==================
// ***************************************************************************************
/// \brief  Sets the advection solver
/// \param  advectionSolver Pointer to AdvectionSolver
/// \param  advectionType Name of AdvcetionSolver
// ***************************************************************************************
void SolverSelection::SetAdvectionSolver(IAdvection **advectionSolver, const std::string& advectionType) {

    if (advectionType == AdvectionMethods::SemiLagrangian) {
        *advectionSolver = new SLAdvect();
    } else {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(SolverSelection).name());
        m_logger->error("Advection method not yet implemented! simulation stopped!");
#endif
        std::exit(1);
        //TODO Error handling
    }
}

// =================== Set diffusion solver ==================
// ***************************************************************************************
/// \brief  Sets the diffusion solver
/// \param  diffusionSolver Pointer to DiffusionSolver
/// \param  diffusionType Name of DiffusionSolver
// ***************************************************************************************
void SolverSelection::SetDiffusionSolver(IDiffusion **diffusionSolver, const std::string& diffusionType) {
    if (diffusionType == DiffusionMethods::Jacobi) {
        *diffusionSolver = new JacobiDiffuse();
    } else if (diffusionType == DiffusionMethods::ColoredGaussSeidel) {
        *diffusionSolver = new ColoredGaussSeidelDiffuse();
    } else if (diffusionType == DiffusionMethods::Explicit) {
        *diffusionSolver = new ExplicitDiffuse();
    } else {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(SolverSelection).name());
        m_logger->error("Diffusion method not yet implemented! Simulation stopped!");
#endif
        std::exit(1);
        //TODO Error handling
    }
}

// =================== Set pressure solver ==================
// ***************************************************************************************
/// \brief  Sets the pressure solver
/// \param  pressureSolver Pointer to PressureSolver
/// \param  pressureType Name of PressureSolver
// ***************************************************************************************
void SolverSelection::SetPressureSolver(IPressure **pressureSolver, const std::string& pressureType, Field *p, Field *rhs) {
    if (pressureType == PressureMethods::VCycleMG) {
        *pressureSolver = new VCycleMG(p, rhs);
    } else {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(SolverSelection).name());
        m_logger->error("Pressure method not yet implemented! Simulation stopped!");
#endif
        std::exit(1);
        //TODO Error handling
    }
}

// =================== Set source solver ==================
// ***************************************************************************************
/// \brief  Sets the source solver
/// \param  sourceSolver Pointer to SourceSolver
/// \param  sourceType Name of SourceSolver
// ***************************************************************************************
void SolverSelection::SetSourceSolver(ISource **sourceSolver, const std::string& sourceType) {
    if (sourceType == SourceMethods::ExplicitEuler) {
        *sourceSolver = new ExplicitEulerSource();
    } else {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(SolverSelection).name());
        m_logger->error("Source method not yet implemented! Simulation stopped!");
#endif
        std::exit(1);
        //TODO Error handling
    }
}

// =================== Set turbulence solver ==================
// ***************************************************************************************
/// \brief  Sets the turbulence solver
/// \param  turbulenceSolver Pointer to TurbulenceSolver
/// \param  turbulenceType Name of TurbulenceSolver
// ***************************************************************************************
void SolverSelection::SetTurbulenceSolver(ITurbulence **turbulenceSolver, const std::string& turbulenceType) {
    if (turbulenceType == TurbulenceMethods::ConstSmagorinsky) {
        *turbulenceSolver = new ConstSmagorinsky();
    } else if (turbulenceType == TurbulenceMethods::DynamicSmagorinsky) {
        *turbulenceSolver = new DynamicSmagorinsky();
    } else {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(SolverSelection).name());
        m_logger->error("Turbulence model is not yet implemented! Simulation stopped!");
#endif
        std::exit(1);
        //TODO Error handling
    }
}
