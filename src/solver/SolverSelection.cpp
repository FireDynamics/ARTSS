/// \file       SolverSelection.cpp
/// \brief      Selects the solver
/// \date       December 18, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "SolverSelection.h"
#include <string>
#include "../advection/SLAdvect.h"
#include "../diffusion/JacobiDiffuse.h"
#include "../diffusion/ColoredGaussSeidelDiffuse.h"
#include "../diffusion/ExplicitDiffuse.h"
#include "../pressure/VCycleMG.h"
#include "../source/ExplicitEulerSource.h"
#include "../turbulence/ConstSmagorinsky.h"
#include "../turbulence/DynamicSmagorinsky.h"


namespace SolverSelection {
static const std::string class_name = "SolverSelection";

// =================== Set advection solver ==================
// ***************************************************************************************
/// \brief  Sets the advection solver
/// \param  advectionSolver Pointer to AdvectionSolver
/// \param  advectionType Name of AdvcetionSolver
// ***************************************************************************************
void SetAdvectionSolver(Settings::Settings const &settings, IAdvection **advectionSolver, const std::string& advectionType) {
    if (advectionType == AdvectionMethods::SemiLagrangian) {
        *advectionSolver = new SLAdvect(settings);
    } else {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(class_name);
        logger->error("Advection method not yet implemented! simulation stopped!");
#endif
        std::exit(1);
        // TODO Error handling
    }
}

// =================== Set diffusion solver ==================
// ***************************************************************************************
/// \brief  Sets the diffusion solver
/// \param  diffusionSolver Pointer to DiffusionSolver
/// \param  diffusionType Name of DiffusionSolver
// ***************************************************************************************
void SetDiffusionSolver(Settings::Settings const &settings, IDiffusion **diffusionSolver, const std::string& diffusionType) {
    if (diffusionType == DiffusionMethods::Jacobi) {
        *diffusionSolver = new JacobiDiffuse(settings);
    } else if (diffusionType == DiffusionMethods::ColoredGaussSeidel) {
        *diffusionSolver = new ColoredGaussSeidelDiffuse(settings);
    } else if (diffusionType == DiffusionMethods::Explicit) {
        *diffusionSolver = new ExplicitDiffuse(settings);
    } else {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(class_name);
        logger->error("Diffusion method not yet implemented! Simulation stopped!");
#endif
        std::exit(1);
        // TODO Error handling
    }
}

// =================== Set pressure solver ==================
// ***************************************************************************************
/// \brief  Sets the pressure solver
/// \param  pressureSolver Pointer to PressureSolver
/// \param  pressureType Name of PressureSolver
// ***************************************************************************************
void SetPressureSolver(Settings::Settings const &settings,
                       IPressure **pressureSolver,
                       const std::string& pressureType) {
    if (pressureType == PressureMethods::VCycleMG) {
        *pressureSolver = new VCycleMG(settings);
    } else {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(class_name);
        logger->error("Pressure method not yet implemented! Simulation stopped!");
#endif
        std::exit(1);
        // TODO Error handling
    }
}

// =================== Set source solver ==================
// ***************************************************************************************
/// \brief  Sets the source solver
/// \param  sourceSolver Pointer to SourceSolver
/// \param  sourceType Name of SourceSolver
// ***************************************************************************************
void SetSourceSolver(Settings::Settings const &settings, ISource **sourceSolver, const std::string& sourceType) {
    if (sourceType == SourceMethods::ExplicitEuler) {
        *sourceSolver = new ExplicitEulerSource(settings);
    } else {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(class_name);
        logger->error("Source method {} not yet implemented! Simulation stopped!", sourceType);
#endif
        std::exit(1);
        // TODO Error handling
    }
}

// =================== Set turbulence solver ==================
// ***************************************************************************************
/// \brief  Sets the turbulence solver
/// \param  turbulence_solver Pointer to TurbulenceSolver
/// \param  turbulenceType Name of TurbulenceSolver
// ***************************************************************************************
void SetTurbulenceSolver(Settings::Settings const &settings, ITurbulence **turbulence_solver, const std::string& turbulenceType) {
    if (turbulenceType == TurbulenceMethods::ConstSmagorinsky) {
        *turbulence_solver = new ConstSmagorinsky(settings);
    } else if (turbulenceType == TurbulenceMethods::DynamicSmagorinsky) {
        *turbulence_solver = new DynamicSmagorinsky(settings);
    } else {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(class_name);
        logger->error("Turbulence model is not yet implemented! Simulation stopped!");
#endif
        std::exit(1);
        // TODO Error handling
    }
}
}  // namespace SolverSelection
