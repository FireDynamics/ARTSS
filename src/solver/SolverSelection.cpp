/// \file 		SolverSelection.cpp
/// \brief 		Selects the solver
/// \date 		December 18, 2018
/// \author 	My Linh Wuerzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include "SolverSelection.h"
#include "../advection/SLAdvect.h"
#include "../advection/ExplicitAdvect.h"
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
/// \param	advectionSolver Pointer to AdvectionSolver
/// \param  advectionType Name of AdvcetionSolver
// ***************************************************************************************
void SolverSelection::SetAdvectionSolver(AdvectionI **advectionSolver, const std::string& advectionType) {
    if (advectionType == AdvectionMethods::SemiLagrangian) {
        *advectionSolver = new SLAdvect();
    } else if (advectionType == AdvectionMethods::Explicit) {
        *advectionSolver = new ExplicitAdvect();
    } else {
        std::cout << "Advection method not yet implemented! Simulation stopped!" << std::endl;
        std::flush(std::cout);
        std::exit(1);
        //TODO Error handling + Logger
    }
}

// =================== Set diffusion solver ==================
// ***************************************************************************************
/// \brief  Sets the diffusion solver
/// \param	diffusionSolver Pointer to DiffusionSolver
/// \param  diffusionType Name of DiffusionSolver
// ***************************************************************************************
void SolverSelection::SetDiffusionSolver(DiffusionI **diffusionSolver, const std::string& diffusionType) {
    if (diffusionType == DiffusionMethods::Jacobi) {
        *diffusionSolver = new JacobiDiffuse();
    } else if (diffusionType == DiffusionMethods::ColoredGaussSeidel) {
        *diffusionSolver = new ColoredGaussSeidelDiffuse();
    } else if (diffusionType == DiffusionMethods::Explicit) {
        *diffusionSolver = new ExplicitDiffuse();
    } else {
        std::cout << "Diffusion method not yet implemented! Simulation stopped!" << std::endl;
        std::flush(std::cout);
        std::exit(1);
        //TODO Error handling + Logger
    }
}

// =================== Set pressure solver ==================
// ***************************************************************************************
/// \brief  Sets the pressure solver
/// \param	pressureSolver Pointer to PressureSolver
/// \param  pressureType Name of PressureSolver
// ***************************************************************************************
void SolverSelection::SetPressureSolver(PressureI **pressureSolver, const std::string& pressureType, Field *p, Field *rhs) {
    if (pressureType == PressureMethods::VCycleMG) {
        *pressureSolver = new VCycleMG(p, rhs);
    } else {
        std::cout << "Pressure method not yet implemented! Simulation stopped!" << std::endl;
        std::flush(std::cout);
        std::exit(1);
        //TODO Error handling + Logger
    }
}

// =================== Set source solver ==================
// ***************************************************************************************
/// \brief  Sets the source solver
/// \param	sourceSolver Pointer to SourceSolver
/// \param  sourceType Name of SourceSolver
// ***************************************************************************************
void SolverSelection::SetSourceSolver(SourceI **sourceSolver, const std::string& sourceType) {
    if (sourceType == SourceMethods::ExplicitEuler) {
        *sourceSolver = new ExplicitEulerSource();
    } else {
        std::cout << "Source method not yet implemented! Simulation stopped!" << std::endl;
        std::flush(std::cout);
        std::exit(1);
        //TODO Error handling + Logger
    }
}

// =================== Set turbulence solver ==================
// ***************************************************************************************
/// \brief  Sets the turbulence solver
/// \param	turbulenceSolver Pointer to TurbulenceSolver
/// \param  turbulenceType Name of TurbulenceSolver
// ***************************************************************************************
void SolverSelection::SetTurbulenceSolver(TurbulenceI **turbulenceSolver, const std::string& turbulenceType) {
    if (turbulenceType == TurbulenceMethods::ConstSmagorinsky) {
        *turbulenceSolver = new ConstSmagorinsky();
    } else if (turbulenceType == TurbulenceMethods::DynamicSmagorinsky) {
        *turbulenceSolver = new DynamicSmagorinsky();
    } else {
        std::cout << "Turbulence model is not yet implemented! Simulation stopped!" << std::endl;
        std::flush(std::cout);
        std::exit(1);
        //TODO Error handling + Logger
    }
}
