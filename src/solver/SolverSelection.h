/// \file       SolverSelection.cpp
/// \brief      Selects the solver
/// \date       December 18, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_SOLVERSELECTION_H_
#define ARTSS_SOLVER_SOLVERSELECTION_H_

#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/ISource.h"
#include "../interfaces/ITurbulence.h"
#include "../utility/Utility.h"

struct AdvectionMethods {
    inline static const std::string SemiLagrangian = "SemiLagrangian";
};

struct DiffusionMethods {
    inline static const std::string ColoredGaussSeidel = "ColoredGaussSeidel";
    inline static const std::string Jacobi = "Jacobi";
    inline static const std::string Explicit = "Explicit";
};

struct PressureMethods {
    inline static const std::string VCycleMG = "VCycleMG";
};

struct SourceMethods {
    inline static const std::string ExplicitEuler = "ExplicitEuler";
    inline static const std::string BuoyancyST_MMS = "BuoyancyST_MMS";
    inline static const std::string Buoyancy = "Buoyancy";
    inline static const std::string Cube = "Cube";
    inline static const std::string GaussST = "GaussST";
    inline static const std::string GaussSC = "GaussSC";
    inline static const std::string Uniform = "Uniform";
    inline static const std::string Zero = "Zero";
};

struct TurbulenceMethods {
    inline static const std::string ConstSmagorinsky = "ConstSmagorinsky";
    inline static const std::string DynamicSmagorinsky = "DynamicSmagorinsky";
};

namespace SolverSelection {
    void SetAdvectionSolver(Settings const &settings, IAdvection **advectionSolver, const std::string& advectionType);

    void SetDiffusionSolver(Settings const &settings, IDiffusion **diffusionSolver, const std::string& diffusionType);

    void SetPressureSolver(Settings const &settings, IPressure **pressureSolver, const std::string& pressureType,
                           const Field &p, const Field &rhs);

    void SetSourceSolver(Settings const &settings, ISource **sourceSolver, const std::string& sourceType);

    void SetTurbulenceSolver(Settings const &settings, ITurbulence **tubulenceSolver, const std::string& turbulenceType);
};

#endif /* ARTSS_SOLVER_SOLVERSELECTION_H_ */

