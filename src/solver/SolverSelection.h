/// \file       SolverSelection.cpp
/// \brief      Selects the solver
/// \date       December 18, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_SOLVERSELECTION_H
#define ARTSS_SOLVER_SOLVERSELECTION_H

#include "../interfaces/IDiffusion.h"
#include "../interfaces/ISource.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/IAdvection.h"
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
    inline static const std::string GaussST = "GaussST";
    inline static const std::string GaussSC = "GaussSC";
    inline static const std::string Uniform = "Uniform";
    inline static const std::string Zero = "Zero";
};

struct TurbulenceMethods {
    inline static const std::string ConstSmagorinsky = "ConstSmagorinsky";
    inline static const std::string DynamicSmagorinsky = "DynamicSmagorinsky";
};

class SolverSelection {
 public:
    static void SetAdvectionSolver(IAdvection **advectionSolver, const std::string& advectionType);

    static void SetDiffusionSolver(IDiffusion **diffusionSolver, const std::string& diffusionType);

    static void SetPressureSolver(IPressure **pressureSolver, const std::string& pressureType, Field *p, Field *rhs);

    static void SetSourceSolver(ISource **sourceSolver, const std::string& sourceType);

    static void SetTurbulenceSolver(ITurbulence **tubulenceSolver, const std::string& turbulenceType);
};


#endif /* ARTSS_SOLVER_SOLVERSELECTION_H */
