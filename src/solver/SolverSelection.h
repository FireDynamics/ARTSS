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
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"
#include "../interfaces/ISourceFunction.h"

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
    inline static const std::string Buoyancy = "Buoyancy";
    inline static const std::string Cube = "Cube";
    inline static const std::string Gauss = "Gauss";
    inline static const std::string Uniform = "Uniform";
    inline static const std::string Zero = "Zero";
};

struct TurbulenceMethods {
    inline static const std::string ConstSmagorinsky = "ConstSmagorinsky";
    inline static const std::string DynamicSmagorinsky = "DynamicSmagorinsky";
};

namespace SolverSelection {
    void SetAdvectionSolver(Settings::Settings const &settings,
                            IAdvection **advection_solver,
                            const std::string& advection_type);

    void SetDiffusionSolver(Settings::Settings const &settings,
                            IDiffusion **diffusion_solver,
                            const std::string& diffusion_type);

    void SetPressureSolver(Settings::Settings const &settings,
                           IPressure **pressure_solver,
                           const std::string& pressure_type);

    void SetSourceSolver(Settings::Settings const &settings,
                         ISource **source_solver,
                         const std::string& source_type);

    void SetTurbulenceSolver(Settings::Settings const &settings,
                             ITurbulence **turbulence_solver,
                             const std::string& turbulence_type);

    void set_source_function(const Settings::Settings &settings,
                             ISourceFunction **source_function,
                             const std::string &source_fct);
};

#endif /* ARTSS_SOLVER_SOLVERSELECTION_H_ */

