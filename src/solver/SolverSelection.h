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
    void set_advection_solver(const Settings::solver::advection_solver &settings,
                              IAdvection **advection_solver);

    void set_diffusion_solver(const Settings::solver::diffusion_solver &settings,
                            IDiffusion **diffusion_solver);

    void set_pressure_solver(const Settings::solver::pressure_solver &settings,
                             IPressure **pressure_solver);

    void set_source_solver(const std::string &source_type,
                           ISource **source_solver,
                           const std::vector<CoordinateAxis> &dir);

    void set_turbulence_solver(const Settings::solver::turbulence_solver &settings,
                               ITurbulence **turbulence_solver);

    void set_temperature_source_function(const Settings::solver::temperature_source &settings,
                                         ISourceFunction **source_function);
    void set_concentration_source_function(const Settings::solver::concentration_source &settings,
                                          ISourceFunction **source_function);
};

#endif /* ARTSS_SOLVER_SOLVERSELECTION_H_ */

