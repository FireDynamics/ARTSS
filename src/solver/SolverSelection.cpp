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
#include "../randomField/UniformRandom.h"
#include "../source/Zero.h"
#include "../source/Cube.h"
#include "../source/BuoyancyMMS.h"
#include "../source/GaussFunction.h"


namespace SolverSelection {
    static const std::string class_name = "SolverSelection";

// =================== Set advection solver ==================
// ***************************************************************************************
/// \brief  Sets the advection solver
/// \param  advection_solver Pointer to AdvectionSolver
/// \param  advection_type Name of AdvectionSolver
// ***************************************************************************************
    void set_advection_solver(const Settings::solver::advection_solver &settings,
                              IAdvection **advection_solver) {
        if (settings.type == AdvectionMethods::SemiLagrangian) {
            *advection_solver = new SLAdvect();
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
/// \param  diffusion_solver Pointer to DiffusionSolver
/// \param  diffusion_type Name of DiffusionSolver
// ***************************************************************************************
    void set_diffusion_solver(const Settings::solver::diffusion_solver &settings,
                              IDiffusion **diffusion_solver) {
        if (settings.type == DiffusionMethods::Jacobi) {
            const auto &jacobi = std::get<Settings::solver::diffusion_solvers::jacobi>(settings.solver.value());
            *diffusion_solver = new JacobiDiffuse(jacobi);
        } else if (settings.type == DiffusionMethods::ColoredGaussSeidel) {
            const auto &cgs = std::get<Settings::solver::diffusion_solvers::colored_gauss_seidel>(settings.solver.value());
            *diffusion_solver = new ColoredGaussSeidelDiffuse(cgs);
        } else if (settings.type == DiffusionMethods::Explicit) {
            *diffusion_solver = new ExplicitDiffuse();
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
/// \param  pressure_solver Pointer to PressureSolver
/// \param  pressure_type Name of PressureSolver
// ***************************************************************************************
    void set_pressure_solver(const Settings::solver::pressure_solver &settings,
                             IPressure **pressure_solver) {
        if (settings.type == PressureMethods::VCycleMG) {
            *pressure_solver = new VCycleMG(settings.solver);
        } else {
#ifndef BENCHMARKING
            auto logger = Utility::create_logger(class_name);
            logger->error("Pressure method not yet implemented! Simulation stopped!");
#endif
            std::exit(1);
            // TODO Error handling
        }
    }

// =================== Set turbulence solver ==================
// ***************************************************************************************
/// \brief  Sets the turbulence solver
/// \param  turbulence_solver Pointer to TurbulenceSolver
/// \param  turbulence_type Name of TurbulenceSolver
// ***************************************************************************************
    void set_turbulence_solver(const Settings::solver::turbulence_solver &settings,
                               ITurbulence **turbulence_solver) {
        std::string turbulence_type = settings.type;
        if (turbulence_type == TurbulenceMethods::ConstSmagorinsky) {
            *turbulence_solver = new ConstSmagorinsky(settings.solver.value());
        } else if (turbulence_type == TurbulenceMethods::DynamicSmagorinsky) {
            *turbulence_solver = new DynamicSmagorinsky();
        } else {
#ifndef BENCHMARKING
            auto logger = Utility::create_logger(class_name);
            logger->error("Turbulence model is not yet implemented! Simulation stopped!");
#endif
            std::exit(1);
            // TODO Error handling
        }
    }

    void add_noise(const Settings::random_parameters &random_parameters, ISourceFunction **source_function) {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(class_name);
        if (!random_parameters.absolute && (random_parameters.range > 1 || random_parameters.range  < 0)) {
            logger->error("range for relative noise has to be between [0:1], given value: {}", random_parameters.range);
            std::exit(1);
        }
#else
        if (!random_parameters.absolute && (random_parameters.range > 1 || random_parameters.range  < 0)) {
            std::cout << "range for relative noise has to be between [0:1], given value: {}" << random_parameters.range << std::endl;
            std::exit(1);
        }
#endif
        IRandomField *noise_maker;
        if (random_parameters.custom_seed) {
            noise_maker = new UniformRandom(random_parameters.range, random_parameters.step_size, random_parameters.seed);
        } else {
            noise_maker = new UniformRandom(random_parameters.range, random_parameters.step_size);
        }
        (*source_function)->set_noise(noise_maker, random_parameters.absolute);
    }

    void set_temperature_source_function(const Settings::solver::temperature_source &settings,
                                         ISourceFunction **source_function) {
        std::string source_fct = settings.temp_fct;
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(class_name);
        logger->debug("create source function {}", source_fct);
#endif
        if (source_fct == SourceMethods::Gauss) {
            const auto &gauss = std::get<Settings::solver::sources::gauss>(settings.temp_function);
            *source_function = new GaussFunction(gauss);
        } else if (source_fct == SourceMethods::Buoyancy) {
            *source_function = new BuoyancyMMS();
        } else if (source_fct == SourceMethods::Cube) {
            const auto &cube = std::get<Settings::solver::sources::cube>(settings.temp_function);
            *source_function = new Cube(cube);
        } else if (source_fct == SourceMethods::Zero) {
            *source_function = new Zero();
        } else {
#ifndef BENCHMARKING
            logger->warn("Source method {} not yet implemented!", source_fct);
#endif
        }
        if (settings.random) {
            add_noise(settings.random_parameters, source_function);
        }
    }

    void set_concentration_source_function(const Settings::solver::concentration_source &settings,
                                           ISourceFunction **source_function) {
        std::string source_fct = settings.con_fct;
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(class_name);
        logger->debug("create source function {}", source_fct);
#endif
        if (source_fct == SourceMethods::Gauss) {
            const auto &gauss = std::get<Settings::solver::sources::gauss>(settings.con_function);
            *source_function = new GaussFunction(gauss);
        } else if (source_fct == SourceMethods::Buoyancy) {
            *source_function = new BuoyancyMMS();
        } else if (source_fct == SourceMethods::Cube) {
            const auto &cube = std::get<Settings::solver::sources::cube>(settings.con_function);
            *source_function = new Cube(cube);
        } else if (source_fct == SourceMethods::Zero) {
            *source_function = new Zero();
        } else {
#ifndef BENCHMARKING
            logger->warn("Source method {} not yet implemented!", source_fct);
#endif
        }
        if (settings.random) {
            add_noise(settings.random_parameters, source_function);
        }
    }


    // =================== Set source solver ==================
    // ***************************************************************************************
    /// \brief  Sets the source solver
    /// \param  source_solver Pointer to SourceSolver
    /// \param  source_type Name of SourceSolver
    // ***************************************************************************************
    void set_source_solver(const std::string &source_type,
                           ISource **source_solver,
                           const std::vector<CoordinateAxis> &dir) {
        if (source_type == SourceMethods::ExplicitEuler) {
            *source_solver = new ExplicitEulerSource(dir);
        } else {
#ifndef BENCHMARKING
            auto logger = Utility::create_logger(class_name);
            logger->error("Source method {} not yet implemented! Simulation stopped!", source_type);
#endif
            std::exit(1);
            // TODO Error handling
        }
    }
}  // namespace SolverSelection
