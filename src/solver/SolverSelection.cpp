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
/// \param  advection_type Name of AdvcetionSolver
// ***************************************************************************************
void SetAdvectionSolver(Settings::Settings const &settings,
                        IAdvection **advection_solver,
                        const std::string& advection_type) {
    if (advection_type == AdvectionMethods::SemiLagrangian) {
        *advection_solver = new SLAdvect(settings);
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
void SetDiffusionSolver(Settings::Settings const &settings,
                        IDiffusion **diffusion_solver,
                        const std::string& diffusion_type) {
    if (diffusion_type == DiffusionMethods::Jacobi) {
        *diffusion_solver = new JacobiDiffuse(settings);
    } else if (diffusion_type == DiffusionMethods::ColoredGaussSeidel) {
        *diffusion_solver = new ColoredGaussSeidelDiffuse(settings);
    } else if (diffusion_type == DiffusionMethods::Explicit) {
        *diffusion_solver = new ExplicitDiffuse(settings);
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
void SetPressureSolver(Settings::Settings const &settings,
                       IPressure **pressure_solver,
                       const std::string &pressure_type) {
    if (pressure_type == PressureMethods::VCycleMG) {
        *pressure_solver = new VCycleMG(settings);
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
/// \param  source_solver Pointer to SourceSolver
/// \param  source_type Name of SourceSolver
// ***************************************************************************************
void SetSourceSolver(Settings::Settings const &settings,
                     ISource **source_solver,
                     const std::string& source_type) {
    if (source_type == SourceMethods::ExplicitEuler) {
        *source_solver = new ExplicitEulerSource(settings);
    } else {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(class_name);
        logger->error("Source method {} not yet implemented! Simulation stopped!", source_type);
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
void SetTurbulenceSolver(Settings::Settings const &settings,
                         ITurbulence **turbulence_solver,
                         const std::string& turbulence_type) {
    if (turbulence_type == TurbulenceMethods::ConstSmagorinsky) {
        *turbulence_solver = new ConstSmagorinsky(settings);
    } else if (turbulence_type == TurbulenceMethods::DynamicSmagorinsky) {
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

    // usuable for concentration and temperature source
    void set_source_function(const Settings::Settings &settings,
                             ISourceFunction **source_function,
                             const std::string &source_fct) {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(class_name);
        logger->debug("create source function {}", source_fct);
#endif
        if (source_fct == SourceMethods::Gauss) {
            real HRR = settings.get_real("solver/temperature/source/HRR");    // heat release rate in [kW]
            real cp = settings.get_real("solver/temperature/source/cp");        // specific heat capacity in [kJ/ kg K]
            real x0 = settings.get_real("solver/temperature/source/x0");
            real y0 = settings.get_real("solver/temperature/source/y0");
            real z0 = settings.get_real("solver/temperature/source/z0");
            real sigma_x = settings.get_real("solver/temperature/source/sigma_x");
            real sigma_y = settings.get_real("solver/temperature/source/sigma_y");
            real sigma_z = settings.get_real("solver/temperature/source/sigma_z");
            real tau = settings.get_real("solver/temperature/source/tau");
            *source_function = new GaussFunction(HRR, cp, x0, y0, z0, sigma_x, sigma_y, sigma_z, tau);
        } else if (source_fct == SourceMethods::BuoyancyST_MMS) {
            *source_function = new BuoyancyMMS(settings);
        } else if (source_fct == SourceMethods::Cube) {
            real x_start = settings.get_real("solver/temperature/source/x_start");
            real y_start = settings.get_real("solver/temperature/source/y_start");
            real z_start = settings.get_real("solver/temperature/source/z_start");
            real x_end = settings.get_real("solver/temperature/source/x_end");
            real y_end = settings.get_real("solver/temperature/source/y_end");
            real z_end = settings.get_real("solver/temperature/source/z_end");
            real val = settings.get_real("solver/temperature/source/value");
            *source_function = new Cube(val, x_start, y_start, z_start, x_end, y_end, z_end);
        } else if (source_fct == SourceMethods::Zero) {
            *source_function = new Zero();
        } else {
#ifndef BENCHMARKING
            logger->warn("Source method {} not yet implemented!", source_fct);
#endif
        }
        bool has_noise = settings.get_bool("solver/temperature/source/random");
        if (has_noise) {
            real range = settings.get_real("solver/temperature/source/random/range");  // +- range of random numbers
            bool has_custom_seed = settings.get_bool("solver/temperature/source/random/custoseed");
            bool has_custom_steps = settings.get_bool("solver/temperature/source/random/custosteps");

            int seed = -1;
            if (has_custom_seed) {
                seed = settings.get_int("solver/temperature/source/random/seed");
            }

            real step_size = 1.0;
            if (has_custom_steps) {
                step_size = settings.get_real("solver/temperature/source/random/step_size");
            }

            IRandomField *noise_maker = new UniformRandom(range, step_size, seed);
            (*source_function)->set_noise(noise_maker);
        }
    }
}  // namespace SolverSelection
