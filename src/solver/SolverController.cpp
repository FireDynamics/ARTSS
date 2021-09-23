/// \file       SolverController.cpp
/// \brief      Manager class for everything regarding solver
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.

#include "AdvectionSolver.h"
#include "AdvectionDiffusionSolver.h"
#include "DiffusionSolver.h"
#include "DiffusionTurbSolver.h"
#include "NSTurbSolver.h"
#include "NSSolver.h"
#include "NSTempSolver.h"
#include "NSTempTurbSolver.h"
#include "NSTempConSolver.h"
#include "NSTempTurbConSolver.h"
#include "PressureSolver.h"
#include "SolverController.h"
#include "SolverSelection.h"
#include "../boundary/BoundaryController.h"
#include "../Domain.h"
#include "../Functions.h"
#include "../source/GaussFunction.h"
#include "../source/BuoyancyMMS.h"
#include "../source/Cube.h"
#include "../source/ExplicitEulerSource.h"
#include "../source/Zero.h"
#include "../utility/Parameters.h"
#include "../randomField/UniformRandom.h"

SolverController::SolverController() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    auto params = Parameters::getInstance();
    m_field_controller = new FieldController();
    std::string string_solver = params->get("solver/description");
    init_solver(string_solver);
#ifndef BENCHMARKING
    m_logger->info("Start initialising....");
#endif
    set_up_sources();
    set_up_fields(string_solver);
    // TODO unclean, first updating device to apply boundary and then updating host to create temporary fields.
    m_field_controller->update_device();
#ifndef BENCHMARKING
    m_logger->debug("set up boundary");
#endif
    m_field_controller->set_up_boundary();
    m_field_controller->update_host();
    m_field_controller->update_data();

    source_velocity = nullptr;
    source_temperature = nullptr;
    source_concentration = nullptr;
}

SolverController::~SolverController() {
    delete m_field_controller;
    delete m_solver;
    delete source_velocity;
    delete source_temperature;
    delete source_concentration;
}

void SolverController::set_up_sources() {
    auto params = Parameters::getInstance();
    // source of temperature
    if (m_has_temperature) {
        // source
        std::string temp_type = params->get("solver/temperature/source/type");
        if (temp_type == SourceMethods::ExplicitEuler) {
            source_temperature = new ExplicitEulerSource();
        } else {
#ifndef BENCHMARKING
            m_logger->critical("Source type {} not yet implemented! Simulation stopped!", temp_type);
#endif
            std::exit(1);
            // TODO Error handling
        }
        // temperature function
        std::string temp_fct = params->get("solver/temperature/source/temp_fct");
        if (temp_fct == SourceMethods::GaussST) {
            real HRR = params->get_real("solver/temperature/source/HRR");    // heat release rate in [kW]
            real cp = params->get_real("solver/temperature/source/cp");        // specific heat capacity in [kJ/ kg K]
            real x0 = params->get_real("solver/temperature/source/x0");
            real y0 = params->get_real("solver/temperature/source/y0");
            real z0 = params->get_real("solver/temperature/source/z0");
            real sigma_x = params->get_real("solver/temperature/source/sigma_x");
            real sigma_y = params->get_real("solver/temperature/source/sigma_y");
            real sigma_z = params->get_real("solver/temperature/source/sigma_z");
            real tau = params->get_real("solver/temperature/source/tau");
            m_source_function_temperature = new GaussFunction(HRR, cp, x0, y0, z0, sigma_x, sigma_y, sigma_z, tau);
        } else if (temp_fct == SourceMethods::BuoyancyST_MMS) {
            m_source_function_temperature = new BuoyancyMMS();
        } else if (temp_fct == SourceMethods::Cube) {
            real x_start = params->get_real("solver/temperature/source/x_start");
            real y_start = params->get_real("solver/temperature/source/y_start");
            real z_start = params->get_real("solver/temperature/source/z_start");
            real x_end = params->get_real("solver/temperature/source/x_end");
            real y_end = params->get_real("solver/temperature/source/y_end");
            real z_end = params->get_real("solver/temperature/source/z_end");
            real val = params->get_real("solver/temperature/source/value");
            m_source_function_temperature = new Cube(val, x_start, y_start, z_start, x_end, y_end, z_end);
        } else if (temp_fct == SourceMethods::Zero) {
            m_source_function_temperature = new Zero();
        } else {
#ifndef BENCHMARKING
            m_logger->warn("Source method {} not yet implemented!", temp_fct);
#endif
        }
        bool has_noise = params->get("solver/temperature/source/random") == XML_TRUE;
        if (has_noise) {
            real range = params->get_real("solver/temperature/source/random/range");  // +- range of random numbers
            bool has_custom_seed = params->get("solver/temperature/source/random/custom_seed") == XML_TRUE;
            bool has_custom_steps = params->get("solver/temperature/source/random/custom_steps") == XML_TRUE;

            int seed = -1;
            if (has_custom_seed) {
                seed = params->get_int("solver/temperature/source/random/seed");
            }

            real step_size = 1.0;
            if (has_custom_steps) {
                step_size = params->get_real("solver/temperature/source/random/step_size");
            }

            IRandomField *noise_maker = new UniformRandom(range, step_size, seed);
            m_source_function_temperature->set_noise(noise_maker);
        }
    }

    if (m_has_concentration) {
        // Source of concentration
        std::string con_type = params->get("solver/concentration/source/type");
        if (con_type == SourceMethods::ExplicitEuler) {
            source_concentration = new ExplicitEulerSource();
        } else {
#ifndef BENCHMARKING
            m_logger->critical("Source type {} not yet implemented! Simulation stopped!", con_type);
#endif
            std::exit(1);
            // TODO Error handling
        }
        // concentration function
        std::string con_fct = params->get("solver/concentration/source/con_fct");
        if (con_fct == SourceMethods::GaussSC) {
            // get parameters for Gauss function
            real HRR = params->get_real("solver/concentration/source/HRR");       // heat release rate in [kW]
            real Hc = params->get_real("solver/concentration/source/Hc");        // heating value in [kJ/kg]
            real Ys = params->get_real("solver/concentration/source/Ys");        // soot yield in [g/g]
            real YsHRR = Ys * HRR;
            real x0 = params->get_real("solver/concentration/source/x0");
            real y0 = params->get_real("solver/concentration/source/y0");
            real z0 = params->get_real("solver/concentration/source/z0");
            real sigma_x = params->get_real("solver/concentration/source/sigma_x");
            real sigma_y = params->get_real("solver/concentration/source/sigma_y");
            real sigma_z = params->get_real("solver/concentration/source/sigma_z");
            real tau = params->get_real("solver/concentration/source/tau");

            m_source_function_concentration = new GaussFunction(YsHRR, Hc, x0, y0, z0, sigma_x, sigma_y, sigma_z, tau);
        } else if (con_fct == SourceMethods::Zero) {
            m_source_function_temperature = new Zero();
        } else {
#ifndef BENCHMARKING
            m_logger->warn("Source method {} not yet implemented!", con_fct);
#endif
        }
    }

    // Source term for momentum
    if (m_has_momentum_source) {
        std::string source_type = params->get("solver/source/type");
        if (source_type == SourceMethods::ExplicitEuler) {
            source_velocity = new ExplicitEulerSource();
        } else {
#ifndef BENCHMARKING
            m_logger->critical("Source function {} not yet implemented! Simulation stopped!", source_type);
#endif
            std::exit(1);
            // TODO(issue 6) Error handling
        }
    }
}

void SolverController::init_solver(const std::string &string_solver) {
    if (string_solver == SolverTypes::AdvectionSolver) {
        m_solver = new AdvectionSolver(m_field_controller);
    } else if (string_solver == SolverTypes::AdvectionDiffusionSolver) {
        m_solver = new AdvectionDiffusionSolver(m_field_controller);
    } else if (string_solver == SolverTypes::DiffusionSolver) {
        m_solver = new DiffusionSolver(m_field_controller);
    } else if (string_solver == SolverTypes::DiffusionTurbSolver) {
        m_solver = new DiffusionTurbSolver(m_field_controller);
    } else if (string_solver == SolverTypes::NSSolver) {
        m_solver = new NSSolver(m_field_controller);
        m_has_momentum_source = true;
    } else if (string_solver == SolverTypes::NSTurbSolver) {
        m_solver = new NSTurbSolver(m_field_controller);
        m_has_momentum_source = true;
        m_has_turbulence = true;
    } else if (string_solver == SolverTypes::NSTempSolver) {
        m_solver = new NSTempSolver(m_field_controller);
        m_has_momentum_source = true;
        m_has_temperature = true;
    } else if (string_solver == SolverTypes::NSTempTurbSolver) {
        m_solver = new NSTempTurbSolver(m_field_controller);
        m_has_momentum_source = true;
        m_has_temperature = true;
        m_has_turbulence = true;
    } else if (string_solver == SolverTypes::NSTempConSolver) {
        m_solver = new NSTempConSolver(m_field_controller);
        m_has_momentum_source = true;
        m_has_temperature = true;
        m_has_concentration = true;
    } else if (string_solver == SolverTypes::NSTempTurbConSolver) {
        m_solver = new NSTempTurbConSolver(m_field_controller);
        m_has_momentum_source = true;
        m_has_temperature = true;
        m_has_concentration = true;
        m_has_turbulence = true;
    } else if (string_solver == SolverTypes::PressureSolver) {
        m_solver = new PressureSolver(m_field_controller);
    } else {
#ifndef BENCHMARKING
        m_logger->error("Solver {} not yet implemented! Simulation stopped!", string_solver);
#else
        std::cout << "Solver not yet implemented! Simulation stopped!" << std::endl;
        std::flush(std::cout);
#endif
        std::exit(1);
        // TODO Error handling
    }
}

void SolverController::set_up_fields(const std::string &string_solver) {
    auto params = Parameters::getInstance();
    std::string string_init_usr_fct = params->get("initial_conditions/usr_fct");
    bool random = params->get("initial_conditions/random") == XML_TRUE;

    if (string_init_usr_fct == FunctionNames::GaussBubble) {
        if (string_solver == SolverTypes::AdvectionSolver) {
            Functions::GaussBubble(m_field_controller->get_field_u(), 0.);
            Functions::GaussBubble(m_field_controller->get_field_v(), 0.);
            Functions::GaussBubble(m_field_controller->get_field_w(), 0.);
        }
    } else if (string_init_usr_fct == FunctionNames::Drift) {
        if (string_solver == SolverTypes::AdvectionSolver ||
            string_solver == SolverTypes::NSSolver ||
            string_solver == SolverTypes::NSTempSolver ||
            string_solver == SolverTypes::NSTempConSolver ||
            string_solver == SolverTypes::NSTempTurbConSolver ||
            string_solver == SolverTypes::NSTempTurbSolver ||
            string_solver == SolverTypes::NSTurbSolver) {
            Functions::Drift(m_field_controller->get_field_u(),
                             m_field_controller->get_field_v(),
                             m_field_controller->get_field_w(),
                             m_field_controller->get_field_p());
        }
        if (string_solver == SolverTypes::NSTempSolver ||
            string_solver == SolverTypes::NSTempConSolver ||
            string_solver == SolverTypes::NSTempTurbConSolver ||
            string_solver == SolverTypes::NSTempTurbSolver) {
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::ExpSinusProd) {
        // Diffusion test case
        if (string_solver == SolverTypes::DiffusionSolver ||
            string_solver == SolverTypes::DiffusionTurbSolver) {
            Functions::ExpSinusProd(m_field_controller->get_field_u(), 0.);
            Functions::ExpSinusProd(m_field_controller->get_field_v(), 0.);
            Functions::ExpSinusProd(m_field_controller->get_field_w(), 0.);
        }
    } else if (string_init_usr_fct == FunctionNames::Hat) {
        if (string_solver == SolverTypes::DiffusionSolver ||
            string_solver == SolverTypes::DiffusionTurbSolver) {
            Functions::Hat(m_field_controller->get_field_u());
            Functions::Hat(m_field_controller->get_field_v());
            Functions::Hat(m_field_controller->get_field_w());
        }
    } else if (string_init_usr_fct == FunctionNames::ExpSinusSum) {
        // Burgers (=nonlinear Advection + Diffusion) test case
        if (string_solver == SolverTypes::AdvectionDiffusionSolver) {
            Functions::ExpSinusSum(m_field_controller->get_field_u(),
                                   m_field_controller->get_field_v(),
                                   m_field_controller->get_field_w(),
                                   0.);
        }
    } else if (string_init_usr_fct == FunctionNames::SinSinSin) {
        if (string_solver == SolverTypes::PressureSolver) {
            // Pressure test case
            Functions::SinSinSin(m_field_controller->get_field_rhs());
        }
    } else if (string_init_usr_fct == FunctionNames::McDermott) {
        // NavierStokes test case: McDermott (no force, no temperature) 2D
        if (string_solver == SolverTypes::NSSolver or \
            string_solver == SolverTypes::NSTurbSolver || \
            string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            Functions::McDermott(m_field_controller->get_field_u(),
                                 m_field_controller->get_field_v(),
                                 m_field_controller->get_field_w(),
                                 m_field_controller->get_field_p(),
                                 0.);
            m_field_controller->get_field_p().set_value(0.);
        }
        if (string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::Vortex) {
        // NavierStokes test case: Vortex (no force, no temperature) 2D
        if (string_solver == SolverTypes::NSSolver || \
            string_solver == SolverTypes::NSTurbSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver || \
            string_solver == SolverTypes::NSTempSolver) {
            Functions::Vortex(m_field_controller->get_field_u(),
                              m_field_controller->get_field_v(),
                              m_field_controller->get_field_w(),
                              m_field_controller->get_field_p());
            m_field_controller->get_field_p().set_value(0.);
        }
        if (string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::VortexY) {
        // NavierStokes test case: Vortex (no force, no temperature) 2D
        if (string_solver == SolverTypes::NSSolver || \
            string_solver == SolverTypes::NSTurbSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver || \
            string_solver == SolverTypes::NSTempSolver) {
            Functions::VortexY(m_field_controller->get_field_u(),
                               m_field_controller->get_field_v(),
                               m_field_controller->get_field_w(),
                               m_field_controller->get_field_p());
            m_field_controller->get_field_p().set_value(0.);
        }
        if (string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::Beltrami) {
        // NavierStokes test case: Beltrami  (no force, no temperature) 3D
        if (string_solver == SolverTypes::NSSolver || \
            string_solver == SolverTypes::NSTurbSolver) {
            Functions::Beltrami(m_field_controller->get_field_u(),
                                m_field_controller->get_field_v(),
                                m_field_controller->get_field_w(),
                                m_field_controller->get_field_p(),
                                0.);
            // m_field_controller->field_p->set_value(0.);
        }
    } else if (string_init_usr_fct == FunctionNames::BuoyancyMMS) {
        // NavierStokesTemp test case
        // MMS
        if (string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            Functions::BuoyancyMMS(m_field_controller->get_field_u(),
                                   m_field_controller->get_field_v(),
                                   m_field_controller->get_field_w(),
                                   m_field_controller->get_field_p(),
                                   m_field_controller->get_field_T(),
                                   0.);
            m_field_controller->get_field_p().set_value(0.);
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::Uniform) {
        // Uniform Temperature unequal to zero
        if (string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            real val = params->get_real("initial_conditions/val");
            Functions::Uniform(m_field_controller->get_field_T(), val);
            if (random) {
                call_random(m_field_controller->get_field_T());
            }
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == "LayersT") {
        if (string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            Functions::Layers(m_field_controller->get_field_T());
            if (random) {
                call_random(m_field_controller->get_field_T());
            }
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::Zero) {
        // NavierStokes test case: Channel Flow (with uniform force in x-direction)
        if ((string_solver == SolverTypes::NSSolver || string_solver == SolverTypes::NSTurbSolver)
            && params->get("solver/source/force_fct") == FunctionNames::Uniform) {
            real val_x = params->get_real("solver/source/val_x");
            real val_y = params->get_real("solver/source/val_y");
            real val_z = params->get_real("solver/source/val_z");
            std::string dir = params->get("solver/source/dir");

            if (dir.find('x') != std::string::npos) {
                Functions::Uniform(m_field_controller->get_field_force_x(), val_x);
            }
            if (dir.find('y') != std::string::npos) {
                Functions::Uniform(m_field_controller->get_field_force_y(), val_y);
            }
            if (dir.find('z') != std::string::npos) {
                Functions::Uniform(m_field_controller->get_field_force_z(), val_z);
            }
        } else {
#ifndef BENCHMARKING
            m_logger->info("Initial values all set to zero!");
#endif
        }
        // Random concentration
        if ((string_solver == SolverTypes::NSTempConSolver ||
             string_solver == SolverTypes::NSTempTurbConSolver)
            && params->get("initial_conditions/con_fct") == FunctionNames::RandomC) {
            real Ca = params->get_real("initial_conditions/Ca");        // ambient concentration
            Functions::Uniform(m_field_controller->get_field_concentration(), Ca);
            call_random(m_field_controller->get_field_concentration());
        }
    } else if (string_init_usr_fct == FunctionNames::Jet) {
        std::string dir = params->get("initial_conditions/dir");
        real value = params->get_real("initial_conditions/value");
        auto domain = Domain::getInstance();
        if (dir == "x") {
            real y1 = params->get_real("initial_conditions/y1");
            real y2 = params->get_real("initial_conditions/y2");
            real z1 = params->get_real("initial_conditions/z1");
            real z2 = params->get_real("initial_conditions/z2");
            size_t index_x1 = domain->get_index_x1();
            size_t index_x2 = domain->get_index_x2();
            size_t index_y1 = Utility::get_index(y1, domain->get_dy(), domain->get_Y1());
            size_t index_y2 = Utility::get_index(y2, domain->get_dy(), domain->get_Y1());
            size_t index_z1 = Utility::get_index(z1, domain->get_dz(), domain->get_Z1());
            size_t index_z2 = Utility::get_index(z2, domain->get_dz(), domain->get_Z1());
            Functions::Jet(m_field_controller->get_field_u(),
                           index_x1, index_x2,
                           index_y1, index_y2,
                           index_z1, index_z2,
                           value);
            if (random) {
                call_random(m_field_controller->get_field_u());
            }
        } else if (dir == "y") {
            real x1 = params->get_real("initial_conditions/x1");
            real x2 = params->get_real("initial_conditions/x2");
            real z1 = params->get_real("initial_conditions/z1");
            real z2 = params->get_real("initial_conditions/z2");
            size_t index_x1 = Utility::get_index(x1, domain->get_dx(), domain->get_X1());
            size_t index_x2 = Utility::get_index(x2, domain->get_dx(), domain->get_X1());
            size_t index_y1 = domain->get_index_y1();
            size_t index_y2 = domain->get_index_y2();
            size_t index_z1 = Utility::get_index(z1, domain->get_dz(), domain->get_Z1());
            size_t index_z2 = Utility::get_index(z2, domain->get_dz(), domain->get_Z1());
            Functions::Jet(m_field_controller->get_field_v(),
                           index_x1, index_x2,
                           index_y1, index_y2,
                           index_z1, index_z2,
                           value);
            if (random) {
                call_random(m_field_controller->get_field_v());
            }
        } else if (dir == "z") {
            real x1 = params->get_real("initial_conditions/x1");
            real x2 = params->get_real("initial_conditions/x2");
            real y1 = params->get_real("initial_conditions/y1");
            real y2 = params->get_real("initial_conditions/y2");
            size_t index_x1 = Utility::get_index(x1, domain->get_dx(), domain->get_X1());
            size_t index_x2 = Utility::get_index(x2, domain->get_dx(), domain->get_X1());
            size_t index_y1 = Utility::get_index(y1, domain->get_dy(), domain->get_Y1());
            size_t index_y2 = Utility::get_index(y2, domain->get_dy(), domain->get_Y1());
            size_t index_z1 = domain->get_index_z1();
            size_t index_z2 = domain->get_index_z2();
            Functions::Jet(m_field_controller->get_field_w(),
                           index_x1, index_x2,
                           index_y1, index_y2,
                           index_z1, index_z2,
                           value);
            if (random) {
                call_random(m_field_controller->get_field_w());
            }
        }
    } else {
#ifndef BENCHMARKING
        m_logger->warn("Unknown initial function {}", string_init_usr_fct);
        m_logger->info("Initial values all set to zero!");
#endif
    }

    // Sight of boundaries
    auto boundary = BoundaryController::getInstance();
    size_t *inner_list = boundary->get_inner_list_level_joined();
    size_t size_inner_list = boundary->get_size_inner_list();

    Field sight = m_field_controller->get_field_sight();
    for (size_t i = 0; i < size_inner_list; i++) {
        size_t idx = inner_list[i];
        sight[idx] = 0.;
    }
}

//======================================= read and call random function ==================================
// ***************************************************************************************
/// \brief  Calls random function and reads necessary input variables
/// \param  field       field as a pointer
// ***************************************************************************************
void SolverController::call_random(Field &field) {
    auto params = Parameters::getInstance();
    real range = params->get_real("initial_conditions/random/range");  // +- range of random numbers
    bool is_absolute = params->get("initial_conditions/random/absolute") == XML_TRUE;
    bool has_custom_seed = params->get("initial_conditions/random/custom_seed") == XML_TRUE;
    bool has_custom_steps = params->get("initial_conditions/random/custom_steps") == XML_TRUE;

    int seed = -1;
    if (has_custom_seed) {
        seed = params->get_int("initial_conditions/random/seed");
    }

    real step_size = 1.0;
    if (has_custom_steps) {
        step_size = params->get_real("initial_conditions/random/step_size");
    }

    Functions::Random(field, range, is_absolute, seed, step_size);
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters temperature source functions
// ***************************************************************************************
void SolverController::temperature_source() {
// Temperature source
    if (Parameters::getInstance()->get("solver/temperature/source/temp_fct") == FunctionNames::BuoyancyST_MMS) {
        Functions::BuoyancyST_MMS(m_field_controller->get_field_source_T(), 0.);
    }
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters force source functions
// ***************************************************************************************
void SolverController::force_source() {
    auto params = Parameters::getInstance();
    // Force
    if (params->get("solver/source/force_fct") == SourceMethods::Buoyancy) {
        std::string dir = params->get("solver/source/dir");
        if (params->get("solver/source/use_init_values") == XML_FALSE) {
            real ambient_temperature_value = params->get_real("solver/source/ambient_temperature_value");
            //m_field_controller->get_field_T().set_value(ambient_temperature_value);
            m_field_controller->get_field_T_ambient().set_value(ambient_temperature_value);
        }

        if (dir.find('x') != std::string::npos) {
            Functions::BuoyancyForce(m_field_controller->get_field_force_x(),
                                     m_field_controller->get_field_T(),
                                     m_field_controller->get_field_T_ambient());
        }
        if (dir.find('y') != std::string::npos) {
            Functions::BuoyancyForce(m_field_controller->get_field_force_y(),
                                     m_field_controller->get_field_T(),
                                     m_field_controller->get_field_T_ambient());
        }
        if (dir.find('z') != std::string::npos) {
            Functions::BuoyancyForce(m_field_controller->get_field_force_z(),
                                     m_field_controller->get_field_T(),
                                     m_field_controller->get_field_T_ambient());
        }
    }
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters momentum source functions
// ***************************************************************************************
void SolverController::momentum_source() {
    // Momentum source
    auto params = Parameters::getInstance();
    std::string dir_vel = params->get("solver/source/dir");
    if (dir_vel.find('x') != std::string::npos) {
        source_velocity->buoyancy_force(m_field_controller->get_field_force_x(),
                                        m_field_controller->get_field_T(),
                                        m_field_controller->get_field_T_ambient());
    }
    if (dir_vel.find('y') != std::string::npos) {
        source_velocity->buoyancy_force(m_field_controller->get_field_force_y(),
                                        m_field_controller->get_field_T(),
                                        m_field_controller->get_field_T_ambient());
    }
    if (dir_vel.find('z') != std::string::npos) {
        source_velocity->buoyancy_force(m_field_controller->get_field_force_z(),
                                        m_field_controller->get_field_T(),
                                        m_field_controller->get_field_T_ambient());
    }
}


//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters such as momentum/temperature/concentration source functions
/// \param  t   time
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void SolverController::update_sources(real t_cur, bool sync) {
    auto params = Parameters::getInstance();

    // Momentum source
    if (m_has_momentum_source) {
        std::string forceFct = params->get("solver/source/force_fct");
        if (forceFct == SourceMethods::Zero || \
            forceFct == SourceMethods::Uniform) {
        } else if (forceFct == SourceMethods::Buoyancy) {
#ifndef BENCHMARKING
            m_logger->info("Update f(T) ...");
#endif
            if (params->get("solver/source/use_init_values") == XML_FALSE) {
                int ambient_temperature_value = params->get_int("solver/source/ambient_temperature_value");
                m_field_controller->get_field_T_ambient().set_value(ambient_temperature_value);
            } else {
                m_field_controller->get_field_T_ambient().copy_data(m_field_controller->get_field_T());
            }
            momentum_source();
        } else {
#ifndef BENCHMARKING
            m_logger->critical("Source function not yet implemented! Simulation stopped!");
#endif
            std::exit(1);
            // TODO(issue 6) Error handling
        }
    }

    // Temperature source
    if (m_has_temperature) {
        m_source_function_temperature->update_source(m_field_controller->get_field_source_T(), t_cur);
        std::string tempFct = params->get("solver/temperature/source/temp_fct");
    }

    // Concentration source
    if (m_has_concentration) {
        m_source_function_concentration->update_source(m_field_controller->get_field_source_concentration(), t_cur);
    }
    if (sync) {
#pragma acc wait
    }
}

void SolverController::solver_do_step(real t, bool sync) {
    m_solver->do_step(t, sync);
}
