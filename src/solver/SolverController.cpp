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
#include "../DomainData.h"
#include "../Functions.h"
#include "../source/GaussFunction.h"
#include "../source/BuoyancyMMS.h"
#include "../source/Cube.h"
#include "../source/ExplicitEulerSource.h"
#include "../source/Zero.h"
#include "../randomField/UniformRandom.h"


SolverController::SolverController(Settings::Settings const &settings) :
        m_settings(settings) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = new FieldController();
    std::string string_solver = m_settings.get("solver/description");
    init_solver(string_solver);
#ifndef BENCHMARKING
    m_logger->info("Start initialising....");
#endif
    set_up_sources();
    set_up_fields(string_solver);
#ifndef BENCHMARKING
    m_logger->debug("set up boundary");
#endif
    m_field_controller->set_up_boundary();
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
#ifndef BENCHMARKING
    m_logger->debug("set up sources");
#endif
    // source of temperature
    if (m_has_temperature) {
        // source
        std::string temp_type = m_settings.get("solver/temperature/source/type");
#ifndef BENCHMARKING
        m_logger->debug("create temperature type function {}", temp_type);
#endif
        if (temp_type == SourceMethods::ExplicitEuler) {
            source_temperature = new ExplicitEulerSource(m_settings);
        } else {
#ifndef BENCHMARKING
            m_logger->critical("Source type {} not yet implemented! Simulation stopped!", temp_type);
#endif
            std::exit(1);
            // TODO Error handling
        }
        // temperature function
        std::string temp_fct = m_settings.get("solver/temperature/source/temp_fct");
#ifndef BENCHMARKING
        m_logger->debug("create temperature source function {}", temp_fct);
#endif
        if (temp_fct == SourceMethods::GaussST) {
            real HRR = m_settings.get_real("solver/temperature/source/HRR");    // heat release rate in [kW]
            real cp = m_settings.get_real("solver/temperature/source/cp");        // specific heat capacity in [kJ/ kg K]
            real x0 = m_settings.get_real("solver/temperature/source/x0");
            real y0 = m_settings.get_real("solver/temperature/source/y0");
            real z0 = m_settings.get_real("solver/temperature/source/z0");
            real sigma_x = m_settings.get_real("solver/temperature/source/sigma_x");
            real sigma_y = m_settings.get_real("solver/temperature/source/sigma_y");
            real sigma_z = m_settings.get_real("solver/temperature/source/sigma_z");
            real tau = m_settings.get_real("solver/temperature/source/tau");
            m_source_function_temperature = new GaussFunction(HRR, cp, x0, y0, z0, sigma_x, sigma_y, sigma_z, tau);
        } else if (temp_fct == SourceMethods::BuoyancyST_MMS) {
            m_source_function_temperature = new BuoyancyMMS(m_settings);
        } else if (temp_fct == SourceMethods::Cube) {
            real x_start = m_settings.get_real("solver/temperature/source/x_start");
            real y_start = m_settings.get_real("solver/temperature/source/y_start");
            real z_start = m_settings.get_real("solver/temperature/source/z_start");
            real x_end = m_settings.get_real("solver/temperature/source/x_end");
            real y_end = m_settings.get_real("solver/temperature/source/y_end");
            real z_end = m_settings.get_real("solver/temperature/source/z_end");
            real val = m_settings.get_real("solver/temperature/source/value");
            m_source_function_temperature = new Cube(val, x_start, y_start, z_start, x_end, y_end, z_end);
        } else if (temp_fct == SourceMethods::Zero) {
            m_source_function_temperature = new Zero();
        } else {
#ifndef BENCHMARKING
            m_logger->warn("Source method {} not yet implemented!", temp_fct);
#endif
        }
        bool has_noise = m_settings.get_bool("solver/temperature/source/random");
        if (has_noise) {
            real range = m_settings.get_real("solver/temperature/source/random/range");  // +- range of random numbers
            bool has_custom_seed = m_settings.get_bool("solver/temperature/source/random/custom_seed");
            bool has_custom_steps = m_settings.get_bool("solver/temperature/source/random/custom_steps");

            int seed = -1;
            if (has_custom_seed) {
                seed = m_settings.get_int("solver/temperature/source/random/seed");
            }

            real step_size = 1.0;
            if (has_custom_steps) {
                step_size = m_settings.get_real("solver/temperature/source/random/step_size");
            }

            IRandomField *noise_maker = new UniformRandom(range, step_size, seed);
            m_source_function_temperature->set_noise(noise_maker);
        }
    }

    if (m_has_concentration) {
        // Source of concentration
        std::string con_type = m_settings.get("solver/concentration/source/type");
        if (con_type == SourceMethods::ExplicitEuler) {
            source_concentration = new ExplicitEulerSource(m_settings);
        } else {
#ifndef BENCHMARKING
            m_logger->critical("Source type {} not yet implemented! Simulation stopped!", con_type);
#endif
            std::exit(1);
            // TODO Error handling
        }
        // concentration function
        std::string con_fct = m_settings.get("solver/concentration/source/con_fct");
        if (con_fct == SourceMethods::GaussSC) {
            // get parameters for Gauss function
            real HRR = m_settings.get_real("solver/concentration/source/HRR");       // heat release rate in [kW]
            real Hc = m_settings.get_real("solver/concentration/source/Hc");        // heating value in [kJ/kg]
            real Ys = m_settings.get_real("solver/concentration/source/Ys");        // soot yield in [g/g]
            real YsHRR = Ys * HRR;
            real x0 = m_settings.get_real("solver/concentration/source/x0");
            real y0 = m_settings.get_real("solver/concentration/source/y0");
            real z0 = m_settings.get_real("solver/concentration/source/z0");
            real sigma_x = m_settings.get_real("solver/concentration/source/sigma_x");
            real sigma_y = m_settings.get_real("solver/concentration/source/sigma_y");
            real sigma_z = m_settings.get_real("solver/concentration/source/sigma_z");
            real tau = m_settings.get_real("solver/concentration/source/tau");

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
        std::string source_type = m_settings.get("solver/source/type");
        if (source_type == SourceMethods::ExplicitEuler) {
            source_velocity = new ExplicitEulerSource(m_settings);
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
#ifndef BENCHMARKING
    m_logger->debug("initialise solver {}", string_solver);
#endif
    if (string_solver == SolverTypes::AdvectionSolver) {
        m_solver = new AdvectionSolver(m_settings, m_field_controller);
    } else if (string_solver == SolverTypes::AdvectionDiffusionSolver) {
        m_solver = new AdvectionDiffusionSolver(m_settings, m_field_controller);
    } else if (string_solver == SolverTypes::DiffusionSolver) {
        m_solver = new DiffusionSolver(m_settings, m_field_controller);
    } else if (string_solver == SolverTypes::DiffusionTurbSolver) {
        m_solver = new DiffusionTurbSolver(m_settings, m_field_controller);
    } else if (string_solver == SolverTypes::NSSolver) {
        m_solver = new NSSolver(m_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (string_solver == SolverTypes::NSTurbSolver) {
        m_solver = new NSTurbSolver(m_settings, m_field_controller);
        m_has_momentum_source = true;
        m_has_turbulence = true;
    } else if (string_solver == SolverTypes::NSTempSolver) {
        m_solver = new NSTempSolver(m_settings, m_field_controller);
        m_has_momentum_source = true;
        m_has_temperature = true;
    } else if (string_solver == SolverTypes::NSTempTurbSolver) {
        m_solver = new NSTempTurbSolver(m_settings, m_field_controller);
        m_has_momentum_source = true;
        m_has_temperature = true;
        m_has_turbulence = true;
    } else if (string_solver == SolverTypes::NSTempConSolver) {
        m_solver = new NSTempConSolver(m_settings, m_field_controller);
        m_has_momentum_source = true;
        m_has_temperature = true;
        m_has_concentration = true;
    } else if (string_solver == SolverTypes::NSTempTurbConSolver) {
        m_solver = new NSTempTurbConSolver(m_settings, m_field_controller);
        m_has_momentum_source = true;
        m_has_temperature = true;
        m_has_concentration = true;
        m_has_turbulence = true;
    } else if (string_solver == SolverTypes::PressureSolver) {
        m_solver = new PressureSolver(m_settings, m_field_controller);
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
    std::string string_init_usr_fct = m_settings.get("initial_conditions/usr_fct");
    bool random = m_settings.get_bool("initial_conditions/random");

    if (string_init_usr_fct == FunctionNames::gauss_bubble) {
        if (string_solver == SolverTypes::AdvectionSolver) {
            real u_lin = m_settings.get_real("initial_conditions/u_lin");
            real v_lin = m_settings.get_real("initial_conditions/v_lin");
            real w_lin = m_settings.get_real("initial_conditions/w_lin");
            real x_shift = m_settings.get_real("initial_conditions/x_shift");
            real y_shift = m_settings.get_real("initial_conditions/y_shift");
            real z_shift = m_settings.get_real("initial_conditions/z_shift");
            real l = m_settings.get_real("initial_conditions/l");

            Functions::gauss_bubble(m_field_controller->get_field_u(), 0.,
                    u_lin, v_lin, w_lin, x_shift, y_shift, z_shift, l);
            Functions::gauss_bubble(m_field_controller->get_field_v(), 0.,
                    u_lin, v_lin, w_lin, x_shift, y_shift, z_shift, l);
            Functions::gauss_bubble(m_field_controller->get_field_w(), 0.,
                    u_lin, v_lin, w_lin, x_shift, y_shift, z_shift, l);
        }
    } else if (string_init_usr_fct == FunctionNames::drift) {
        if (string_solver == SolverTypes::AdvectionSolver ||
            string_solver == SolverTypes::NSSolver ||
            string_solver == SolverTypes::NSTempSolver ||
            string_solver == SolverTypes::NSTempConSolver ||
            string_solver == SolverTypes::NSTempTurbConSolver ||
            string_solver == SolverTypes::NSTempTurbSolver ||
            string_solver == SolverTypes::NSTurbSolver) {

            real u_lin = m_settings.get_real("initial_conditions/u_lin");
            real v_lin = m_settings.get_real("initial_conditions/v_lin");
            real w_lin = m_settings.get_real("initial_conditions/w_lin");
            real pa = m_settings.get_real("initial_conditions/pa");

            Functions::drift(m_field_controller->get_field_u(),
                             m_field_controller->get_field_v(),
                             m_field_controller->get_field_w(),
                             m_field_controller->get_field_p(),
                             u_lin, v_lin, w_lin, pa);
        }
        if (string_solver == SolverTypes::NSTempSolver ||
            string_solver == SolverTypes::NSTempConSolver ||
            string_solver == SolverTypes::NSTempTurbConSolver ||
            string_solver == SolverTypes::NSTempTurbSolver) {
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::exp_sinus_prod) {
        // Diffusion test case
        if (string_solver == SolverTypes::DiffusionSolver ||
            string_solver == SolverTypes::DiffusionTurbSolver) {
            real nu = m_settings.get_real("physical_parameters/nu");
            real l = m_settings.get_real("initial_conditions/l");

            Functions::exp_sinus_prod(m_field_controller->get_field_u(), 0., nu, l);
            Functions::exp_sinus_prod(m_field_controller->get_field_v(), 0., nu, l);
            Functions::exp_sinus_prod(m_field_controller->get_field_w(), 0., nu, l);
        }
    } else if (string_init_usr_fct == FunctionNames::hat) {
        if (string_solver == SolverTypes::DiffusionSolver ||
            string_solver == SolverTypes::DiffusionTurbSolver) {
            real start_x = m_settings.get_real("initial_conditions/x1");
            real end_x = m_settings.get_real("initial_conditions/x2");
            real start_y = m_settings.get_real("initial_conditions/y1");
            real end_y = m_settings.get_real("initial_conditions/y2");
            real start_z = m_settings.get_real("initial_conditions/z1");
            real end_z = m_settings.get_real("initial_conditions/z2");
            real val_in = m_settings.get_real("initial_conditions/val_in");
            real val_out = m_settings.get_real("initial_conditions/val_out");

            Functions::hat(m_field_controller->get_field_u(),
                    start_x, end_x, start_y, end_y, start_z, end_z, val_in, val_out);
            Functions::hat(m_field_controller->get_field_v(),
                    start_x, end_x, start_y, end_y, start_z, end_z, val_in, val_out);
            Functions::hat(m_field_controller->get_field_w(),
                    start_x, end_x, start_y, end_y, start_z, end_z, val_in, val_out);
        }
    } else if (string_init_usr_fct == FunctionNames::exp_sinus_sum) {
        // Burgers (=nonlinear Advection + Diffusion) test case
        if (string_solver == SolverTypes::AdvectionDiffusionSolver) {
            real nu = m_settings.get_real("physical_parameters/nu");
            Functions::exp_sinus_sum(m_field_controller->get_field_u(),
                                     m_field_controller->get_field_v(),
                                     m_field_controller->get_field_w(),
                                     0., nu);
        }
    } else if (string_init_usr_fct == FunctionNames::sin_sin_sin) {
        if (string_solver == SolverTypes::PressureSolver) {
            // Pressure test case
            m_field_controller->get_field_p().set_value(0.);
            real l = m_settings.get_real("initial_conditions/l");
            Functions::sin_sin_sin(m_field_controller->get_field_rhs(), l);
        }
    } else if (string_init_usr_fct == FunctionNames::mcdermott) {
        // NavierStokes test case: McDermott (no force, no temperature) 2D
        if (string_solver == SolverTypes::NSSolver or \
            string_solver == SolverTypes::NSTurbSolver || \
            string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            real nu = m_settings.get_real("physical_parameters/nu");
            real A = m_settings.get_real("initial_conditions/A");

            Functions::mcdermott(m_field_controller->get_field_u(),
                                 m_field_controller->get_field_v(),
                                 m_field_controller->get_field_w(),
                                 m_field_controller->get_field_p(),
                                 0., nu, A);
            m_field_controller->get_field_p().set_value(0.);
        }
        if (string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::vortex) {
        // NavierStokes test case: Vortex (no force, no temperature) 2D
        if (string_solver == SolverTypes::NSSolver || \
            string_solver == SolverTypes::NSTurbSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver || \
            string_solver == SolverTypes::NSTempSolver) {
            real u_lin = m_settings.get_real("initial_conditions/u_lin");
            real v_lin = m_settings.get_real("initial_conditions/v_lin");
            real pa = m_settings.get_real("initial_conditions/pa");
            real rhoa = m_settings.get_real("initial_conditions/rhoa");
            Functions::vortex(m_field_controller->get_field_u(),
                              m_field_controller->get_field_v(),
                              m_field_controller->get_field_w(),
                              m_field_controller->get_field_p(),
                              u_lin, v_lin, pa, rhoa);
            m_field_controller->get_field_p().set_value(0.);
        }
        if (string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::vortex_y) {
        // NavierStokes test case: Vortex (no force, no temperature) 2D
        if (string_solver == SolverTypes::NSSolver || \
            string_solver == SolverTypes::NSTurbSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver || \
            string_solver == SolverTypes::NSTempSolver) {
            real u_lin = m_settings.get_real("initial_conditions/u_lin");
            real v_lin = m_settings.get_real("initial_conditions/v_lin");
            real pa = m_settings.get_real("initial_conditions/pa");
            real rhoa = m_settings.get_real("initial_conditions/rhoa");

            Functions::vortex_y(m_field_controller->get_field_u(),
                                m_field_controller->get_field_v(),
                                m_field_controller->get_field_w(),
                                m_field_controller->get_field_p(),
                                u_lin, v_lin, pa, rhoa);
            m_field_controller->get_field_p().set_value(0.);
        }
        if (string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::beltrami) {
        // NavierStokes test case: Beltrami  (no force, no temperature) 3D
        if (string_solver == SolverTypes::NSSolver || \
            string_solver == SolverTypes::NSTurbSolver) {
            real a = m_settings.get_real("initial_conditions/a");  // 0.25 * M_PI;
            real d = m_settings.get_real("initial_conditions/d");  // 0.5 * M_PI;
            real nu = m_settings.get_real("physical_parameters/nu");  // 1;

            Functions::beltrami(m_field_controller->get_field_u(),
                                m_field_controller->get_field_v(),
                                m_field_controller->get_field_w(),
                                m_field_controller->get_field_p(),
                                0., a, d, nu);
            // m_field_controller->field_p->set_value(0.);
        }
    } else if (string_init_usr_fct == FunctionNames::buoyancy_mms) {
        // NavierStokesTemp test case
        // MMS
        if (string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            real g = m_settings.get_real("physical_parameters/g");  // 1;
            real nu = m_settings.get_real("physical_parameters/nu");  // 1;
            real beta = m_settings.get_real("physical_parameters/beta");
            real rhoa = m_settings.get_real("initial_conditions/rhoa");

            Functions::buoyancy_mms(m_field_controller->get_field_u(),
                                    m_field_controller->get_field_v(),
                                    m_field_controller->get_field_w(),
                                    m_field_controller->get_field_p(),
                                    m_field_controller->get_field_T(),
                                    0., nu, beta, g, rhoa);
            m_field_controller->get_field_p().set_value(0.);
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::uniform) {
        // Uniform Temperature unequal to zero
        if (string_solver == SolverTypes::NSTempSolver || \
            string_solver == SolverTypes::NSTempConSolver || \
            string_solver == SolverTypes::NSTempTurbConSolver || \
            string_solver == SolverTypes::NSTempTurbSolver) {
            real val = m_settings.get_real("initial_conditions/val");
            Functions::uniform(m_field_controller->get_field_T(), val);
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
            std::string const log_level;
            std::string const log_file;

            int n_layers = m_settings.get_int("initial_conditions/n_layers");
            std::string const dir = m_settings.get("initial_conditions/dir");
            CoordinateAxis axis = Axis::match_axis(dir);
            real *borders = new real[n_layers + 1];
            real *values = new real[n_layers];

            for (int l = 1; l < n_layers; ++l) {
                std::string val_bord_l = "initial_conditions/border_";
                val_bord_l += std::to_string(l);
                borders[l] = m_settings.get_real(val_bord_l);
            }

            for (int l = 0; l < n_layers; ++l) {
                std::string val_out_l = "initial_conditions/value_";
                val_out_l += std::to_string(l + 1);
                values[l] = m_settings.get_real(val_out_l);
            }

            Functions::layers(m_field_controller->get_field_T(),
                              n_layers, axis, borders, values);

            if (random) {
                call_random(m_field_controller->get_field_T());
            }
            force_source();
            temperature_source();
        }
    } else if (string_init_usr_fct == FunctionNames::zero) {
        // NavierStokes test case: Channel Flow (with uniform force in x-direction)
        if ((string_solver == SolverTypes::NSSolver || string_solver == SolverTypes::NSTurbSolver)
            && m_settings.get("solver/source/force_fct") == FunctionNames::uniform) {
            real val_x = m_settings.get_real("solver/source/val_x");
            real val_y = m_settings.get_real("solver/source/val_y");
            real val_z = m_settings.get_real("solver/source/val_z");
            std::string dir = m_settings.get("solver/source/dir");

            if (dir.find('x') != std::string::npos) {
                Functions::uniform(m_field_controller->get_field_force_x(), val_x);
            }
            if (dir.find('y') != std::string::npos) {
                Functions::uniform(m_field_controller->get_field_force_y(), val_y);
            }
            if (dir.find('z') != std::string::npos) {
                Functions::uniform(m_field_controller->get_field_force_z(), val_z);
            }
        } else {
#ifndef BENCHMARKING
            m_logger->info("Initial values all set to zero!");
#endif
        }
    } else if (string_init_usr_fct == FunctionNames::jet) {
        std::string dir = m_settings.get("initial_conditions/dir");
        real value = m_settings.get_real("initial_conditions/value");
        auto domain = DomainData::getInstance();
        if (dir == "x") {
            real y1 = m_settings.get_real("initial_conditions/y1");
            real y2 = m_settings.get_real("initial_conditions/y2");
            real z1 = m_settings.get_real("initial_conditions/z1");
            real z2 = m_settings.get_real("initial_conditions/z2");
            size_t index_x1 = domain->get_index_x1();
            size_t index_x2 = domain->get_index_x2();
            size_t index_y1 = Utility::get_index(y1, domain->get_dy(), domain->get_Y1());
            size_t index_y2 = Utility::get_index(y2, domain->get_dy(), domain->get_Y1());
            size_t index_z1 = Utility::get_index(z1, domain->get_dz(), domain->get_Z1());
            size_t index_z2 = Utility::get_index(z2, domain->get_dz(), domain->get_Z1());
            Functions::jet(m_field_controller->get_field_u(),
                           index_x1, index_x2,
                           index_y1, index_y2,
                           index_z1, index_z2,
                           value);
            if (random) {
                call_random(m_field_controller->get_field_u());
            }
        } else if (dir == "y") {
            real x1 = m_settings.get_real("initial_conditions/x1");
            real x2 = m_settings.get_real("initial_conditions/x2");
            real z1 = m_settings.get_real("initial_conditions/z1");
            real z2 = m_settings.get_real("initial_conditions/z2");
            size_t index_x1 = Utility::get_index(x1, domain->get_dx(), domain->get_X1());
            size_t index_x2 = Utility::get_index(x2, domain->get_dx(), domain->get_X1());
            size_t index_y1 = domain->get_index_y1();
            size_t index_y2 = domain->get_index_y2();
            size_t index_z1 = Utility::get_index(z1, domain->get_dz(), domain->get_Z1());
            size_t index_z2 = Utility::get_index(z2, domain->get_dz(), domain->get_Z1());
            Functions::jet(m_field_controller->get_field_v(),
                           index_x1, index_x2,
                           index_y1, index_y2,
                           index_z1, index_z2,
                           value);
            if (random) {
                call_random(m_field_controller->get_field_v());
            }
        } else if (dir == "z") {
            real x1 = m_settings.get_real("initial_conditions/x1");
            real x2 = m_settings.get_real("initial_conditions/x2");
            real y1 = m_settings.get_real("initial_conditions/y1");
            real y2 = m_settings.get_real("initial_conditions/y2");
            size_t index_x1 = Utility::get_index(x1, domain->get_dx(), domain->get_X1());
            size_t index_x2 = Utility::get_index(x2, domain->get_dx(), domain->get_X1());
            size_t index_y1 = Utility::get_index(y1, domain->get_dy(), domain->get_Y1());
            size_t index_y2 = Utility::get_index(y2, domain->get_dy(), domain->get_Y1());
            size_t index_z1 = domain->get_index_z1();
            size_t index_z2 = domain->get_index_z2();
            Functions::jet(m_field_controller->get_field_w(),
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
    size_t *domain_list = boundary->get_domain_inner_list_level_joined();
    size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);

    Field &sight = m_field_controller->get_field_sight();
    sight.update_host();
    for (size_t i = 0; i < size_domain_list; i++) {
        size_t idx = domain_list[i];
        sight[idx] = 0.;
    }
}

//======================================= read and call random function ==================================
// ***************************************************************************************
/// \brief  Calls random function and reads necessary input variables
/// \param  field       field as a pointer
// ***************************************************************************************
void SolverController::call_random(Field &field) {
    real range = m_settings.get_real("initial_conditions/random/range");  // +- range of random numbers
    bool is_absolute = m_settings.get_bool("initial_conditions/random/absolute");
    bool has_custom_seed = m_settings.get_bool("initial_conditions/random/custom_seed");
    bool has_custom_steps = m_settings.get_bool("initial_conditions/random/custom_steps");

    int seed = -1;
    if (has_custom_seed) {
        seed = m_settings.get_int("initial_conditions/random/seed");
    }

    real step_size = 1.0;
    if (has_custom_steps) {
        step_size = m_settings.get_real("initial_conditions/random/step_size");
    }

    Functions::random(field, range, is_absolute, seed, step_size);
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters temperature source functions
// ***************************************************************************************
void SolverController::temperature_source() {
// Temperature source
    if (m_settings.get("solver/temperature/source/temp_fct") == FunctionNames::buoyancy_st_mms) {
        real g = m_settings.get_real("physical_parameters/g");  // 1;
        real nu = m_settings.get_real("physical_parameters/nu");  // 1;
        real beta = m_settings.get_real("physical_parameters/beta");
        real kappa = m_settings.get_real("physical_parameters/kappa");
        real rhoa = m_settings.get_real("initial_conditions/rhoa");

        Functions::buoyancy_st_mms(m_field_controller->get_field_source_T(),
                                   0., nu, beta, kappa, g, rhoa);
    }
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters force source functions, once at initialisation
// ***************************************************************************************
void SolverController::force_source() {
    // Force
    if (m_settings.get("solver/source/force_fct") == SourceMethods::Buoyancy) {
        std::string dir = m_settings.get("solver/source/dir");
        if (!m_settings.get_bool("solver/source/use_init_values")) {
            real ambient_temperature_value = m_settings.get_real("solver/source/ambient_temperature_value");
            //m_field_controller->get_field_T().set_value(ambient_temperature_value);
            m_field_controller->get_field_T_ambient().set_value(ambient_temperature_value);
        }

        real beta = m_settings.get_real("physical_parameters/beta");
        real g = m_settings.get_real("physical_parameters/g");

        if (dir.find('x') != std::string::npos) {
            Functions::buoyancy_force(m_field_controller->get_field_force_x(),
                                      m_field_controller->get_field_T(),
                                      m_field_controller->get_field_T_ambient(),
                                      beta, g);
        }
        if (dir.find('y') != std::string::npos) {
            Functions::buoyancy_force(m_field_controller->get_field_force_y(),
                                      m_field_controller->get_field_T(),
                                      m_field_controller->get_field_T_ambient(),
                                      beta, g);
        }
        if (dir.find('z') != std::string::npos) {
            Functions::buoyancy_force(m_field_controller->get_field_force_z(),
                                      m_field_controller->get_field_T(),
                                      m_field_controller->get_field_T_ambient(),
                                      beta, g);
        }
    }
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters momentum source functions
// ***************************************************************************************
void SolverController::momentum_source() {
    // Momentum source
    std::string dir_vel = m_settings.get("solver/source/dir");
    if (dir_vel.find('x') != std::string::npos) {
        source_velocity->buoyancy_force(m_settings,
                                        m_field_controller->get_field_force_x(),
                                        m_field_controller->get_field_T(),
                                        m_field_controller->get_field_T_ambient());
    }
    if (dir_vel.find('y') != std::string::npos) {
        source_velocity->buoyancy_force(m_settings,
                                        m_field_controller->get_field_force_y(),
                                        m_field_controller->get_field_T(),
                                        m_field_controller->get_field_T_ambient());
    }
    if (dir_vel.find('z') != std::string::npos) {
        source_velocity->buoyancy_force(m_settings,
                                        m_field_controller->get_field_force_z(),
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
    // Momentum source
    if (m_has_momentum_source) {
        std::string forceFct = m_settings.get("solver/source/force_fct");
        if (forceFct == SourceMethods::Zero || \
            forceFct == SourceMethods::Uniform) {
        } else if (forceFct == SourceMethods::Buoyancy) {
#ifndef BENCHMARKING
            m_logger->info("Update f(T) ...");
#endif
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
