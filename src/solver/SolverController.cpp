/// \file       SolverController.cpp
/// \brief      Manager class for everything regarding solver
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.

#include "SolverController.h"

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
#include "SolverSelection.h"
#include "../domain/DomainController.h"
#include "../domain/DomainData.h"
#include "../Functions.h"


SolverController::SolverController(Settings::Settings const &settings, const Settings::Settings_new &settings_new) :
        m_settings(settings), m_settings_new(settings_new) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = new FieldController();
    init_solver(m_settings_new.solver_parameters);
#ifndef BENCHMARKING
    m_logger->info("Start initialising....");
#endif
    set_up_fields(m_settings_new.solver_parameters.description);
#ifndef BENCHMARKING
    m_logger->debug("set up boundary");
#endif
    m_field_controller->set_up_boundary();
    update_sources(0, true);
    m_field_controller->update_data();
}

SolverController::~SolverController() {
    delete m_field_controller;
    delete m_solver;
}

void SolverController::init_solver(const Settings::solver_parameters &solver_settings) {
#ifndef BENCHMARKING
    m_logger->debug("initialise solver {}", solver_settings.description);
#endif
    if (solver_settings.description == SolverTypes::AdvectionSolver) {
        auto ic = std::get<Settings::initial_conditions::gauss_bubble>(m_settings_new.initial_conditions_parameters.ic.value());
        m_solver = new AdvectionSolver(solver_settings, m_field_controller, ic.velocity_lin);
    } else if (solver_settings.description == SolverTypes::AdvectionDiffusionSolver) {
        m_solver = new AdvectionDiffusionSolver(solver_settings, m_settings, m_field_controller);
    } else if (solver_settings.description == SolverTypes::DiffusionSolver) {
        m_solver = new DiffusionSolver(solver_settings, m_field_controller);
    } else if (solver_settings.description == SolverTypes::DiffusionTurbSolver) {
        m_solver = new DiffusionTurbSolver(solver_settings, m_settings, m_field_controller);
    } else if (solver_settings.description == SolverTypes::NSSolver) {
        m_solver = new NSSolver(solver_settings, m_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (solver_settings.description == SolverTypes::NSTurbSolver) {
        m_solver = new NSTurbSolver(solver_settings, m_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (solver_settings.description == SolverTypes::NSTempSolver) {
        m_solver = new NSTempSolver(solver_settings, m_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (solver_settings.description == SolverTypes::NSTempTurbSolver) {
        m_solver = new NSTempTurbSolver(solver_settings, m_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (solver_settings.description == SolverTypes::NSTempConSolver) {
        m_solver = new NSTempConSolver(solver_settings, m_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (solver_settings.description == SolverTypes::NSTempTurbConSolver) {
        m_solver = new NSTempTurbConSolver(solver_settings, m_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (solver_settings.description == SolverTypes::PressureSolver) {
        m_solver = new PressureSolver(solver_settings, m_field_controller);
    } else {
#ifndef BENCHMARKING
        m_logger->error("Solver {} not yet implemented! Simulation stopped!", solver_settings.description);
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
            CoordinateAxis axis = Mapping::match_axis(dir);
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
        auto domain_data = DomainData::getInstance();
        if (dir == "x") {
            real y1 = m_settings.get_real("initial_conditions/y1");
            real y2 = m_settings.get_real("initial_conditions/y2");
            real z1 = m_settings.get_real("initial_conditions/z1");
            real z2 = m_settings.get_real("initial_conditions/z2");
            size_t index_x1 = domain_data->get_index_x1();
            size_t index_x2 = domain_data->get_index_x2();
            size_t index_y1 = Utility::get_index(y1, domain_data->get_dy(), domain_data->get_Y1());
            size_t index_y2 = Utility::get_index(y2, domain_data->get_dy(), domain_data->get_Y1());
            size_t index_z1 = Utility::get_index(z1, domain_data->get_dz(), domain_data->get_Z1());
            size_t index_z2 = Utility::get_index(z2, domain_data->get_dz(), domain_data->get_Z1());
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
            size_t index_x1 = Utility::get_index(x1, domain_data->get_dx(), domain_data->get_X1());
            size_t index_x2 = Utility::get_index(x2, domain_data->get_dx(), domain_data->get_X1());
            size_t index_y1 = domain_data->get_index_y1();
            size_t index_y2 = domain_data->get_index_y2();
            size_t index_z1 = Utility::get_index(z1, domain_data->get_dz(), domain_data->get_Z1());
            size_t index_z2 = Utility::get_index(z2, domain_data->get_dz(), domain_data->get_Z1());
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
            size_t index_x1 = Utility::get_index(x1, domain_data->get_dx(), domain_data->get_X1());
            size_t index_x2 = Utility::get_index(x2, domain_data->get_dx(), domain_data->get_X1());
            size_t index_y1 = Utility::get_index(y1, domain_data->get_dy(), domain_data->get_Y1());
            size_t index_y2 = Utility::get_index(y2, domain_data->get_dy(), domain_data->get_Y1());
            size_t index_z1 = domain_data->get_index_z1();
            size_t index_z2 = domain_data->get_index_z2();
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
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    Field &sight = m_field_controller->get_field_sight();
    sight.update_host();
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
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
/// \brief  Updates time dependent parameters force source functions, once at initialisation
// ***************************************************************************************
void SolverController::force_source() {
    // Force
    if (m_settings.get("solver/source/force_fct") == SourceMethods::Buoyancy) {
        std::string dir = m_settings.get("solver/source/dir");
        if (!m_settings.get_bool("solver/source/use_init_values")) {
            real ambient_temperature_value = m_settings.get_real("solver/source/ambient_temperature_value");
            m_field_controller->get_field_T_ambient().set_value(ambient_temperature_value);
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
        ISource::buoyancy_force(m_settings,
                                m_field_controller->get_field_force_x(),
                                m_field_controller->get_field_T(),
                                m_field_controller->get_field_T_ambient());
    }
    if (dir_vel.find('y') != std::string::npos) {
        ISource::buoyancy_force(m_settings,
                                m_field_controller->get_field_force_y(),
                                m_field_controller->get_field_T(),
                                m_field_controller->get_field_T_ambient());
    }
    if (dir_vel.find('z') != std::string::npos) {
        ISource::buoyancy_force(m_settings,
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
        if (forceFct == SourceMethods::Buoyancy) {
#ifndef BENCHMARKING
            m_logger->info("Update f(T) ...");
#endif
            momentum_source();
        }
    }
    m_solver->update_source(t_cur);
    if (sync) {
#pragma acc wait
    }
}

void SolverController::solver_do_step(real t, bool sync) {
    m_solver->do_step(t, sync);
}
