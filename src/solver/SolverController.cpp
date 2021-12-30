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


SolverController::SolverController(const Settings::Settings_new &settings_new) :
        m_settings_new(settings_new) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_field_controller = new FieldController();
    init_solver(m_settings_new.solver_parameters);
#ifndef BENCHMARKING
    m_logger->info("Start initialising....");
#endif
    set_up_fields(m_settings_new.solver_parameters.description, m_settings_new.initial_conditions_parameters);
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
        m_solver = new AdvectionDiffusionSolver(solver_settings, m_field_controller);
    } else if (solver_settings.description == SolverTypes::DiffusionSolver) {
        m_solver = new DiffusionSolver(solver_settings, m_field_controller);
    } else if (solver_settings.description == SolverTypes::DiffusionTurbSolver) {
        m_solver = new DiffusionTurbSolver(solver_settings, m_field_controller);
    } else if (solver_settings.description == SolverTypes::NSSolver) {
        m_solver = new NSSolver(solver_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (solver_settings.description == SolverTypes::NSTurbSolver) {
        m_solver = new NSTurbSolver(solver_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (solver_settings.description == SolverTypes::NSTempSolver) {
        m_solver = new NSTempSolver(solver_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (solver_settings.description == SolverTypes::NSTempTurbSolver) {
        m_solver = new NSTempTurbSolver(solver_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (solver_settings.description == SolverTypes::NSTempConSolver) {
        m_solver = new NSTempConSolver(solver_settings, m_field_controller);
        m_has_momentum_source = true;
    } else if (solver_settings.description == SolverTypes::NSTempTurbConSolver) {
        m_solver = new NSTempTurbConSolver(solver_settings, m_field_controller);
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

void SolverController::set_up_fields(const std::string &string_solver, const Settings::initial_conditions_parameters &ic_settings) {
    std::string string_init_usr_fct = ic_settings.usr_fct;
    if (string_init_usr_fct == FunctionNames::gauss_bubble) {
        auto gauss_bubble = std::get<Settings::initial_conditions::gauss_bubble>(ic_settings.ic.value());
        Functions::gauss_bubble(m_field_controller->get_field_u(), 0., gauss_bubble);
        Functions::gauss_bubble(m_field_controller->get_field_v(), 0., gauss_bubble);
        Functions::gauss_bubble(m_field_controller->get_field_w(), 0., gauss_bubble);
    } else if (string_init_usr_fct == FunctionNames::drift) {
        auto drift = std::get<Settings::initial_conditions::drift>(ic_settings.ic.value());
        Functions::drift(m_field_controller->get_field_u(),
                         m_field_controller->get_field_v(),
                         m_field_controller->get_field_w(),
                         m_field_controller->get_field_p(),
                         drift);
    } else if (string_init_usr_fct == FunctionNames::exp_sinus_prod) {
        // Diffusion test case
        auto exp_sinus_prod = std::get<Settings::initial_conditions::exp_sinus_prod>(ic_settings.ic.value());
        Functions::exp_sinus_prod(m_field_controller->get_field_u(), 0., exp_sinus_prod);
        Functions::exp_sinus_prod(m_field_controller->get_field_v(), 0., exp_sinus_prod);
        Functions::exp_sinus_prod(m_field_controller->get_field_w(), 0., exp_sinus_prod);
    } else if (string_init_usr_fct == FunctionNames::hat) {
        auto hat = std::get<Settings::initial_conditions::hat>(ic_settings.ic.value());
        Functions::hat(m_field_controller->get_field_u(), hat);
        Functions::hat(m_field_controller->get_field_v(), hat);
        Functions::hat(m_field_controller->get_field_w(), hat);
    } else if (string_init_usr_fct == FunctionNames::exp_sinus_sum) {
        // Burgers (=nonlinear Advection + Diffusion) test case
        Functions::exp_sinus_sum(m_field_controller->get_field_u(),
                                 m_field_controller->get_field_v(),
                                 m_field_controller->get_field_w(),
                                 0.);
    } else if (string_init_usr_fct == FunctionNames::sin_sin_sin) {
        // Pressure test case
        auto sin_sin_sin = std::get<Settings::initial_conditions::sin_sin_sin>(ic_settings.ic.value());
        m_field_controller->get_field_p().set_value(0.);
        Functions::sin_sin_sin(m_field_controller->get_field_rhs(), sin_sin_sin);
    } else if (string_init_usr_fct == FunctionNames::mcdermott) {
        // NavierStokes test case: McDermott (no force, no temperature) 2D
        auto mc_dermott = std::get<Settings::initial_conditions::mc_dermott>(ic_settings.ic.value());
        Functions::mcdermott(m_field_controller->get_field_u(),
                             m_field_controller->get_field_v(),
                             m_field_controller->get_field_w(),
                             m_field_controller->get_field_p(),
                             0., mc_dermott);
        m_field_controller->get_field_p().set_value(0.);
    } else if (string_init_usr_fct == FunctionNames::vortex) {
        // NavierStokes test case: Vortex (no force, no temperature) 2D
        auto vortex = std::get<Settings::initial_conditions::vortex>(ic_settings.ic.value());
        Functions::vortex(m_field_controller->get_field_u(),
                          m_field_controller->get_field_v(),
                          m_field_controller->get_field_w(),
                          m_field_controller->get_field_p(),
                          vortex);
        m_field_controller->get_field_p().set_value(0.);
    } else if (string_init_usr_fct == FunctionNames::vortex_y) {
        // NavierStokes test case: Vortex (no force, no temperature) 2D
        auto vortex = std::get<Settings::initial_conditions::vortex>(ic_settings.ic.value());
        Functions::vortex_y(m_field_controller->get_field_u(),
                            m_field_controller->get_field_v(),
                            m_field_controller->get_field_w(),
                            m_field_controller->get_field_p(),
                            vortex);
    } else if (string_init_usr_fct == FunctionNames::beltrami) {
        // NavierStokes test case: Beltrami  (no force, no temperature) 3D
        auto beltrami = std::get<Settings::initial_conditions::beltrami>(ic_settings.ic.value());

        Functions::beltrami(m_field_controller->get_field_u(),
                            m_field_controller->get_field_v(),
                            m_field_controller->get_field_w(),
                            m_field_controller->get_field_p(),
                            0., beltrami);
        // m_field_controller->field_p->set_value(0.);
    } else if (string_init_usr_fct == FunctionNames::uniform) {
        // Uniform Temperature unequal to zero
        auto uniform = std::get<Settings::initial_conditions::uniform>(ic_settings.ic.value());
        Functions::uniform(m_field_controller->get_field_T(), uniform);
        if (ic_settings.random) {
            Functions::random(m_field_controller->get_field_T(), ic_settings.random_parameters.value());
        }
    } else if (string_init_usr_fct == "LayersT") {
        auto layers = std::get<Settings::initial_conditions::layers_temperature>(ic_settings.ic.value());
        Functions::layers(m_field_controller->get_field_T(), layers);

        if (ic_settings.random) {
            Functions::random(m_field_controller->get_field_T(), ic_settings.random_parameters.value());
        }
    } else if (string_init_usr_fct == FunctionNames::zero) {
#ifndef BENCHMARKING
            m_logger->info("Initial values all set to zero!");
#endif
    } else if (string_init_usr_fct == FunctionNames::jet) {
        auto jet = std::get<Settings::initial_conditions::jet>(ic_settings.ic.value());
        Functions::jet(m_field_controller->get_field_u(), jet);
        if (ic_settings.random) {
            Functions::random(m_field_controller->get_field_u(), ic_settings.random_parameters.value());
        }
    } else {
#ifndef BENCHMARKING
        m_logger->warn("Unknown initial function {}", string_init_usr_fct);
        m_logger->info("Initial values all set to zero!");
#endif
    }
    if (string_solver == SolverTypes::NSTempSolver ||
        string_solver == SolverTypes::NSTempConSolver ||
        string_solver == SolverTypes::NSTempTurbConSolver ||
        string_solver == SolverTypes::NSTempTurbSolver) {
        force_source();
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

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters force source functions, once at initialisation
// ***************************************************************************************
void SolverController::force_source() {
    // Force
    std::string force_fct = m_settings_new.solver_parameters.source.force_fct;
    if (force_fct == SourceMethods::Buoyancy) {
        auto buoyancy = std::get<Settings::solver::source_solvers::buoyancy>(m_settings_new.solver_parameters.source.force_function);
        if (!buoyancy.use_init_values) {
            m_field_controller->get_field_T_ambient().set_value(buoyancy.ambient_temperature_value.value());
        }
    } else if (force_fct == SourceMethods::Uniform) {
        auto uniform = std::get<Settings::solver::source_solvers::uniform>(m_settings_new.solver_parameters.source.force_function);
        std::vector<CoordinateAxis> dir = m_settings_new.solver_parameters.source.direction;
        bool force_x = std::find(dir.begin(), dir.end(), CoordinateAxis::X) != dir.end();
        if (force_x) {
            Settings::initial_conditions::uniform uni{uniform.velocity_value[CoordinateAxis::X]};
            Functions::uniform(m_field_controller->get_field_force_x(), uni);
        }
        bool force_y = std::find(dir.begin(), dir.end(), CoordinateAxis::Y) != dir.end();
        if (force_y) {
            Settings::initial_conditions::uniform uni{uniform.velocity_value[CoordinateAxis::Y]};
            Functions::uniform(m_field_controller->get_field_force_y(), uni);
        }
        bool force_z = std::find(dir.begin(), dir.end(), CoordinateAxis::Z) != dir.end();
        if (force_z) {
            Settings::initial_conditions::uniform uni{uniform.velocity_value[CoordinateAxis::Z]};
            Functions::uniform(m_field_controller->get_field_force_z(), uni);
        }
    }

}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates time dependent parameters momentum source functions
// ***************************************************************************************
void SolverController::momentum_source() {
    // Momentum source
    std::vector<CoordinateAxis> dir = m_settings_new.solver_parameters.source.direction;
    bool force_x = std::find(dir.begin(), dir.end(), CoordinateAxis::X) != dir.end();
    if (force_x) {
        ISource::buoyancy_force(m_field_controller->get_field_force_x(),
                                m_field_controller->get_field_T(),
                                m_field_controller->get_field_T_ambient());
    }
    bool force_y = std::find(dir.begin(), dir.end(), CoordinateAxis::Y) != dir.end();
    if (force_y) {
        ISource::buoyancy_force(m_field_controller->get_field_force_y(),
                                m_field_controller->get_field_T(),
                                m_field_controller->get_field_T_ambient());
    }
    bool force_z = std::find(dir.begin(), dir.end(), CoordinateAxis::Z) != dir.end();
    if (force_z) {
        ISource::buoyancy_force(m_field_controller->get_field_force_z(),
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
        if (m_settings_new.solver_parameters.source.force_fct == SourceMethods::Buoyancy) {
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
