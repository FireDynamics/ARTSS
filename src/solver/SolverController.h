/// \file       SolverController.h
/// \brief      Manager class for everything regarding solver
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.

#ifndef ARTSS_SOLVER_SOLVERCONTROLLER_H_
#define ARTSS_SOLVER_SOLVERCONTROLLER_H_

#include <string>

#include "../interfaces/ISolver.h"
#include "../interfaces/ISource.h"
#include "../interfaces/ISourceFunction.h"
#include "../field/FieldController.h"
#include "../utility/Utility.h"
#include "../utility/settings/Settings.h"


class SolverController {
 public:
    explicit SolverController(Settings::Settings const &settings, const Settings::Settings_new &settings_new);
    ~SolverController();

    void solver_do_step(real t, bool sync);
    void update_sources(real t_cur, bool sync);

    FieldController *get_field_controller() { return m_field_controller; }

 private:
    void set_up_source(const Settings::solver::source_solver &source_settings);
    void init_solver(const std::string& string_solver);
    void set_up_fields(const std::string& string_solver);
    void call_random(Field &field);

    void force_source();
    void momentum_source();

    Settings::Settings const &m_settings;
    const Settings::Settings_new &m_settings_new;

    FieldController *m_field_controller;
    ISolver *m_solver;

    ISource *source_temperature;
    ISource *source_velocity;
    ISource *source_concentration;

    bool m_has_momentum_source = false;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_SOLVER_SOLVERCONTROLLER_H_ */

