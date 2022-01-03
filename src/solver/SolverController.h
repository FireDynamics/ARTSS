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
    explicit SolverController(const Settings::Settings &settings);
    ~SolverController();

    void solver_do_step(real t, bool sync);
    void update_sources(real t_cur, bool sync);

    FieldController *get_field_controller() { return m_field_controller; }

 private:
    void init_solver(const Settings::solver_parameters &solver_settings);
    void set_up_fields(const std::string &solver_description, const Settings::initial_conditions_parameters &ic_settings);

    void force_source();
    void momentum_source();

    const Settings::Settings &m_settings;

    FieldController *m_field_controller;
    ISolver *m_solver;

    bool m_has_momentum_source = false;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

};

#endif /* ARTSS_SOLVER_SOLVERCONTROLLER_H_ */

