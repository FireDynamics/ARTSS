/// \file       SolverController.h
/// \brief      Manager class for everything regarding solver
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_SOLVER_SOLVERCONTROLLER_H
#define ARTSS_SOLVER_SOLVERCONTROLLER_H

#include "../interfaces/ISolver.h"
#include "../interfaces/ISourceFunction.h"
#include "../field/FieldController.h"
#include "../utility/Utility.h"

class SolverController {
public:
    SolverController();
    ~SolverController();

    void solver_do_step(real t, bool sync);
    void update_sources(real t_cur, bool sync);

    FieldController* get_field_controller() { return m_field_controller; };

private:
    void set_up_sources(const std::string &string_solver);
    void init_solver(const std::string& string_solver);
    void set_up_fields(const std::string& string_solver);
    void call_random(Field &field);

    void force_source();
    void temperature_source();
    void momentum_source();

    FieldController *m_field_controller;
    ISolver *m_solver;
    ISourceFunction *m_source_function_concentration;
    ISourceFunction *m_source_function_temperature;

    ISource *source_temperature;
    ISource *source_velocity;
    ISource *source_concentration;

    bool m_has_temperature = false;
    bool m_has_momentum_source = false;
    bool m_has_turbulence = false;
    bool m_has_concentration = false;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};


#endif /* ARTSS_SOLVER_SOLVERCONTROLLER_H */
