/// \file       TimeIntegration.h
/// \brief      Runs the time loop
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_TIMEINTEGRATION_H_
#define ARTSS_TIMEINTEGRATION_H_

#include "interfaces/ISolver.h"
#include "interfaces/ISource.h"
#include "utility/GlobalMacrosTypes.h"
#include "solver/SolverController.h"
#include "adaption/Adaption.h"

#ifndef BENCHMARKING
#include "analysis/Analysis.h"
#include "analysis/Solution.h"
#include "visualisation/Visual.h"
#else
// only needed if no logger will be available
#include <iostream>
#endif

class TimeIntegration {
 public:
    explicit TimeIntegration(SolverController *sc);
    void run();

 private:
    real m_dt;
    real m_t_end;
    real m_t_cur;

    FieldController *m_field_controller;
    SolverController *m_solver_controller;
    Adaption *m_adaption;
#ifndef BENCHMARKING
    Visual *m_visual;
    Solution *m_solution;
    Analysis *m_analysis;
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_TIMEINTEGRATION_H_ */
