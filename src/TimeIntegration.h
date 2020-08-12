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

#ifndef BENCHMARKING

#include "analysis/Analysis.h"
#include "analysis/Solution.h"
#include "visualisation/Visual.h"
#include "adaption/Adaption.h"

#endif

class TimeIntegration {
public:
    explicit TimeIntegration(SolverController *sc);
    void run();

private:
    void set_up();

    real m_dt;
    real m_t_end;
    real m_t_cur;

    FieldController *m_field_controller;
    SolverController *m_solver_controller;
    Adaption *m_adaption;
    Visual *m_visual;
    Solution *m_solution;
    Analysis *m_analysis;
};

#endif /* ARTSS_TIMEINTEGRATION_H_ */
