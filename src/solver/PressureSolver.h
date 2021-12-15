/// \file       PressureSolver.h
/// \brief      Defines the steps to solve the pressure Poisson equation
/// \date       Sep 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_PRESSURESOLVER_H_
#define ARTSS_SOLVER_PRESSURESOLVER_H_

#include "SolverSelection.h"
#include "../pressure/VCycleMG.h"
#include "../boundary/DomainData.h"
#include "../boundary/BoundaryData.h"
#include "../interfaces/ISolver.h"
#include "../interfaces/IPressure.h"
#include "../field/FieldController.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"

class PressureSolver : public ISolver {
public:
    PressureSolver(Settings::Settings const &settings, FieldController *field_controller);
    ~PressureSolver();
    void do_step(real t, bool sync) override;

private:
    Settings::Settings const &m_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    FieldController *m_field_controller;
    IPressure *pres;
    void control();
};

#endif /* ARTSS_SOLVER_PRESSURESOLVER_H_ */
