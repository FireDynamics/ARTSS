/// \file       DiffusionTurbSolver.h
/// \brief      Defines the steps to solve the turbulent diffusion equation
/// \date       Aug 18, 2016
/// \author     Suryanarayana Maddu
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_DIFFUSIONTURBSOLVER_H_
#define ARTSS_SOLVER_DIFFUSIONTURBSOLVER_H_

#include "SolverSelection.h"
#include "../domain/DomainData.h"
#include "../field/FieldController.h"
#include "../interfaces/ISolver.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/ITurbulence.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"
#include "../utility/Utility.h"

class DiffusionTurbSolver: public ISolver {
 public:
    explicit DiffusionTurbSolver(Settings::Settings const &settings, FieldController *field_controller);
    ~DiffusionTurbSolver();

    void do_step(real t, bool sync) override;

 private:
    Settings::Settings const &m_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    FieldController *m_field_controller;

    IDiffusion *dif;
    ITurbulence *mu_tub;

    void control();
};

#endif /* ARTSS_SOLVER_DIFFUSIONTURBSOLVER_H_ */
