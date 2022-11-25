/// \file       DiffusionSolver.h
/// \brief      Defines the steps to solve the diffusion equation
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_DIFFUSIONSOLVER_H_
#define ARTSS_SOLVER_DIFFUSIONSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IDiffusion.h"
#include "../field/FieldController.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"

class DiffusionSolver: public ISolver {
 public:
    DiffusionSolver(const Settings::solver_parameters &solver_settings, FieldController *field_controller);
    ~DiffusionSolver() override;

    void do_step(real t, bool sync) override;
    void update_source(real) override {};
    void update_obstacle_change() override {};
    void replace_heat_source(const Settings::solver::temperature_source &temperature_source) override {};

 private:
    const Settings::solver_parameters &m_solver_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    FieldController *m_field_controller;
    IDiffusion *dif;

    void control();
};

#endif /* ARTSS_SOLVER_DIFFUSIONSOLVER_H_ */
