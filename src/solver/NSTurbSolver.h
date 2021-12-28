/// \file       NSTurbSolver.h
/// \brief      Defines the steps to solve advection, diffusion, pressure and add sources (with LES turbulence)
/// \date       Feb 15, 2017
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_NSTURBSOLVER_H_
#define ARTSS_SOLVER_NSTURBSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/ISource.h"
#include "../interfaces/ITurbulence.h"
#include "../field/FieldController.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"

class NSTurbSolver : public ISolver {
 public:
    NSTurbSolver(const Settings::solver_parameters &solver_settings, Settings::Settings const &setings, FieldController *field_controller);
    ~NSTurbSolver();

    void do_step(real t, bool sync) override;
    void update_source(real) override {};

 private:
    Settings::Settings const &m_settings;
    const Settings::solver_parameters &m_solver_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    IAdvection *adv_vel;
    IDiffusion *dif_vel;
    IPressure *pres;
    ISource *sou_vel;
    ITurbulence *mu_tub;

    FieldController *m_field_controller;

    void control();

    std::string m_force_function;
};

#endif /* ARTSS_SOLVER_NSTURBSOLVER_H_ */
