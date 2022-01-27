/// \file       NSSolver.h
/// \brief      Defines the (fractional) steps to solve the incompressible Navier-Stokes equations
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved
#ifndef ARTSS_SOLVER_NSSOLVER_H_
#define ARTSS_SOLVER_NSSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/ISource.h"
#include "../field/FieldController.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"

#ifdef BENCHMARKING
#include "../utility/Utility.h"
#endif


class NSSolver : public ISolver {
 public:
    NSSolver(const Settings::solver_parameters &solver_settings, FieldController *field_controller);
    ~NSSolver() override;

    void do_step(real t, bool sync) override;
    void update_source(real) override {};
    void replace_heat_source(const Settings::solver::temperature_source &temperature_source) override {};

 private:
    const Settings::solver_parameters &m_solver_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    FieldController *m_field_controller;

    IAdvection *adv_vel;
    IDiffusion *dif_vel;
    IPressure *pres;
    ISource *sou_vel;

    bool m_add_source;

    void control();
};

#endif /* ARTSS_SOLVER_NSSOLVER_H_ */
