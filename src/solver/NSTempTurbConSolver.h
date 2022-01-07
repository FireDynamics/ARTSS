/// \file       NSTempTurbConSolver.h
/// \brief      Defines the (fractional) steps to solve the incompressible Navier-Stokes equations with force f(T), turbulence and concentration C
/// \date       Oct 02, 2017
/// \author     Küsters
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_SOLVER_NSTEMPTURBCONSOLVER_H_
#define ARTSS_SOLVER_NSTEMPTURBCONSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/ISource.h"
#include "../interfaces/ITurbulence.h"
#include "../field/FieldController.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"
#include "../interfaces/ISourceFunction.h"

class NSTempTurbConSolver : public ISolver {
 public:
    NSTempTurbConSolver(const Settings::solver_parameters &solver_settings, FieldController *field_controller);
    ~NSTempTurbConSolver() override;

    void do_step(real t, bool sync) override;

    void update_source(real) override;
    void replace_heat_source(const Settings::solver::temperature_source &temperature_source) override;
private:
    const Settings::solver_parameters &m_solver_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    FieldController *m_field_controller;

    IAdvection *adv_vel;
    IDiffusion *dif_vel;
    IAdvection *adv_temp;
    IDiffusion *dif_temp;
    IAdvection *adv_con;
    IDiffusion *dif_con;
    IPressure *pres;
    ISource *sou_vel;
    ISource *sou_temp;
    ISource *sou_con;
    ITurbulence *mu_tub;
    ISourceFunction *m_source_function_concentration;
    ISourceFunction *m_source_function_temperature;

    void control();

    bool m_add_source;
    bool m_add_temp_source;
    bool m_add_con_source;
};

#endif /* ARTSS_SOLVER_NSTEMPTURBCONSOLVER_H_ */
