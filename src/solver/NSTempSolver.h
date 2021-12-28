/// \file       NSTempSolver.h
/// \brief      Defines the steps to solve advection, diffusion, pressure and add sources (dependent on T), solves temperature equation
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_NSTEMPSOLVER_H_
#define ARTSS_SOLVER_NSTEMPSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/ISource.h"
#include "../field/FieldController.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"
#include "../interfaces/ISourceFunction.h"

class NSTempSolver : public ISolver {
public:
    NSTempSolver(const Settings::solver_parameters &solver_settings, Settings::Settings const &settings, FieldController *field_controller);
    ~NSTempSolver();

    void do_step(real t, bool sync) override;
    void update_source(real) override;
 private:
    Settings::Settings const &m_settings;
    const Settings::solver_parameters &m_solver_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    FieldController *m_field_controller;

    IAdvection *adv_vel;
    IDiffusion *dif_vel;
    IAdvection *adv_temp;
    IDiffusion *dif_temp;
    IPressure *pres;
    ISource *sou_vel;
    ISource *sou_temp;

    ISourceFunction *m_source_function_temperature;
    std::string m_dir_vel;

    void control();

    std::string m_forceFct;
    bool m_has_dissipation;
    std::string m_tempFct;
};

#endif /* ARTSS_SOLVER_NSTEMPSOLVER_H_ */
