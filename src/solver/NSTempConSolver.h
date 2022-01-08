/// \file       NSTempConSolver.cpp 
/// \brief      Navier-Stokes Solver with force f(T)
/// \details    Defines the steps to solve advection, diffusion, pressure and add sources (dependent on T), solves temperature and concentration equation
/// \date       Sep 27, 2017
/// \author     Küsters
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_SOLVER_NSTEMPCONSOLVER_H_
#define ARTSS_SOLVER_NSTEMPCONSOLVER_H_

#include <string>
#include <vector>

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

class NSTempConSolver: public ISolver {
 public:
    NSTempConSolver(const Settings::solver_parameters &solver_settings, FieldController *field_controller);
    ~NSTempConSolver() override;

    void do_step(real t, bool sync) override;
    void update_source(real) override;
    void replace_heat_source(const Settings::solver::temperature_source &temperature_source, real t_cur) override;

private:
    const Settings::solver_parameters &m_solver_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
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

    FieldController *m_field_controller;
    ISourceFunction *m_source_function_concentration;
    ISourceFunction *m_source_function_temperature;

    bool m_add_source;
    bool m_add_temp_source;
    bool m_add_con_source;
    void control();
};

#endif /* ARTSS_SOLVER_NSTEMPCONSOLVER_H_ */
