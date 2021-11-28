/// \file       NSTempConSolver.cpp 
/// \brief      Navier-Stokes Solver with force f(T)
/// \details    Defines the steps to solve advection, diffusion, pressure and add sources (dependent on T), solves temperature and concentration equation
/// \date       Sep 27, 2017
/// \author     Küsters
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_SOLVER_NSTEMPCONSOLVER_H_
#define ARTSS_SOLVER_NSTEMPCONSOLVER_H_


#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/ISource.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../field/FieldController.h"
#include "../utility/Utility.h"

class NSTempConSolver: public ISolver {
 public:
    NSTempConSolver(Settings const &settings, FieldController *field_controller);
    ~NSTempConSolver() override;

    void do_step(real t, bool sync) override;

private:
    Settings const &m_settings;
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

    std::string m_dir_vel;

    void control();

    std::string m_forceFct;
    bool m_hasDissipation;
    std::string m_tempFct;
    std::string m_conFct;
};

#endif /* ARTSS_SOLVER_NSTEMPCONSOLVER_H_ */
