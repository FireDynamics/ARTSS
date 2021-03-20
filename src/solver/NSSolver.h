/// \file       NSTempTurbConSolver.h
/// \brief      Defines the (fractional) steps to solve the incompressible
///             Navier-Stokes equations with force f(T), turbulence and
///             concentration C
/// \date       Oct 02, 2017
/// \author     Küsters
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_SOLVER_NSSOLVER2_H_
#define ARTSS_SOLVER_NSSOLVER2_H_

#include <memory>
#include <string>

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/ISource.h"
#include "../interfaces/ITurbulence.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"
#include "../field/FieldController.h"

template <bool t_temp, bool t_turb, bool t_con>
class NSSolver2 : public ISolver {
 public:
    explicit NSSolver2(FieldController *field_controller);
    ~NSSolver2();

    void do_step(real t, bool sync) override;

 private:
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

    real m_nu;
    real m_kappa;
    real m_gamma;
    std::string m_dir_vel;

    static void control();

    bool m_has_turbulence_temperature = false;
    bool m_has_turbulence_concentration = false;
    bool m_has_dissipation = false;
    std::string m_forceFct;
    std::string m_tempFct;
    std::string m_conFct;
};

#endif /* ARTSS_SOLVER_NSSOLVER2_H_ */
